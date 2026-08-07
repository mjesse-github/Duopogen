#!/usr/bin/env python3
"""
glimpse.py -- GLIMPSE2 imputation/phasing backend for duopogen.

An alternative to the Beagle 4.1 path in Duopogen.py. Same inputs (per-
chromosome filtered BAM lists from preProcess, a 1KG3 phased panel), same
kind of output (phased per-chromosome VCFs), different engine.

Why it exists: Beagle 4.1 rebuilds its HMM from the 3202-sample panel for
every window of every run, and the window count grows with donor count.
GLIMPSE2 binarises the panel ONCE per chunk and reuses it for every donor
and every run, and conditions on a PBWT-selected subset of reference
haplotypes rather than all of them.

Three deliberate differences from the Beagle path, all accuracy-relevant:

  1. GLs are computed at ALL panel sites via `bcftools call -C alleles`,
     including sites where every read matches the reference. The Beagle
     path drops those (the `grep -v '<*>'`), discarding real hom-ref
     evidence. This is likely part of why its het/hom-alt ratio runs low.
  2. A genetic map is used. The Beagle runs log "No genetic map is
     specified: using 1 cM = 1 Mb", which is a poor approximation.
  3. Output covers every panel site with an INFO score, not only observed
     ones -- Beagle runs with impute=false.

Command syntax follows the GLIMPSE2 tutorial:
https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/
"""

import os
import shutil
import logging
from multiprocessing import Pool

from germline import CHROMS, runCMD, print_parameters_given

logger = logging.getLogger(__name__)

GLIMPSE_BINS = ("GLIMPSE2_chunk", "GLIMPSE2_split_reference",
                "GLIMPSE2_phase", "GLIMPSE2_ligate")


def _which(prog):
	p = shutil.which(prog)
	assert p is not None, \
		"{} not found on PATH -- is duopogen_env active and GLIMPSE2 installed?".format(prog)
	return p


def check_dependencies_glimpse(args):
	for prog in ("bcftools", "bgzip", "tabix") + GLIMPSE_BINS:
		logger.debug("Using {} -> {}".format(prog, _which(prog)))
	assert os.path.isdir(args.map_dir), \
		"Genetic map directory {} not found. GLIMPSE2 ships these under " \
		"maps/genetic_maps.b38/".format(args.map_dir)


def validate_user_setting_glimpse(args):
	assert os.path.isfile(args.reference), \
		"Reference FASTA {} not found".format(args.reference)
	assert os.path.isfile(args.reference + ".fai"), \
		"{}.fai not found -- run `samtools faidx {}`".format(
			args.reference, args.reference)
	assert os.path.isdir(args.imputation_panel), \
		"Panel directory {} not found".format(args.imputation_panel)

	for chrom in _chroms(args):
		gmap = os.path.join(args.map_dir, "{}.b38.gmap.gz".format(chrom))
		assert os.path.isfile(gmap), "Genetic map not found: {}".format(gmap)
		panel = _panel_vcf(args, chrom)
		assert os.path.isfile(panel), "Panel file not found: {}".format(panel)
		if args.step in ("varScan", "phase", "ligate", "all"):
			bamlst = args.out + "/Bam/" + chrom + ".filter.bam.lst"
			assert os.path.isfile(bamlst), \
				"{} not found -- did preProcess complete?".format(bamlst)


def _chroms(args):
	"""Which chromosomes to process. --chroms overrides the default 1-22."""
	if getattr(args, "chroms", None):
		return [c.strip() for c in args.chroms.split(",") if c.strip()]
	return CHROMS


def _panel_vcf(args, chrom):
	return (args.imputation_panel
		+ "CCDG_14151_B01_GRM_WGS_2020-08-05_" + chrom
		+ ".filtered.shapeit2-duohmm-phased.vcf.gz")


def _paths(args, chrom):
	out = os.path.abspath(args.out)
	pdir = os.path.join(os.path.abspath(args.panel_dir), chrom)
	return dict(
		out=out,
		pdir=pdir,
		bcf=os.path.join(pdir, "panel.{}.bcf".format(chrom)),
		sites_vcf=os.path.join(pdir, "sites.{}.vcf.gz".format(chrom)),
		sites_tsv=os.path.join(pdir, "sites.{}.tsv.gz".format(chrom)),
		chunks=os.path.join(pdir, "chunks.{}.txt".format(chrom)),
		binpfx=os.path.join(pdir, "split", "panel"),
		gmap=os.path.join(os.path.abspath(args.map_dir),
			"{}.b38.gmap.gz".format(chrom)),
		bamlst=os.path.join(out, "Bam", "{}.filter.bam.lst".format(chrom)),
		gl=os.path.join(out, "glimpse", "{}.gl.vcf.gz".format(chrom)),
		impdir=os.path.join(out, "glimpse", "impute"),
		phased=os.path.join(out, "glimpse", "{}.phased.bcf".format(chrom)),
		script=os.path.join(out, "Script", "glimpse_{}".format(chrom)),
	)


def _read_chunks(path):
	"""GLIMPSE2_chunk output: ID, CHR, input-region (with buffer),
	output-region (without). Returns a list of (id, chr, irg, org)."""
	rows = []
	with open(path) as f:
		for line in f:
			fields = line.split()
			if len(fields) >= 4:
				rows.append((fields[0], fields[1], fields[2], fields[3]))
	assert rows, "No chunks parsed from {}".format(path)
	return rows


# ---------------------------------------------------------------- steps

def _prep_panel_cmd(args, chrom):
	"""QC the panel to biallelic SNPs, extract sites, chunk, binarise.
	Runs once per panel; the .bin files are reused by every later run."""
	p = _paths(args, chrom)
	os.makedirs(os.path.join(p["pdir"], "split"), exist_ok=True)

	excl = ""
	if args.exclude_samples:
		# For benchmarking against a sample that is IN the panel (e.g.
		# NA12878 and her parents), you MUST drop them or concordance is
		# meaningless. See the GLIMPSE2 tutorial, section 2.2.
		excl = " -s ^" + args.exclude_samples

	cmd = "set -euo pipefail\n"
	cmd += ("bcftools norm -m -any {panel} -Ou --threads 2"
		" | bcftools view -m 2 -M 2 -v snps{excl} --threads 2 -Ob -o {bcf}\n"
		"bcftools index -f {bcf} --threads 2\n").format(
			panel=_panel_vcf(args, chrom), excl=excl, bcf=p["bcf"])

	cmd += ("bcftools view -G -Oz -o {sv} {bcf}\n"
		"bcftools index -f {sv}\n").format(sv=p["sites_vcf"], bcf=p["bcf"])

	# The TSV drives `bcftools call -C alleles`: it pins REF,ALT per site so
	# a PL is emitted even where no non-reference read was seen.
	cmd += ("bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {sv}"
		" | bgzip -c > {st}\n"
		"tabix -s1 -b2 -e2 -f {st}\n").format(sv=p["sites_vcf"], st=p["sites_tsv"])

	cmd += ("GLIMPSE2_chunk --input {sv} --region {chrom}"
		" --map {gmap} --sequential --output {chunks} --threads {t}\n").format(
			sv=p["sites_vcf"], chrom=chrom, gmap=p["gmap"],
			chunks=p["chunks"], t=args.bin_threads)

	# One .bin per chunk, named ${prefix}_${CHR}_${REGS}_${REGE}.bin
	cmd += ('while IFS="" read -r LINE || [ -n "$LINE" ]; do\n'
		'  IRG=$(echo $LINE | cut -d" " -f3)\n'
		'  ORG=$(echo $LINE | cut -d" " -f4)\n'
		'  GLIMPSE2_split_reference --reference {bcf} --map {gmap}'
		' --input-region ${{IRG}} --output-region ${{ORG}}'
		' --output {binpfx} --threads {t}\n'
		'done < {chunks}\n').format(bcf=p["bcf"], gmap=p["gmap"],
			binpfx=p["binpfx"], chunks=p["chunks"], t=args.bin_threads)
	return cmd


def _varscan_cmd(args, chrom):
	"""Genotype likelihoods at every panel site, hom-ref included."""
	p = _paths(args, chrom)
	os.makedirs(os.path.join(p["out"], "glimpse"), exist_ok=True)

	# -I skips indels (bcftools GLs for indels are not reliable, and the
	#    panel is filtered to SNPs above anyway)
	# -E recomputes BAQ; drop it if you A/B and find it costs more than it
	#    buys on ATAC data
	# --ns 0 is inherited from the Beagle path: it disables bcftools'
	#    default skipping of UNMAP,SECONDARY,QCFAIL,DUP. Keeping it makes
	#    this comparable to the Beagle run, but including PCR duplicates
	#    makes GLs overconfident -- worth revisiting once the comparison
	#    is done.
	cmd = "set -euo pipefail\n"
	cmd += ("( bcftools mpileup -f {ref} -I -E -a 'FORMAT/DP'"
		" -T {sv} -r {chrom} -b {bamlst}"
		" -q {q} -Q {Q} --ns 0 -Ou"
		" | bcftools call -Aim -C alleles -T {st} -Oz -o {gl}"
		" ) 2> {log}\n"
		"bcftools index -f {gl}\n").format(
			ref=args.reference, sv=p["sites_vcf"], chrom=chrom,
			bamlst=p["bamlst"], q=args.min_mq, Q=args.min_bq,
			st=p["sites_tsv"], gl=p["gl"],
			log=p["gl"].replace(".vcf.gz", ".log"))

	cmd += ('n=$(bcftools view -H {gl} | wc -l)\n'
		'[ "$n" -gt 0 ] || {{ echo "ERROR: {chrom} produced 0 GL records" >&2;'
		' exit 1; }}\n').format(gl=p["gl"], chrom=chrom)
	return cmd


def _phase_cmds(args, chrom):
	"""One job per chunk. Returns a list of shell commands."""
	p = _paths(args, chrom)
	os.makedirs(p["impdir"], exist_ok=True)
	assert os.path.isfile(p["chunks"]), \
		"{} not found -- run -s prepPanel first".format(p["chunks"])

	cmds = []
	for _, chrm, irg, org in _read_chunks(p["chunks"]):
		regs, rege = irg.split(":")[1].split("-")
		binfile = "{}_{}_{}_{}.bin".format(p["binpfx"], chrm, regs, rege)
		outbcf = os.path.join(p["impdir"],
			"imputed_{}_{}_{}.bcf".format(chrm, regs, rege))
		cmds.append("set -euo pipefail\n"
			+ "GLIMPSE2_phase --input-gl {gl} --reference {bin}"
			  " --output {out} --threads {t}\n".format(
				gl=p["gl"], bin=binfile, out=outbcf, t=args.phase_threads))
	return cmds


def _ligate_cmd(args, chrom):
	p = _paths(args, chrom)
	lst = os.path.join(p["impdir"], "list.{}.txt".format(chrom))
	# -1v sorts numerically so chunks ligate in genomic order.
	return ("set -euo pipefail\n"
		"ls -1v {impdir}/imputed_{chrom}_*.bcf > {lst}\n"
		"GLIMPSE2_ligate --input {lst} --output {phased} --threads {t}\n"
		"bcftools index -f {phased}\n").format(
			impdir=p["impdir"], chrom=chrom, lst=lst,
			phased=p["phased"], t=args.bin_threads)


def _run(cmds, nthreads, label):
	"""Write each command to its own script and run the batch through a Pool."""
	if not cmds:
		return
	logger.info("%s: %d job(s) on %d worker(s)", label, len(cmds), nthreads)
	joblst = []
	for path, cmd in cmds:
		os.makedirs(os.path.dirname(path), exist_ok=True)
		with open(path, "w") as f:
			f.write(cmd)
		joblst.append("bash " + path)
	with Pool(processes=nthreads) as pool:
		pool.map(runCMD, joblst)


# ---------------------------------------------------------------- entry

def germline_glimpse(args):
	logger.info("Performing germline variant calling (GLIMPSE2 backend)...")
	print_parameters_given(args)

	logger.info("Checking dependencies...")
	check_dependencies_glimpse(args)
	logger.info("Checking existence of essential resource files...")
	validate_user_setting_glimpse(args)

	out = os.path.abspath(args.out)
	os.makedirs(out + "/Script", exist_ok=True)
	os.makedirs(out + "/glimpse", exist_ok=True)
	chroms = _chroms(args)

	# Steps run in sequence, not as one big joblist: `phase` needs the
	# chunk file that `prepPanel` produces, so the scripts cannot all be
	# generated upfront the way the Beagle path does it.

	if args.step in ("prepPanel", "all"):
		todo = []
		for c in chroms:
			p = _paths(args, c)
			if args.force_panel or not os.path.isfile(p["chunks"]):
				todo.append((p["script"] + "_prepPanel.sh",
					_prep_panel_cmd(args, c)))
			else:
				logger.info("%s: binary panel already built, skipping "
					"(use --force-panel to rebuild)", c)
		_run(todo, args.nthreads, "prepPanel")

	if args.step in ("varScan", "all"):
		_run([(_paths(args, c)["script"] + "_varScan.sh",
			_varscan_cmd(args, c)) for c in chroms],
			args.nthreads, "varScan")

	if args.step in ("phase", "all"):
		todo = []
		for c in chroms:
			p = _paths(args, c)
			for i, cmd in enumerate(_phase_cmds(args, c)):
				todo.append((p["script"] + "_phase_{:03d}.sh".format(i), cmd))
		_run(todo, args.nthreads, "phase")

	if args.step in ("ligate", "all"):
		_run([(_paths(args, c)["script"] + "_ligate.sh", _ligate_cmd(args, c))
			for c in chroms], args.nthreads, "ligate")

	logger.info("GLIMPSE2 germline calling complete.")


def add_glimpse_subparser(subparsers, common_parser, argparse_mod):
	gg = subparsers.add_parser('germline-glimpse', parents=[common_parser],
		help='Germline calling with GLIMPSE2 instead of Beagle 4.1',
		formatter_class=argparse_mod.ArgumentDefaultsHelpFormatter)
	gg.add_argument('-o', '--out', required=True,
		help="Output directory (same one preProcess wrote to)")
	gg.add_argument('-g', '--reference', required=True,
		help="Genome reference FASTA (needs a .fai)")
	gg.add_argument('-p', '--imputation-panel', required=True,
		help="Directory of per-chromosome 1KG3 phased panel VCFs")
	gg.add_argument('--map-dir', required=True,
		help="Directory of GLIMPSE b38 genetic maps (chrN.b38.gmap.gz)")
	gg.add_argument('--panel-dir', required=True,
		help="Where to build the binary reference panel. Built once and "
		     "reused across every donor and every run -- put it somewhere "
		     "persistent, not in the per-run output directory.")
	gg.add_argument('-s', '--step', required=True,
		choices=['prepPanel', 'varScan', 'phase', 'ligate', 'all'],
		help="Run step by step")
	gg.add_argument('-c', '--chroms', required=False, default=None,
		help="Comma-separated subset, e.g. chr21,chr22. Default: chr1-chr22")
	gg.add_argument('-t', '--nthreads', required=False, type=int, default=1,
		help="Concurrent jobs. GLIMPSE2 is far lighter on memory than "
		     "Beagle -- budget ~2-4 GB per unit, not 20.")
	gg.add_argument('--phase-threads', required=False, type=int, default=1,
		help="Threads inside each GLIMPSE2_phase job. Keep "
		     "(nthreads x phase-threads) at or below your core allocation.")
	gg.add_argument('--bin-threads', required=False, type=int, default=2,
		help="Threads for panel prep, chunking and ligation")
	gg.add_argument('--min-mq', required=False, type=int, default=20,
		help="mpileup -q, minimum mapping quality")
	gg.add_argument('--min-bq', required=False, type=int, default=20,
		help="mpileup -Q, minimum base quality")
	gg.add_argument('--exclude-samples', required=False, default=None,
		help="Comma-separated panel samples to drop, e.g. "
		     "NA12878,NA12891,NA12892. REQUIRED for any benchmark against "
		     "a sample that is itself in the panel.")
	gg.add_argument('--force-panel', required=False, action='store_true',
		help="Rebuild the binary panel even if it already exists")
	gg.set_defaults(func=germline_glimpse)
	return gg