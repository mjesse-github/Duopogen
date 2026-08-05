#!/usr/bin/env python3
"""
Monopogen.py -- main interface for duopogen (germline-only Monopogen fork).

Subcommands: preProcess, germline.

Changes from upstream are tagged in comments:
  [FIX]   correctness bug in upstream
  [PERF]  speed / memory
  [3.14]  needed for Python 3.14
  [ENV]   binaries now come from PATH, not a vendored apps/ dir
"""

import argparse
import sys
import os
import logging
import shutil
import multiprocessing as mp
from multiprocessing import Pool

from germline import *

# [3.14] Python 3.14 changed the POSIX default start method from fork to
#        forkserver (bpo/gh-84559). Upstream sets `samtools`, `bcftools`,
#        `beagle` etc. as module globals inside main() and assumes Pool
#        workers inherit them. The germline path happens to survive
#        forkserver -- runCMD takes a string, BamFilter takes a dict with
#        paths already resolved -- but that is luck, not design. Pin fork.
try:
	mp.set_start_method('fork', force=True)
except RuntimeError:
	pass

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def germline(args):
	logger.info("Performing germline variant calling...")
	print_parameters_given(args)

	logger.info("Checking existence of essential resource files...")
	validate_user_setting_germline(args)

	logger.info("Checking dependencies...")
	check_dependencies(args)

	out = os.path.abspath(args.out)
	os.makedirs(out + "/germline", exist_ok=True)
	os.makedirs(out + "/Script", exist_ok=True)

	joblst = []
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if not record or not record[0]:
				continue
			if len(record) == 1:
				jobid = record[0]
			else:
				jobid = record[0] + ":" + record[1] + "-" + record[2]

			bam_filter = out + "/Bam/" + record[0] + ".filter.bam.lst"
			with open(bam_filter) as p:
				N_sample = sum(1 for _ in p)

			imputation_vcf = (args.imputation_panel
				+ "CCDG_14151_B01_GRM_WGS_2020-08-05_" + record[0]
				+ ".filtered.shapeit2-duohmm-phased.vcf.gz")
			assert os.path.isfile(imputation_vcf), \
				"Panel file not found: {}".format(imputation_vcf)

			# ---- variant scan -------------------------------------------
			# [FIX]  Upstream wrote `-b" + bam_filter` with no space. getopt
			#        tolerates -bFILE, but it is one edit from breaking.
			# [PERF] -d 10000000 -> 1000. The 10M cap is what produces
			#        upstream's own "(mpileup) Max depth is above 1M.
			#        Potential memory hog!" warning. At scRNA/scATAC depths
			#        nothing is lost, and it removes the main OOM source
			#        when several regions run concurrently.
			cmd1 = (samtools + " mpileup -b " + bam_filter
				+ " -f " + args.reference + " -r " + jobid
				+ " -q 20 -Q 20 --incl-flags 0 --excl-flags 0 -t DP -d 1000 -v ")
			cmd1 += (" | " + bcftools + " view "
				+ " | " + bcftools + ' filter -e \'REF !~ "^[ATGC]$"\' '
				+ " | " + bcftools + " norm -m-both -f " + args.reference)
			cmd1 += (" | grep -v \"<X>\" | " + bgzip + " -c > "
				+ out + "/germline/" + jobid + ".gl.vcf.gz")

			# ---- imputation ---------------------------------------------
			# [PERF] chrom= now carries the CHUNK interval, not just the
			#        chromosome. This is the single biggest germline win.
			#        resource/GRCh38.region.lst already splits at 50 Mb, but
			#        upstream passed chrom=chr1 for all five chr1 chunks, so
			#        Beagle parsed the ENTIRE chr1 reference panel five
			#        times. Restricting it cuts reference parsing roughly
			#        5x on the large chromosomes and is what makes the
			#        smaller heap below safe.
			#        (Beagle 4.1 accepts chrom=<CHROM>[:<START>-<END>].)
			#
			# [PERF] -Xmx20g -> 8g. 20g was sized for a whole-chromosome
			#        window. At 50 Mb with 3202 reference samples, 8g is
			#        ample -- and it is the difference between 5 and 12
			#        concurrent JVMs in a 100 GB allocation.
			#
			# [PERF] nthreads 24 -> 2. Beagle 4.1 parallelises over TARGET
			#        samples; with one donor, 24 threads bought nothing and
			#        oversubscribed the node once -t > 1. Keep
			#        (nthreads x -t) at or below your -c allocation.
			#
			# NOT changed: modelscale=2, niterations=0, impute=false.
			# All three affect accuracy, not just speed. Do not touch them
			# without a WGS-truth comparison.
			cmd3 = (java + " -Xmx8g -jar " + BEAGLE_JAR
				+ " gl=" + out + "/germline/" + jobid + ".gl.vcf.gz"
				+ " ref=" + imputation_vcf
				+ "  chrom=" + jobid
				+ " out=" + out + "/germline/" + jobid + ".gp "
				+ "impute=false  modelscale=2  nthreads=2  gprobs=true  niterations=0")

			# [FIX] Upstream used `zless -S`, which is absent on most compute
			#       nodes and was only reachable because apps/ shipped a gzip
			#       symlink. zcat is universally present.
			#       (Both branches of upstream's N_sample if/elif were
			#       identical, so the branch is gone -- N_sample is kept
			#       only for the log line below.)
			cmd4 = ("zcat " + out + "/germline/" + jobid + ".gp.vcf.gz > "
				+ out + "/germline/" + jobid + ".germline.vcf")

			cmd5 = (java + " -Xmx8g -jar " + BEAGLE_JAR
				+ " gt=" + out + "/germline/" + jobid + ".germline.vcf"
				+ " ref=" + imputation_vcf
				+ "  chrom=" + jobid
				+ " out=" + out + "/germline/" + jobid + ".phased "
				+ "impute=false  modelscale=2  nthreads=2  gprobs=true  niterations=0")
			cmd5 += "\nrm -f " + out + "/germline/" + jobid + ".germline.vcf"

			logger.debug("region {} ({} sample(s))".format(jobid, N_sample))

			# [FIX] Upstream opened one script file per region inside the
			#       loop and called f_out.close() ONCE, after the loop --
			#       so only the last file was explicitly closed. CPython's
			#       refcounting papered over it; the failure mode is a
			#       truncated runGermline_*.sh that then looks like a
			#       Beagle error. `with` closes each one.
			script = out + "/Script/runGermline_" + jobid + ".sh"
			with open(script, "w") as f_out:
				f_out.write("set -euo pipefail\n")
				if args.step in ("varScan", "all"):
					f_out.write(cmd1 + "\n")
				if args.step in ("varImpute", "all"):
					f_out.write(cmd3 + "\n")
					f_out.write(cmd4 + "\n")
				if args.step in ("varPhasing", "all"):
					f_out.write(cmd5 + "\n")

			joblst.append("bash " + script)

	if args.norun == "TRUE":
		logger.info("Generated {} job scripts in {}/Script/ (not run)."
			.format(len(joblst), out))
		return

	# runCMD now asserts on nonzero exit, so a failure here raises rather
	# than silently producing a missing file three steps downstream.  [FIX]
	with Pool(processes=args.nthreads) as pool:
		pool.map(runCMD, joblst)


def preProcess(args):
	logger.info("Performing data preprocess before variant calling...")
	print_parameters_given(args)

	assert os.path.isfile(args.bamFile), \
		"The bam list file {} cannot be found!".format(args.bamFile)
	out = os.path.abspath(args.out)
	os.makedirs(out + "/Bam", exist_ok=True)

	samples = []
	para_lst = []
	with open(args.bamFile) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if not record or not record[0]:
				continue
			assert len(record) == 2, \
				("Every line needs exactly 2 comma-delimited columns "
				 "(id,bam). Offending sample: {}".format(record[0]))
			assert os.path.isfile(record[1]), \
				"Bam file {} cannot be found!".format(record[1])
			assert os.path.isfile(record[1] + ".bai"), \
				"Bam file {} has not been indexed!".format(record[1])
			assert record[0] not in samples, \
				"Duplicate sample id {} in {}".format(record[0], args.bamFile)
			samples.append(record[0])

			logger.debug("PreProcessing sample {}".format(record[0]))
			for chrom in CHROMS:
				para_lst.append(dict(
					chr=chrom,
					out=out,
					id=record[0],
					bamFile=record[1],
					max_mismatch=args.max_mismatch,
					samtools=samtools,
					# [PERF] samtools does its own BGZF threading now. Keep
					#        (bam_threads x nthreads) near your -c count.
					bam_threads=args.bam_threads,
				))

	with Pool(processes=args.nthreads) as pool:
		pool.map(BamFilter, para_lst)

	# [PERF] uses CHROMS instead of range(1, 23) -- see germline.py
	for chrom in CHROMS:
		with open(out + "/Bam/" + chrom + ".filter.bam.lst", "w") as bamlist:
			for s in samples:
				bamlist.write("{}/Bam/{}_{}.filter.bam\n".format(out, s, chrom))


def main():
	parser = argparse.ArgumentParser(
		description="duopogen: germline SNV calling from single cell sequencing",
		epilog="Typical modules: preProcess, germline",
		formatter_class=argparse.RawTextHelpFormatter)

	subparsers = parser.add_subparsers(title='Available subcommands',
		dest="subcommand")
	common_parser = argparse.ArgumentParser(add_help=False)

	# [ENV] -a/--app-path is gone from both subparsers. Binaries come from
	#       PATH; the Beagle jar is located relative to this file.

	pp = subparsers.add_parser('preProcess', parents=[common_parser],
		help='Filter reads with high alignment mismatch, split by chromosome',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	pp.add_argument('-b', '--bamFile', required=True,
		help="CSV list of samples: one 'id,/abs/path.bam' per line")
	pp.add_argument('-o', '--out', required=True, help="Output directory")
	pp.add_argument('-m', '--max-mismatch', required=False, type=int, default=3,
		help="Maximal alignment mismatch allowed in one read")
	pp.add_argument('-t', '--nthreads', required=False, type=int, default=1,
		help="Number of concurrent per-chromosome filter jobs")
	pp.add_argument('--bam-threads', required=False, type=int, default=2,
		help="BGZF threads passed to samtools within each job")
	pp.set_defaults(func=preProcess)

	gl = subparsers.add_parser('germline', parents=[common_parser],
		help='Germline variant discovery and genotype calling',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	gl.add_argument('-r', '--region', required=True,
		help="Region file: 'chr1' or 'chr1,1,50000001' per line")
	gl.add_argument('-s', '--step', required=True, default="all",
		choices=['varScan', 'varImpute', 'varPhasing', 'all'],
		help="Run germline variant calling step by step")
	gl.add_argument('-o', '--out', required=True, help="Output directory")
	gl.add_argument('-g', '--reference', required=True,
		help="Genome reference FASTA used for alignment")
	gl.add_argument('-p', '--imputation-panel', required=True,
		help="Directory of per-chromosome 1KG3 phased panel VCFs")
	gl.add_argument('-m', '--max-softClipped', required=False, type=int, default=1,
		help="Maximal soft-clipping allowed in one read")
	gl.add_argument('-t', '--nthreads', required=False, type=int, default=1,
		help="Number of concurrent region jobs. Each spawns a JVM at -Xmx8g, "
		     "so budget ~8 GB per unit of -t.")
	gl.add_argument('-n', '--norun', required=False, default="FALSE",
		choices=['TRUE', 'FALSE'],
		help="Generate job scripts only; do not run them")
	gl.set_defaults(func=germline)

	args = parser.parse_args()
	if args.subcommand is None:
		print("Please specify one subcommand! Exiting!")
		print("-" * 80)
		parser.print_help()
		sys.exit(1)

	# [ENV] Resolve tools from the active environment rather than from a
	#       vendored apps/ directory. Upstream shipped samtools 1.2 (2015)
	#       and bcftools 1.8 (2018) alongside libcrypto.so.1.0.0, which
	#       users then had to prepend to LD_LIBRARY_PATH -- shadowing the
	#       env's own OpenSSL. Deleting apps/ removes that hazard entirely.
	global out, samtools, bcftools, bgzip, java
	out = os.path.abspath(args.out)

	def _which(prog):
		p = shutil.which(prog)
		assert p is not None, \
			"{} not found on PATH -- is duopogen_env active?".format(prog)
		return p

	samtools = _which("samtools")
	bcftools = _which("bcftools")
	bgzip    = _which("bgzip")
	java     = _which("java")

	args.func(args)
	logger.info("Success! See instructions above.")


if __name__ == "__main__":
	main()