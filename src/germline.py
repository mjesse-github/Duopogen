#!/usr/bin/env python3
"""
germline.py -- germline SNV calling helpers for duopogen.

Fork of Monopogen (KChen-lab), germline-only. Somatic modules removed.

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
import subprocess
import pysam
import multiprocessing as mp
from multiprocessing import Pool

# [ENV] Dropped: pandas, numpy, gzip, glob, re, time, VariantFile.
#       None are used once the monovar tree is deleted. numpy+pandas alone
#       cost ~1.5 s of interpreter startup, paid once per region job.

# [ENV] The only thing still vendored is the Beagle 4.1 jar. Beagle 5.x
#       dropped gl= (genotype-likelihood) input, which this pipeline
#       depends on, so bioconda's `beagle` package is NOT a substitute.
#       Everything else resolves from the active conda env via PATH.
APP_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "apps"))
BEAGLE_JAR = os.path.join(APP_PATH, "beagle.27Jul16.86a.jar")

# [PERF] Upstream hardcodes `for chr in range(1, 23)` in two places, which
#        silently drops chrX/chrY. Add "chrX" here (and a matching panel
#        file + region.lst entry) if you want it. Note chrX needs ploidy
#        handling that Beagle will not infer on its own.
CHROMS = ["chr" + str(i) for i in range(1, 23)]

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg == "func":
			continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))


def validate_user_setting_germline(args):
	assert os.path.isfile(args.reference), \
		"The genome reference fasta file {} cannot be found!".format(args.reference)
	assert os.path.isdir(args.imputation_panel), \
		"Imputation panel directory {} cannot be found!".format(args.imputation_panel)
	assert os.path.isfile(args.region), \
		"The region file {} cannot be found!".format(args.region)

	# [PERF] uses CHROMS instead of range(1, 23)
	for chrom in CHROMS:
		bamFile = args.out + "/Bam/" + chrom + ".filter.bam.lst"
		assert os.path.isfile(bamFile), \
			"{} not found -- did preProcess complete?".format(bamFile)
		with open(bamFile) as f_in:
			for line in f_in:
				line = line.strip()
				assert os.path.isfile(line), \
					"The bam file {} cannot be found!".format(line)
				assert (os.path.isfile(line + ".bai")
					or os.path.isfile(line + ".csi")), \
					"Index for {} cannot be found!".format(line)

	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			assert len(record) == 3 or len(record) == 1, \
				("Every line needs exactly 3 comma-delimited columns "
				 "(chr1,1,100000) or 1 (chr1). Offending line: {}".format(line))


def check_dependencies(args):
	# [ENV] Upstream asserted vcftools and picard.jar exist in apps/. Neither
	#       is ever invoked anywhere in the codebase -- they were pure
	#       friction. Removed. The remaining four are genuinely called.
	for prog in ("bgzip", "bcftools", "samtools", "java"):
		path = shutil.which(prog)
		assert path is not None, \
			"Program {} not found on PATH -- is duopogen_env active?".format(prog)
		logger.debug("Using {} -> {}".format(prog, path))
	assert os.path.isfile(BEAGLE_JAR), \
		"Beagle jar not found at {}".format(BEAGLE_JAR)


def runCMD(cmd):
	# [FIX] Upstream read:
	#           output = os.system(cmd)
	#           if output == 0:
	#               return(region)      <- NameError, `region` is undefined
	#       It only triggered on SUCCESS, and lived in somatic.py where a
	#       working duplicate shadowed it. Deleting somatic.py exposes it.
	#       Upstream also swallowed nonzero exits entirely, so a failed
	#       Beagle run surfaced as a missing file several steps later.
	output = os.system(cmd)
	assert output == 0, "command failed (exit {}): {}".format(output, cmd)
	return cmd


def addChr(in_bam, samtools):
	"""Prefix contig names with 'chr'. Only runs for BAMs aligned to a
	bare-numbered reference. CellRanger-ARC output is already chr-prefixed,
	so this is a no-op path for most users -- left on the slow pysam
	implementation deliberately."""
	prefix = 'chr'
	out_bam = in_bam + "tmp.bam"
	input_bam = pysam.AlignmentFile(in_bam, "rb")
	new_head = input_bam.header.to_dict()
	for seq in new_head['SQ']:
		seq['SN'] = prefix + seq['SN']
	with pysam.AlignmentFile(out_bam, "wb", header=new_head) as outf:
		for read in input_bam.fetch():
			a = pysam.AlignedSegment(outf.header)
			a.query_name = read.query_name
			a.query_sequence = read.query_sequence
			a.reference_name = prefix + read.reference_name
			a.flag = read.flag
			a.reference_start = read.reference_start
			a.mapping_quality = read.mapping_quality
			a.cigar = read.cigar
			a.next_reference_id = read.next_reference_id
			a.next_reference_start = read.next_reference_start
			a.template_length = read.template_length
			a.query_qualities = read.query_qualities
			a.tags = read.tags
			outf.write(a)
	input_bam.close()
	subprocess.run([samtools, "index", out_bam], check=True)
	os.replace(out_bam, in_bam)
	os.replace(out_bam + ".bai", in_bam + ".bai")


def BamFilter(myargs):
	"""Keep reads with < max_mismatch alignment mismatches, per chromosome,
	and stamp a read group derived from the sample id."""
	bamFile      = myargs["bamFile"]
	search_chr   = myargs["chr"]
	samtools     = myargs["samtools"]
	chrom        = search_chr
	sample_id    = myargs["id"]
	max_mismatch = myargs["max_mismatch"]
	out          = myargs["out"]
	threads      = myargs.get("bam_threads", 2)

	os.makedirs(out + "/Bam", exist_ok=True)

	# Detect whether the BAM uses chr-prefixed contigs.
	with pysam.AlignmentFile(bamFile, "rb") as probe:
		has_chr = any(c.startswith("chr") for c in probe.references)
	if not has_chr:
		logger.info("Contigs lack a 'chr' prefix; will add it after filtering.")
		search_chr = search_chr[3:]

	outbam = "{}/Bam/{}_{}.filter.bam".format(out, sample_id, chrom)

	# The RG that mpileup uses to name the sample. Upstream derived SM/ID
	# from the BAM *filename*, which collides when several donors share a
	# filename (e.g. every donor's atac_possorted_bam.bam). Using the id
	# column from bam.lst instead.  [FIX]
	rg = "ID:{0}\\tSM:{0}\\tLB:0.1\\tPL:ILLUMINA\\tPU:{0}".format(sample_id)

	# [PERF] Upstream looped over every read in Python via pysam, rebuilding
	#        each record, then wrote it back out -- roughly 1-3 us/read plus
	#        full decompress/recompress in the interpreter. samtools does the
	#        same work in C with threaded BGZF: expect ~20-50x.
	#
	#        The filter expression also encodes the [FIX] for upstream's
	#        uninitialised `val`: upstream used two independent `if`s, so a
	#        read carrying neither NM nor nM inherited the PREVIOUS read's
	#        mismatch count (or raised NameError on record one). Here a read
	#        with neither tag matches neither clause and is dropped.
	#
	#        Verify against the old path once with:
	#            samtools view -c old.bam ; samtools view -c new.bam
	expr = "([NM] < {0}) || ([nM] < {0})".format(max_mismatch)
	view = [samtools, "view", "-@", str(threads), "-b",
	        "-e", expr, bamFile, search_chr]
	addrg = [samtools, "addreplacerg", "-@", str(threads),
	         "-r", rg, "-o", outbam, "-"]

	p1 = subprocess.Popen(view, stdout=subprocess.PIPE)
	p2 = subprocess.Popen(addrg, stdin=p1.stdout)
	p1.stdout.close()
	p2.communicate()
	assert p1.wait() == 0, "samtools view failed for {} {}".format(sample_id, chrom)
	assert p2.returncode == 0, "samtools addreplacerg failed for {} {}".format(sample_id, chrom)

	subprocess.run([samtools, "index", "-@", str(threads), outbam], check=True)
	if not has_chr:
		addChr(outbam, samtools)
	return outbam


def robust_get_tag(read, tag_name):
	try:
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"