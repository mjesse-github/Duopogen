#!/bin/bash
#SBATCH --job-name=duo_germline
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH -c 32
#SBATCH --account=ACCOUNT
#SBATCH --partition=PARTITION
#SBATCH --output=logs/slurm-%x.%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail
conda activate duopogen_env

DUO="<where duopogen is installed>/Duopogen"
RES="<where resources are>/resource"
WORKDIR="<where workdir is>/duopogen_workdir"
POOL="${WORKDIR}/_mono_pooled"
mkdir -p "${POOL}"

#all the cores are arbitrarily set here, most likely is not optimal

# One bam.lst with every donor. The id column becomes the VCF sample name --
# this is why the @RG rewrite in BamFilter had to stop deriving SM from the
# filename: every donor's source BAM was called atac_possorted_bam.bam.
: > "${POOL}/bam.lst"
while read -r donor; do
    bam="$(realpath "${WORKDIR}/${donor}/input/${donor}.bam")"
    [[ -f "$bam" ]] || { echo "MISSING merged BAM for ${donor}" >&2; exit 1; }
    printf '%s,%s\n' "$donor" "$bam" >> "${POOL}/bam.lst"
done < donors.txt
echo "Pooling $(wc -l < "${POOL}/bam.lst") donors"

# -t 10: 10 concurrent per-chromosome filter jobs x 2 BGZF threads = 20 cores.
python "${DUO}/src/Duopogen.py" preProcess \
    -b "${POOL}/bam.lst" -o "${POOL}/duopogen" \
    -t 10 --bam-threads 2

# hg38.fa = GRCh38 no-alt analysis set, chr-prefixed (UCSC-style contig names).
# NOT the UCSC browser download (keeps ALT contigs -> ambiguous multi-mapping,
# breaks genotyping at HLA/immune loci). NOT bare-numbered NCBI default either
# (contig names wouldn't match CellRanger-ARC's chr1/chr2/... BAMs).

# -t 14: each region job spawns a JVM at -Xmx8g -> 112 GB peak against the
# 120 GB request, with mpileup at -d 1000 alongside. Beagle nthreads=2 inside
# each, so 14 x 2 = 28 cores of the 32 requested.
python "${DUO}/src/Duopogen.py" germline \
    -t 14 \
    -r "${DUO}/resource/GRCh38.region.lst" \
    -p "${RES}/1kg3_panel/" \
    -g "${RES}/hg38.fa" \
    -s all -o "${POOL}/duopogen"

echo "=== germline done ==="
