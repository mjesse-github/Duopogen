#!/bin/bash
#SBATCH --job-name=duo_merge
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH -c 40
#SBATCH --account=ACCOUNT
#SBATCH --partition=PARTITION
#SBATCH --output=logs/slurm-%x.%A_%a.out
#SBATCH --mail-type=END,FAIL
# sbatch --array=1-$(wc -l < donors.txt)%10 02_merge_array.sh

set -euo pipefail
conda activate duopogen_env

WORKDIR="<where workdir is>/duopogen_workdir"
DONOR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" donors.txt)
ddir="${WORKDIR}/${DONOR}"
merged="${ddir}/input/${DONOR}.bam"     # named per donor, not atac_possorted_bam.bam
mapfile -t bams < "${ddir}/input/.source_bams.txt"

echo "[$DONOR] merging ${#bams[@]} BAM(s)"

if [[ ${#bams[@]} -eq 1 ]]; then
    #if one sample for donor, just symlink to avoid samtools merge overhead
    ln -sf "$(realpath "${bams[0]}")"     "$merged"
    ln -sf "$(realpath "${bams[0]}").bai" "${merged}.bai"
else
    samtools merge -@ 40 -f -c -p --write-index -o "$merged" "${bams[@]}"
fi

samtools quickcheck -v "$merged" && echo "[$DONOR] merge OK"
