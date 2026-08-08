#!/bin/bash
#Run from login

# Build donors.txt and per-donor source BAM lists from direct directory contents.
set -euo pipefail

# Configuration
FRAGS_DIR="Where fragmetns are for donor"
WORKDIR="<where workdir is>/duopogen_workdir"
DONOR_NAME="<Donor ID, should be unique>"

declare -A donor_bams

# Find all matching BAMs directly in the frags directory
shopt -s nullglob
bams=("${FRAGS_DIR}"/*.bam "${FRAGS_DIR}"/*/atac_possorted_bam.bam)
shopt -u nullglob

if [[ ${#bams[@]} -eq 0 ]]; then
    echo "WARNING: No BAM files found in ${FRAGS_DIR}" >&2
fi

for bam in "${bams[@]}"; do
    if [[ ! -f "${bam}.bai" ]]; then
        echo "WARNING: ${bam} is not indexed" >&2
    fi
    donor_bams["$DONOR_NAME"]+="$bam "
done

mkdir -p "$WORKDIR"
: > donors.txt
for donor in $(printf '%s\n' "${!donor_bams[@]}" | sort); do   # sorted = stable array indices
    echo "$donor" >> donors.txt
    mkdir -p "${WORKDIR}/${donor}/input"
    printf '%s\n' ${donor_bams[$donor]} > "${WORKDIR}/${donor}/input/.source_bams.txt"
done

echo "Wrote donors.txt with $(wc -l < donors.txt) donors"
echo "Total source BAMs: $(cat ${WORKDIR}/*/input/.source_bams.txt | wc -l)"