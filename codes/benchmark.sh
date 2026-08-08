#!/usr/bin/env bash
#
# benchmark.sh -- score Beagle and GLIMPSE2 callsets against GIAB NA12878.
#
#   bash benchmark.sh
#
# Edit the five paths below. Runs on a login node in ~15 min.

set -euo pipefail

W=<where workdir is>/duopogen_workdir
T=<where resources are>/resource/giab
G=$T/HG001_GRCh38_1_22_v4.2.1_benchmark
SAMPLE=GM
TOOLS="beagle glimpse"

cd "$W"

# ---- truth ------------------------------------------------------------
# Restricted to GIAB high-confidence regions (outside them GIAB makes no
# claim), biallelic SNPs only (the panel is SNP-only, GIAB includes indels),
# and renamed HG001 -> $SAMPLE so bcftools can pair the samples.
if [[ ! -f truth.vcf.gz ]]; then
    echo "HG001 ${SAMPLE}" > rename.txt
    bcftools view -R $G.bed -v snps -m2 -M2 -Ou $G.vcf.gz \
      | bcftools reheader -s rename.txt \
      | bcftools view -Oz -o truth.vcf.gz
    bcftools index -t truth.vcf.gz
fi
NTRUTH=$(bcftools index -n truth.vcf.gz)

# ---- per tool ---------------------------------------------------------
for f in $TOOLS; do
    # same filter as truth, or you are comparing different site sets
    if [[ ! -f $f.hc.vcf.gz ]]; then
        bcftools view -R $G.bed -v snps -m2 -M2 -Oz -o $f.hc.vcf.gz $f.all.vcf.gz
        bcftools index -t $f.hc.vcf.gz
    fi
    # sites the tool recovered, and non-ref calls (for the false-positive check)
    eval "N_$f=\$(bcftools isec -n=2 -w1 truth.vcf.gz $f.hc.vcf.gz | grep -vc '^#')"
    eval "ALT_$f=\$(bcftools view -H -i 'GT=\"alt\"' $f.hc.vcf.gz | wc -l)"
    # discordance and dosage r2
    eval "GC_$f=\$(bcftools stats -s $SAMPLE truth.vcf.gz $f.hc.vcf.gz 2>/dev/null | awk '\$1==\"GCsS\"')"
done

# ---- fair comparison: score each tool on the OTHER's site set ----------
# GLIMPSE reports every panel site; Beagle only observed ones. Raw
# discordance therefore flatters Beagle, which is graded on the easier
# subset it chose to answer. This restricts both to the shared sites.
if [[ ! -f shared.vcf.gz ]]; then
    bcftools isec -n=2 -w1 -Oz -o shared.sites.vcf.gz beagle.hc.vcf.gz glimpse.hc.vcf.gz
    bcftools index -t shared.sites.vcf.gz
fi
for f in $TOOLS; do
    if [[ ! -f $f.shared.vcf.gz ]]; then
        bcftools view -R shared.sites.vcf.gz -Oz -o $f.shared.vcf.gz $f.hc.vcf.gz
        bcftools index -t $f.shared.vcf.gz
    fi
    eval "SH_$f=\$(bcftools stats -s $SAMPLE truth.vcf.gz $f.shared.vcf.gz 2>/dev/null | awk '\$1==\"GCsS\"')"
done

# ---- report -----------------------------------------------------------
echo
echo "GIAB high-confidence biallelic SNPs: $NTRUTH"
echo
printf '%-10s %10s %8s %10s %8s %12s %10s\n' \
    TOOL RECOVERED RECALL NRD r2 "NRD(shared)" "NONREF"
for f in $TOOLS; do
    n=$(eval echo \$N_$f); alt=$(eval echo \$ALT_$f)
    gc=$(eval echo \$GC_$f); sh=$(eval echo \$SH_$f)
    printf '%-10s %10d %7.1f%% %9s%% %8s %11s%% %10d\n' \
        "$f" "$n" "$(echo "scale=4; 100*$n/$NTRUTH" | bc)" \
        "$(echo $gc | cut -d' ' -f4)" "$(echo $gc | awk '{print $NF}')" \
        "$(echo $sh | cut -d' ' -f4)" "$alt"
done
echo
echo "RECALL      how many true variants each tool recovered"
echo "NRD         non-ref discordance on the sites it called (lower better)"
echo "NRD(shared) same, restricted to sites BOTH tools called -- the fair one"
echo "NONREF      non-ref calls in HC regions; well above $NTRUTH suggests false positives"