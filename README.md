# Duopogen

Germline SNV calling from ATAC and RNA sequencing, using linkage disequilibrium
from a reference panel to genotype sites individual reads cannot resolve.

A fork of [Monopogen](https://github.com/KChen-lab/Monopogen) (Dou et al.,
*Nat Biotechnol* 2024), reduced to the germline module and modernised:

- Python 3.14, current pysam / samtools / bcftools
- No vendored binaries except the Beagle 4.1 jar
- Two imputation backends: **GLIMPSE2** (recommended) and **Beagle 4.1**
  (upstream's, kept for reproducibility)
- Fails loudly instead of silently producing empty output

The somatic module, the MonoVar tree, and the R scripts are removed. For
somatic SNV calling, use upstream.

## Benchmark

**NB! The benchmark is preliminary!**

GM12878 bulk ATAC, genome-wide, scored against GIAB HG001 v4.2.1
high-confidence biallelic SNPs (n = 3,253,394). NA12878 and both parents were
removed from the reference panel first.

| | Beagle 4.1 | GLIMPSE2 |
|---|---|---|
| Variants recovered | 2,737,339 | **3,181,744** |
| Recall | 84.1% | **97.8%** |
| Non-ref discordance, shared sites | 0.514% | **0.500%** |
| Non-ref calls in HC regions | 2,731,302 | 3,159,988 |

**GLIMPSE2 is marginally more accurate on identical sites and recovers 444,405
more true variants.** 

Do not compare raw discordance or dosage r² between backends directly: they
call different site sets, so Beagle's headline r² of 0.979 against GLIMPSE2's
0.953 reflects site difficulty, not accuracy. Report shared-site discordance,
or report r² alongside recall.

---

## Contents

1. [Installation](#installation)
2. [Reference data](#reference-data)
3. [Pipeline](#pipeline)
4. [Pooled vs per-donor](#pooled-vs-per-donor)
5. [Validation](#validation)
6. [Troubleshooting](#troubleshooting)
7. [What changed from upstream](#what-changed-from-upstream)

---

## Installation

### 1. Clone

```bash
git clone https://github.com/mjesse-github/Duopogen.git
cd Duopogen
```

### 2. Conda environment

```bash
conda env create -f env.yml
conda activate duopogen_env
```

Two pins are load-bearing:

**`python_abi=3.14=*_cp314`** forces the standard (GIL) CPython build. Without
it the solver may pick the free-threaded variant, which has no pysam wheel —
pysam then compiles from source for several minutes and runs slower
single-threaded. The tell is `cp314t` ABI tags in pip output instead of
`cp314`, or pip printing `Building wheel for pysam`.

**`samtools`, `bcftools`, `htslib` all at 1.23**, matching the htslib that
pysam 0.24 bundles. Skew here causes BGZF and EOF-marker problems that surface
much later as truncated `.vcf.gz`.

Verify:

```bash
python -c "import sysconfig; print('free-threaded:', bool(sysconfig.get_config_var('Py_GIL_DISABLED')))"
python -c "import pysam; print(pysam.__version__, pysam.__samtools_version__)"
python src/Duopogen.py --help
```

Want `free-threaded: False`, `0.24.0 1.23.x`, and three subcommands.

### 3. GLIMPSE2

Deliberately **not** in `env.yml`. The bioconda package is a C++ binary linked
against specific htslib and boost versions; co-installing it can break the
htslib pin above or produce binaries that solve cleanly but fail at runtime
looking for the wrong shared object.

```bash
GLIMPSE_DIR=/path/to/persistent/glimpse2
mkdir -p ${GLIMPSE_DIR}/bin && cd ${GLIMPSE_DIR}/bin

for t in chunk split_reference phase ligate concordance; do
    wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_${t}_static
    mv GLIMPSE2_${t}_static GLIMPSE2_${t}
    chmod +x GLIMPSE2_${t}
done
```

Genetic maps are not in the release archive — they live in the repo:

```bash
cd $(dirname ${GLIMPSE_DIR})
git clone --depth 1 https://github.com/odelaneau/GLIMPSE.git glimpse2-src
ls glimpse2-src/maps/genetic_maps.b38/     # chr1.b38.gmap.gz ... chr22.b38.gmap.gz

conda activate duopogen_env
export PATH="${GLIMPSE_DIR}/bin:$PATH"
```

Verify — this also confirms the binaries actually link, which is the failure
mode the conda route would have given you:

```bash
for t in chunk split_reference phase ligate concordance; do
    printf '%-18s ' "$t"; GLIMPSE2_${t} --help >/dev/null 2>&1 && echo OK || echo FAIL
done
```

Static binaries are **Linux x86_64 only** — upstream publishes no arm64 or
macOS assets. Elsewhere, build from source or use a *separate* conda env:

```bash
conda create -n glimpse2 -c conda-forge -c bioconda "glimpse-bio>=2"
export PATH="$(conda info --base)/envs/glimpse2/bin:$PATH"
```

---

## Reference data

### Genome reference

GRCh38 **no-alt analysis set**, chr-prefixed, with decoy and HLA contigs.

Not the UCSC browser download — it keeps ALT contigs, which cause ambiguous
multi-mapping and break genotyping at HLA and other immune loci. Not the
bare-numbered NCBI default either — contig names would not match
CellRanger-ARC's `chr1`, `chr2`, … BAMs.

```bash
samtools faidx hg38.fa      # both backends require the .fai
```

### Imputation panel

1000 Genomes 30x NYGC phased panel, one VCF per chromosome, named as the code
expects:

```
CCDG_14151_B01_GRM_WGS_2020-08-05_<chrom>.filtered.shapeit2-duohmm-phased.vcf.gz
```

From the EBI 1000 Genomes FTP, `data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/`.
~520 MB per chromosome, ~11 GB total.

> **Benchmarking against a sample that is in the panel?** GM12878 is NA12878,
> and she and both parents (NA12891, NA12892) are among the 3,202. Strip the
> trio from a *copy* of the panel before measuring anything, or concordance is
> near-perfect by construction:
>
> ```bash
> bash strip_panel_trio.sh /path/to/panel_copy 8
> bcftools query -l /path/to/panel_copy/CCDG_..._chr22...vcf.gz | wc -l   # 3199
> ```

### Truth set (validation only)

```bash
BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
wget ${BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget ${BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget ${BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

The sample column is `HG001`, not `NA12878` — rename it to match your `bam.lst`
id before comparing. The `.bed` is not optional: outside it GIAB makes no
claim, so any concordance measured there is noise.

---

## Pipeline

Three stages, under `pipeline/`.

### 01 — build the donor manifest

Writes `donors.txt` and a per-donor list of source BAMs. Adjust the globbing to
your directory layout.

```bash
bash 01_prepare.sh
```

For a single-sample benchmark that is one donor with one or more BAMs; for a
cohort, one line per donor with however many libraries each.

### 02 — merge per donor

```bash
sbatch --array=1-$(wc -l < donors.txt)%10 02_merge_array.sh
```

Multiple libraries per donor are merged; a single library is symlinked, since
CellRanger-ARC output is already coordinate-sorted and re-sorting a 100+ GB BAM
produces a byte-identical file.

Two details that cost real debugging time:

- `samtools merge --write-index` writes **CSI**, not BAI. The code accepts
  both, but if you add tooling that assumes `.bai`, index explicitly with
  `samtools index`.
- The merged file is `${DONOR}.bam`, not `atac_possorted_bam.bam`. Every donor
  sharing a filename is what made pooling unsafe under upstream's read-group
  logic, which derived the sample name from the filename.

### 03 — call germline variants

The **id column of `bam.lst` becomes the VCF sample name**, so it must be
unique. Do not rely on filenames being distinct.

**GLIMPSE2 backend (recommended):**

```bash
python "${DUO}/src/Duopogen.py" preProcess \
    -b "${POOL}/bam.lst" -o "${POOL}/duopogen" \
    -t 10 --bam-threads 2

python "${DUO}/src/Duopogen.py" germline-glimpse \
    -o "${POOL}/duopogen" -g "${RES}/hg38.fa" -p "${RES}/1kg3_panel/" \
    --map-dir "${GLIMPSE}/maps/genetic_maps.b38" \
    --panel-dir "${RES}/glimpse_panel_bin" \
    -s all -t 14 --phase-threads 2
```

Steps run in order: `prepPanel` → `varScan` → `phase` → `ligate`.

**`--panel-dir` must be persistent and outside your run directory.** The
binarised reference panel is built once and reused by every donor and every
future run — that reuse is most of the speed advantage. Panel prep takes a few
hours genome-wide and is skipped on later runs unless you pass `--force-panel`.

Budget ~2–4 GB per unit of `-t`, not 20. Add `-c chr21,chr22` while testing.

**Beagle 4.1 backend (optional):**

```bash
python "${DUO}/src/Duopogen.py" germline \
    -o "${POOL}/duopogen" -g "${RES}/hg38.fa" -p "${RES}/1kg3_panel/" \
    -r "${DUO}/resource/GRCh38.region.lst" \
    -s all -t 14
```

Each region spawns a JVM at `-Xmx20g`, so budget ~20 GB per unit of `-t`.
`region.lst` accepts whole chromosomes (`chr1`) or chunks
(`chr1,1,50000001`); chunking raises parallelism above the 22-job ceiling and
lets Beagle restrict its panel read to the chunk.

Kept for reproducibility against upstream. Slower and less complete — see the
benchmark above.

### Concatenate

```bash
cd ${POOL}/duopogen/glimpse
bcftools concat -Oz -o $W/glimpse.all.vcf.gz chr{1..22}.phased.bcf
bcftools index -t $W/glimpse.all.vcf.gz
```

`chr{1..22}` expands numerically. A glob gives chr1, chr10, chr11, … chr2,
which produces an out-of-order file tabix refuses to index.

For the Beagle backend, substitute `chr{1..22}.phased.vcf.gz`.

---

## Pooled vs per-donor

Both work. **Per-donor is preferred with the GLIMPSE2 backend.**

### Per-donor (recommended)

One array job per donor, each writing a single-donor `bam.lst`:

```bash
DONOR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" donors.txt)
d="${WORKDIR}/${DONOR}"
printf '%s,%s\n' "$DONOR" "${d}/input/${DONOR}.bam" > "${d}/bam.lst"

python "${DUO}/src/Duopogen.py" preProcess -b "${d}/bam.lst" -o "${d}/duopogen" -t 10
python "${DUO}/src/Duopogen.py" germline-glimpse -o "${d}/duopogen" \
    -g "${RES}/hg38.fa" -p "${RES}/1kg3_panel/" \
    --map-dir "${GLIMPSE}/maps/genetic_maps.b38" \
    --panel-dir "${RES}/glimpse_panel_bin" -s all -t 8
```

```bash
sbatch --array=1-$(wc -l < donors.txt)%20 03_germline_per_donor.sh
```

### Pooled

```bash
: > "${POOL}/bam.lst"
while read -r donor; do
    printf '%s,%s\n' "$donor" "${WORKDIR}/${donor}/input/${donor}.bam" >> "${POOL}/bam.lst"
done < donors.txt
```

then one `preProcess` and one `germline-glimpse` against that list.

### Why per-donor wins with GLIMPSE2

Pooling exists to share a cost. With Beagle that cost was the reference model,
rebuilt from the 3,202-sample panel for every window of every run — so N
separate runs paid it N times and pooling paid it once.

GLIMPSE2 writes that model to disk. The `.bin` files under `--panel-dir` are
built once and read by every run, pooled or not. **There is no shared cost left
for pooling to capture**, and `GLIMPSE2_phase` handles samples independently, so
cost is linear in donors either way.

What pooling then costs you:

- **Serialisation.** One job is one job; N jobs fill the queue and use every
  core you can get. Cohort wall clock becomes roughly one donor's runtime.
- **All-or-nothing failure.** A bad BAM at donor 97 takes down all 150.
- **I/O pressure.** `bcftools mpileup -b` opens and streams every BAM in the
  list concurrently. At 150 donors that is where you hit a wall first, and it
  is file handles and I/O rather than CPU.
- **Awkward reruns.** Adding three donors later means re-running mpileup on
  all 153.

What pooling still buys: a single multi-sample VCF, so no merge afterwards. If
you want a joint callset for QTL work that is a genuine convenience — but
merging per-donor output is straightforward:

```bash
bcftools merge -Oz -o cohort.vcf.gz */glimpse.all.vcf.gz
bcftools index -t cohort.vcf.gz
```

**Rule of thumb:** ≤10 donors, either is fine — pool if you want joint output.
Above that, run per-donor, or batch 10–20 per job to get parallelism and
restartability without producing 150 files to merge.

### Beagle is different, and worse

Pooling helped Beagle in principle but **backfired in practice**. With 10
donors pooled, nearly every panel site gained coverage in some donor, so the
target marker set expanded to fill the reference panel — `target markers:
49958` against `reference markers: 50000`. More markers means more windows, and
each window rebuilds the model. Measured: ~2 h per window, ~18 windows for
chr1, paid twice. chr1 alone exceeded a 24 h walltime.

If you use the Beagle backend at cohort scale, split by step (`-s varScan`,
then `-s varImpute`, then `-s varPhasing`) with `afterok` dependencies so each
job stays under the walltime cap.

---

## Validation

```bash
bcftools view -H sample.phased.vcf.gz | wc -l
bcftools view -H sample.phased.vcf.gz | cut -f10 | cut -d: -f1 | sort | uniq -c
```

Het / hom-alt around 1.3–1.6 for a European sample; `0|1` and `1|0` roughly
balanced; all genotypes using `|` not `/` — a `/` means phasing did not run.

Against a truth set:

```bash
bash benchmark.sh
```

Reports recall and discordance separately, plus discordance restricted to sites
both backends called. That restriction matters: GLIMPSE2 reports every panel
site while Beagle's `impute=false` reports only observed ones, so raw
discordance grades Beagle on the easier subset it chose to answer.

---

## What changed from upstream

**Removed:** the somatic module, the R scripts, the orphaned MonoVar tree (11
files reachable only from `monovar.py`, which nothing imports), and every
vendored binary in `apps/` except the Beagle jar. Upstream shipped samtools 1.2
(2015), bcftools 1.8 (2018) and `libcrypto.so.1.0.0`, which users had to prepend
to `LD_LIBRARY_PATH` — shadowing the environment's own OpenSSL. `vcftools` and
`picard.jar` were asserted to exist but never invoked.

**Fixed:** an uninitialised variable letting a read inherit the previous read's
mismatch count; a `runCMD` that raised `NameError` on success and swallowed
failures; script files opened per region but closed once; a `zless` dependency
absent from compute nodes; missing `##contig` headers; hardcoded `range(1, 23)`
silently dropping chrX; read groups derived from filenames rather than sample
ids.

**Added:** the GLIMPSE2 backend; genotype likelihoods computed at all panel
sites via `bcftools call -C alleles`, keeping hom-ref evidence the Beagle path
discards; a genetic map, where the Beagle runs log `No genetic map is specified:
using 1 cM = 1 Mb`.

**Changed:** binaries resolve from `PATH`; the read filter runs in samtools
rather than a Python loop; `multiprocessing` pins the `fork` start method, since
Python 3.14 changed the POSIX default to `forkserver`.

**Unchanged:** `modelscale=2`, `niterations=0`, `impute=false` in
the Beagle calls, and `--ns 0` in both pileups. The Beagle parameters affect
accuracy, not just speed. `--ns 0` disables bcftools' default skipping of
`UNMAP,SECONDARY,QCFAIL,DUP` — inherited from upstream, and worth revisiting for
single-cell data, where duplicate rates are far higher than in bulk and
including duplicates makes genotype likelihoods overconfident.

---

## Citation

> Dou J, Tan Y, Kock KH, *et al.* Single-nucleotide variant calling in
> single-cell sequencing data with Monopogen. *Nature Biotechnology* 42,
> 803–812 (2024).

If using the GLIMPSE2 backend:

> Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. Imputation of
> low-coverage sequencing data from 150,119 UK Biobank genomes.
> *Nature Genetics* 55, 1088–1090 (2023).