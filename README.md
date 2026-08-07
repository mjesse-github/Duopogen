# Duopogen

Germline SNV calling from single-cell / bulk ATAC and RNA sequencing, using
linkage disequilibrium from an external reference panel to genotype sites that
individual reads cannot resolve on their own.

A fork of [Monopogen](https://github.com/KChen-lab/Monopogen) (Dou et al.,
*Nat Biotechnol* 2024), reduced to the germline module and modernised:

- Runs on Python 3.14 with current pysam / samtools / bcftools
- No vendored binaries except the Beagle 4.1 jar — everything else resolves
  from the active conda environment
- Two interchangeable imputation backends: **Beagle 4.1** (upstream's) and
  **GLIMPSE2** (faster, and better suited to many donors)
- Fails loudly instead of silently producing empty output

The somatic module, the MonoVar tree, and the R scripts are **removed**. If you
need somatic SNV calling, use upstream.

---

## Contents

1. [Installation](#installation)
2. [Reference data](#reference-data)
3. [Running the pipeline](#running-the-pipeline)
4. [Choosing a backend](#choosing-a-backend)
5. [Output](#output)
6. [Troubleshooting](#troubleshooting)
7. [What changed from upstream](#what-changed-from-upstream)

---

## Installation

### 1. Clone

```bash
git clone https://github.com/mjesse-github/Duopogen.git
cd Duopogen
```

### 2. Create the conda environment

```bash
conda env create -f env.yml
conda activate duopogen_env
```

Two pins in `env.yml` are load-bearing and should not be relaxed casually:

- **`python_abi=3.14=*_cp314`** forces the standard (GIL) CPython build.
  Without it the solver may select the free-threaded variant, which has no
  pysam wheel — pysam then compiles from source for several minutes, and the
  result is slower single-threaded. You can spot the wrong build by `cp314t`
  ABI tags in the pip output instead of `cp314`.
- **`samtools`, `bcftools` and `htslib` all pinned to 1.23**, matching the
  htslib that pysam 0.24 bundles. Version skew here causes BGZF and
  EOF-marker problems that surface much later as truncated `.vcf.gz`.

Verify:

```bash
python -c "import sysconfig; print('free-threaded:', bool(sysconfig.get_config_var('Py_GIL_DISABLED')))"
# want: free-threaded: False

python -c "import pysam; print(pysam.__version__, pysam.__samtools_version__)"
# want: 0.24.0 1.23.x

python src/Duopogen.py --help
# want: three subcommands (preProcess, germline, germline-glimpse)
```

If pip printed `Building wheel for pysam`, you are on the free-threaded build.
Remove the env and recreate it.

### 3. Install GLIMPSE2 (optional — only for the `germline-glimpse` backend)

GLIMPSE2 is deliberately **not** in `env.yml`. Its conda package is a C++
binary linked against a specific htslib and boost; co-installing it can either
break the htslib pin above or produce binaries that solve cleanly but fail at
runtime looking for the wrong shared object.

Use the static binaries instead:

```bash
bash install_glimpse2.sh /path/to/persistent/glimpse2
```

This downloads five static binaries (~41 MB) and the b38 genetic maps (~23 MB),
verifies each binary actually links, and prints the two lines you need. Total
footprint ~62 MB. Pass a release tag as the second argument to pin a different
version; the default is `v2.0.1`.

Then in your job scripts, after activating the env:

```bash
export PATH="/path/to/persistent/glimpse2/bin:$PATH"
```

The static binaries are **Linux x86_64 only** — upstream publishes no arm64 or
macOS release assets. On other platforms, build GLIMPSE2 from source or use a
separate conda environment:

```bash
conda create -n glimpse2 -c conda-forge -c bioconda "glimpse-bio>=2"
export PATH="$(conda info --base)/envs/glimpse2/bin:$PATH"
```

Keep it in its own environment so it cannot perturb the pins in
`duopogen_env`.

---

## Reference data

Three things are needed, none of them shipped with this repo.

### Genome reference

GRCh38 **no-alt analysis set**, chr-prefixed, with decoy and HLA contigs.

Not the UCSC browser download (it keeps ALT contigs, which cause ambiguous
multi-mapping and break genotyping at HLA and other immune loci), and not the
bare-numbered NCBI default (contig names would not match CellRanger-ARC's
`chr1`, `chr2`, … BAMs).

Must be indexed:

```bash
samtools faidx hg38.fa       # both backends need the .fai
```

### Imputation panel

The 1000 Genomes 30x NYGC phased panel, one VCF per chromosome, named as
upstream expects:

```
CCDG_14151_B01_GRM_WGS_2020-08-05_<chrom>.filtered.shapeit2-duohmm-phased.vcf.gz
```

Available from the
[EBI 1000 Genomes FTP](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/).
Roughly 520 MB per chromosome.

> **If you are benchmarking against a sample that is itself in the panel**
> — GM12878/NA12878 is, along with both her parents — you must remove the
> trio from the panel first, or the concordance you measure is meaningless.
> See [Choosing a backend](#choosing-a-backend).

### Genetic maps (GLIMPSE2 backend only)

Installed by `install_glimpse2.sh` into `<install_dir>/maps/genetic_maps.b38/`.

The Beagle backend runs without a map and logs
`No genetic map is specified: using 1 cM = 1 Mb`. That is a poor approximation
and one of the reasons the GLIMPSE2 path is expected to be more accurate.

---

## Running the pipeline

### Input

A CSV with one line per sample: `id,/absolute/path/to.bam`

```
DONOR_A,/data/DONOR_A.bam
DONOR_B,/data/DONOR_B.bam
```

BAMs must be coordinate-sorted and indexed (`.bai` or `.csi` — both work).
The **id column becomes the VCF sample name**, so it must be unique. Do not
rely on filenames being distinct.

If a donor has several libraries, merge them first. `samtools merge` requires
sorted inputs and emits sorted output, so no re-sort is needed:

```bash
samtools merge -@ 8 -c -o DONOR_A.bam lib1.bam lib2.bam
samtools index -@ 8 DONOR_A.bam
```

### Step 1 — preProcess

Filters reads by alignment mismatch, splits by chromosome, stamps a read group
derived from the sample id.

```bash
python src/Duopogen.py preProcess \
    -b bam.lst \
    -o /path/to/workdir \
    -t 10 --bam-threads 2
```

`-t` is concurrent per-chromosome jobs; `--bam-threads` is BGZF threads inside
each. Keep the product at or below your core allocation.

**Pool all donors into one `bam.lst`.** The panel-model cost is shared across
whatever is in a single run, so a pooled run is cheaper than N separate ones.

### Step 2a — germline (Beagle 4.1 backend)

```bash
python src/Duopogen.py germline \
    -o /path/to/workdir \
    -g hg38.fa \
    -p /path/to/1kg3_panel/ \
    -r resource/GRCh38.region.lst \
    -s all \
    -t 10
```

Each region spawns a JVM at `-Xmx20g`, so budget ~20 GB per unit of `-t`.
Steps can be run separately with `-s varScan`, `-s varImpute`, `-s varPhasing`.
`-n TRUE` generates the job scripts without running them — useful for
inspecting what will execute, and the fastest way to catch a configuration
error.

`resource/GRCh38.region.lst` accepts whole chromosomes (`chr1`) or chunks
(`chr1,1,50000001`). Chunking increases parallelism and lets Beagle restrict
its panel read to the chunk; whole chromosomes limit you to 22 concurrent jobs.

### Step 2b — germline-glimpse (GLIMPSE2 backend)

```bash
python src/Duopogen.py germline-glimpse \
    -o /path/to/workdir \
    -g hg38.fa \
    -p /path/to/1kg3_panel/ \
    --map-dir /path/to/glimpse2/maps/genetic_maps.b38 \
    --panel-dir /path/to/persistent/glimpse_panel \
    -s all \
    -t 8 --phase-threads 2
```

Steps run in order: `prepPanel` → `varScan` → `phase` → `ligate`.

**`--panel-dir` should be persistent and outside your run directory.** The
binarised reference panel is built once and reused by every donor and every
subsequent run — that reuse is most of the speed advantage. `prepPanel` skips
chromosomes that are already built unless you pass `--force-panel`.

Memory per job is far lower than Beagle: budget ~2–4 GB per unit of `-t`, not
20.

Use `-c chr21,chr22` to restrict to a subset while testing.

---

## Choosing a backend

| | Beagle 4.1 | GLIMPSE2 |
|---|---|---|
| Panel model | rebuilt every window, every run | binarised once, reused |
| Genetic map | none | yes |
| Hom-ref evidence | discarded | kept |
| Output sites | observed only (`impute=false`) | all panel sites, with INFO |
| Scaling in donors | window count grows with donors | sublinear |
| Vendored? | yes, the 4.1 jar | no, installed separately |

Beagle 4.1 is upstream's choice and is kept for reproducibility. It is the only
Beagle that accepts genotype-likelihood input — Beagle 5.x dropped `gl=` — so
the vendored jar cannot simply be upgraded.

**GLIMPSE2 is the better choice for anything beyond a handful of donors**, and
the gap widens with cohort size. The published accuracy comparisons are for
uniform low-coverage WGS rather than ATAC, so if accuracy matters for your
analysis, benchmark before committing:

```bash
# 1. rebuild the panel WITHOUT the sample you are validating against
python src/Duopogen.py germline-glimpse -s prepPanel -c chr20,chr22 \
    --exclude-samples NA12878,NA12891,NA12892 ...

# 2. run both backends on the same input

# 3. score against a truth set, r2 by MAF bin
GLIMPSE2_concordance --input concordance.lst --min-val-dp 8 \
    --min-val-gl 0.9999 --output concordance/out
```

Step 1 is not optional when the validation sample is in the panel. Without it
both tools score near-perfectly by construction and the comparison tells you
nothing.

---

## Output

```
<workdir>/
  Bam/                          # filtered per-donor per-chromosome BAMs
    <chrom>.filter.bam.lst
    <id>_<chrom>.filter.bam
  Script/                       # generated job scripts, kept for inspection
  germline/                     # Beagle backend
    <region>.gl.vcf.gz          # genotype likelihoods
    <region>.gp.vcf.gz          # imputed, unphased
    <region>.phased.vcf.gz      # imputed and phased, with contig headers
    <region>.varScan.log
  glimpse/                      # GLIMPSE2 backend
    <chrom>.gl.vcf.gz
    impute/imputed_<chrom>_*.bcf
    <chrom>.phased.bcf
```

Multi-donor output is one multi-sample VCF per region. Split per donor with:

```bash
bcftools view -s DONOR_A -Oz -o DONOR_A.chr1.vcf.gz germline/chr1.phased.vcf.gz
```

Sanity checks worth running on a first output:

```bash
bcftools view -H <out>.phased.vcf.gz | wc -l                    # nonzero
bcftools view -H <out>.phased.vcf.gz | cut -f10 | cut -d: -f1 | sort | uniq -c
```

Het / hom-alt should land somewhere around 1.3–1.6 for a European sample, and
`0|1` and `1|0` should be roughly balanced. All genotypes should use `|`, not
`/` — a `/` means phasing did not run.

---

## Troubleshooting

**`No VCF records found` from Beagle** — the `.gl.vcf.gz` is header-only.
Check `<region>.varScan.log`, then bisect the `varScan` pipeline stage by
stage. The usual causes are an empty filtered BAM or a `bcftools` filter that
matched nothing.

**`java.lang.StackOverflowError` in `dag.MergeableDag.similar`** — Beagle's
DAG merge is recursive and the default 1 MB JVM thread stack is not enough with
many target samples. Add `-Xss16m` to the java calls in `Duopogen.py`.

**`filtered BAM is empty`** — most often the source BAM carries no `NM`/`nM`
tag. htslib filter expressions treat a missing tag as false, so `[NM] < 3`
would drop every read. The code probes for the tag and skips the filter with a
warning rather than silently discarding everything; run
`samtools calmd -b in.bam ref.fa` first if you want the filter applied.

**`Contig 'chrN' is not defined in the header`** — Beagle emits no `##contig`
lines. The phasing step repairs this with `bcftools reheader --fai`; older
output can be fixed in place the same way.

**`bcftools filter -e 'REF !~ ...'` drops everything** — a regression in
bcftools 1.23: the `!~` regex operator excludes every record and emits
per-record `pass=N [BASE]` debug to stderr. Correct in 1.19. This code uses the
non-regex form `-i 'REF="A" || REF="C" || REF="G" || REF="T"'`, which is
verified byte-identical in output.

**`mpileup: invalid option -- 't'`** — you are on old code. samtools 1.23
removed VCF/BCF output from `samtools mpileup`; the current `varScan` step uses
`bcftools mpileup`. Note that bcftools writes the non-reference symbolic allele
as `<*>` where old samtools wrote `<X>`.

**A step exits 0 but produces nothing** — this was upstream's dominant failure
mode. The fork asserts on nonzero exit and guards against empty output at the
BAM-filter and GL stages, but if you add a stage, add a guard with it.

---

## What changed from upstream

**Removed:** the somatic module, `LDrefinement.R` and the other R scripts, the
orphaned MonoVar tree (11 files reachable only from `monovar.py`, which nothing
imports), and every vendored binary in `apps/` except the Beagle jar. Upstream
shipped samtools 1.2 (2015), bcftools 1.8 (2018) and `libcrypto.so.1.0.0`,
which users had to prepend to `LD_LIBRARY_PATH` — shadowing the environment's
own OpenSSL. `vcftools` and `picard.jar` were asserted to exist but never
invoked.

**Fixed:** an uninitialised variable in the mismatch filter that let a read
inherit the previous read's mismatch count; a `runCMD` that raised `NameError`
on success and swallowed failures; script files that were opened per region but
closed only once; a `zless` dependency absent from most compute nodes; missing
`##contig` headers; hardcoded `range(1, 23)` that silently dropped chrX.

**Changed:** binaries resolve from `PATH` rather than `apps/`; the read filter
runs in samtools rather than a Python loop; `multiprocessing` pins the `fork`
start method, since Python 3.14 changed the POSIX default to `forkserver`;
Beagle's `chrom=` carries the chunk interval so the panel read is restricted.

**Unchanged deliberately:** `modelscale=2`, `niterations=0`, `impute=false` in
the Beagle calls. These affect accuracy, not just speed. Do not touch them
without a truth-set comparison.

---

## Citation

If you use this, cite the original Monopogen paper:

> Dou J, Tan Y, Kock KH, *et al.* Single-nucleotide variant calling in
> single-cell sequencing data with Monopogen. *Nature Biotechnology* 42,
> 803–812 (2024).

And, if you use the GLIMPSE2 backend:

> Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. Imputation of
> low-coverage sequencing data from 150,119 UK Biobank genomes.
> *Nature Genetics* 55, 1088–1090 (2023).