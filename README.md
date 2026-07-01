
# 🧬 Kabuki WGS Pipeline

### A reproducible whole-genome sequencing pipeline for detecting disease-causing variants associated with Kabuki syndrome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Shell](https://img.shields.io/badge/Shell-Bash-blue.svg)](pipeline.sh)
[![Reference Genome](https://img.shields.io/badge/Reference-GRCh37%20%7C%20hg19-green.svg)](https://genome.ucsc.edu/)
[![Conda](https://img.shields.io/badge/Install-Conda-brightgreen.svg)](environment.yml)

</div>

---

## Table of Contents

* [Background](#-background)
* [Pipeline Overview](#-pipeline-overview)
* [Requirements](#-requirements)
* [Installation](#-installation)
* [Usage](#-usage)
* [Pipeline Steps](#-pipeline-steps)
* [Output Files](#-output-files)
* [Configuration](#-configuration)
* [Troubleshooting](#-troubleshooting)
* [References](#-references)

---

## 🔬 Background

**Kabuki syndrome** is a rare congenital disorder (\~1 in 32,000 births) characterised by intellectual disability, distinctive facial features, and post-natal growth deficiency. It is caused by mutations in one of two chromatin-remodelling genes:

||*KMT2D*|*KDM6A*|
|-|:-:|:-:|
|**Syndrome**|Kabuki type 1|Kabuki type 2|
|**Chromosome**|12q13.12|Xp11.3|
|**Protein**|H3K4 methyltransferase|H3K27 demethylase|
|**Inheritance**|Autosomal dominant|X-linked dominant|
|**OMIM**|[#147920](https://www.omim.org/entry/147920)|[#300867](https://www.omim.org/entry/300867)|

Most causative variants are **de novo** — absent from population databases — which is the key insight behind the filtering strategy: after calling all variants, those with a dbSNP rsID are removed, leaving only novel candidates.

---

## 🔄 Pipeline Overview

```
FASTQ reads (NA12878 / HG001)
        │
        ▼
 ┌─────────────┐
 │  BWA MEM   │  ← align to hg19 (chr12 + chrX)
 └──────┬──────┘
        │ SAM
        ▼
 ┌─────────────┐
 │  samtools  │  ← SAM → BAM, coordinate sort, index
 └──────┬──────┘
        │ sorted BAM
        ▼
 ┌──────────────────┐
 │ Picard           │  ← flag PCR duplicates
 │ MarkDuplicates   │
 └────────┬─────────┘
          │ deduplicated BAM
          ▼
 ┌──────────────────────────────────────────────┐
 │               Variant Calling                │
 │                                              │
 │  bcftools mpileup   freebayes   GATK HC      │
 └──────────────────────────────────────────────┘
          │ VCF × 3
          ▼
 ┌────────────────────────┐
 │  bedtools              │  ← Jaccard + intersection
 │  jaccard / intersect   │
 └────────────┬───────────┘
              │
              ▼
 ┌────────────────────────┐
 │  SnpSift annotate      │  ← tag known variants (dbSNP rsIDs)
 └────────────┬───────────┘
              │
              ▼
 ┌────────────────────────┐
 │  snpEff GRCh37.75      │  ← predict functional consequence
 └────────────┬───────────┘
              │
              ▼
 ┌──────────────────────────────────────────────────┐
 │  SnpSift filter                                   │
 │  HIGH or MODERATE impact  +  no rsID  +  QUAL>30  │
 └──────────────────────────────────────────────────┘
              │
              ▼
    candidates.vcf  +  vep\\\_submission.vcf
```

---

## 📦 Requirements

|Tool|Version|Purpose|
|-|:-:|-|
|[BWA](https://github.com/lh3/bwa)|0.7.17|Short-read alignment|
|[samtools](https://www.htslib.org)|1.17|BAM processing and QC|
|[bcftools](https://www.htslib.org)|1.17|Pileup-based variant calling|
|[freebayes](https://github.com/freebayes/freebayes)|1.3.6|Haplotype-based variant calling|
|[GATK](https://gatk.broadinstitute.org)|4.4|HaplotypeCaller variant calling|
|[Picard](https://broadinstitute.github.io/picard)|3.0|PCR duplicate marking|
|[bedtools](https://bedtools.readthedocs.io)|2.31|Variant set comparison|
|[snpEff + SnpSift](https://pcingola.github.io/SnpEff)|5.1 / 5.3|Annotation and filtering|
|[htslib](https://www.htslib.org)|1.17|tabix + bgzip|

> **Java Note:**  
> GATK 4 requires **Java 17**, while **SnpSift 5.3+** requires **Java 21**.  
> The pipeline automatically detects and uses a compatible Java installation for each tool, allowing both applications to run correctly even when the Conda environment uses Java 17.

**System Requirements**

- **Operating System:** Linux
- **Memory:** Minimum 8 GB RAM (16 GB recommended)
- **Disk Space:** Approximately 5 GB (excluding downloaded datasets)

---

## ⚡ Installation

### 1\. Clone

```bash
git clone https://github.com/YOUR\\\_USERNAME/kabuki-wgs-pipeline.git
cd kabuki-wgs-pipeline
```

### 2\. Create conda environment

```bash
conda env create -f environment.yml
conda activate kabuki-wgs
```

### 3\. Download snpEff database (one-time, \~1.5 GB)

```bash
snpEff download GRCh37.75
```

---

## 🚀 Usage

```bash
conda activate kabuki-wgs
bash pipeline.sh
```

**Options:**

```bash
bash pipeline.sh --threads 8
bash pipeline.sh --workdir /data/my-analysis
bash pipeline.sh --threads 8 --workdir /data/out
```

**Bring your own reads** — place FASTQs here and the download step is skipped:

```
kabuki-wgs-output/data/
├── patientA\\\_1.fq.gz
├── patientA\\\_2.fq.gz
├── patientB\\\_1.fq.gz
└── patientB\\\_2.fq.gz
```

---

## 🔍 Pipeline Steps

### Step 0 · Dependency Check + Java Detection

Checks all required tools are on `PATH`. Then independently locates a compatible Java binary for both SnpSift (requires Java 21) and snpEff (may also require Java 21), separate from the default `java` used by GATK4 (requires Java 17). This solves the common conflict where a conda environment can only hold one Java version at a time.

---

### Step 1 · Directory Setup

```
kabuki-wgs-output/
├── data/           ← FASTQ reads
├── ref/
│   ├── genome/     ← hg19 FASTA + indexes
│   ├── dbsnp/      ← dbSNP VCF (bgzipped + tabixed)
│   └── intervals/  ← KMT2D + KDM6A BED file
├── results/
│   ├── patientA/
│   └── patientB/
└── logs/
```

---

### Step 2 · Reference Genome

Downloads **chromosome 12 (chr12)** and **chromosome X (chrX)** from the **UCSC hg19 (GRCh37)** reference genome, merges the two chromosomes into a single reference FASTA file, and builds the required indexes for **BWA**, **samtools**, and **Picard**.

> [!TIP]
> Using only **chr12** and **chrX** reduces the reference genome size from approximately **3 GB** to **550 MB** without affecting the analysis, as the pipeline focuses exclusively on the **KMT2D** and **KDM6A** genes located on these chromosomes.

---

### Step 3 · Input Reads

Streams NA12878 (HG001) reads for only the KMT2D and KDM6A regions (\~50–100 MB vs the full \~100 GB BAM) using two fallback sources.

**Target regions (hg19):**

|Gene|Region|Syndrome|
|-|-|-|
|*KMT2D*|chr12:49,400,000–49,500,000|Kabuki type 1|
|*KDM6A*|chrX:44,900,000–45,100,000|Kabuki type 2|

---

### Step 4 · Alignment — BWA MEM

```bash
bwa mem \
    -t 4 \
    -R "@RG\tID:patientA\tSM:patientA\tPL:ILLUMINA\tLB:lib1" \
    hg19_chr12_chrX.fasta \
    patientA_1.fq.gz \
    patientA_2.fq.gz \
    > patientA.sam
```

---

### Step 5 · SAM → Sorted BAM

Converts, coordinate-sorts, and indexes the alignment. The original SAM is deleted to save disk space.

---

### Step 6 · PCR Duplicate Marking

Picard MarkDuplicates flags reads that are PCR duplicates (identical read pairs from the same original DNA molecule) so variant callers can exclude them from allele frequency calculations.

---

### Step 7 · Indel Realignment

**GATK4 (default):** handled internally by HaplotypeCaller — no separate step needed.  
**GATK3 (optional):** set `GATK3\\\_JAR` in `config/pipeline.conf` and uncomment the GATK3 block.

---

### Step 8 · Alignment QC

`samtools idxstats` — mapped reads per chromosome  
`samtools flagstat` — total reads, duplicate rate, mapping rate, pairing statistics

---

### Steps 9–11 · Variant Calling (three callers)

Three independent callers run in parallel for concordance comparison:

|Caller|Approach|Notes|
|-|-|-|
|**bcftools**|Pileup + multinomial likelihood|Fast, conservative|
|**freebayes**|Haplotype-based|Better indel sensitivity|
|**GATK HaplotypeCaller**|Local de novo assembly|Gold standard for clinical WGS|

> ⚠️ \\\*\\\*Ploidy:\\\*\\\* bcftools defaults to diploid. Use `--ploidy GRCh37` for correct chrX calls in biological males (hemizygous).

---

### Step 12 · Caller Concordance

Measures agreement between callers using the **Jaccard statistic**:

```
Jaccard = |intersection| / |union|   (0 = no overlap, 1 = perfect agreement)
```

A full three-way intersection table is written to `caller\\\_comparison.txt`.

---

### Step 13 · dbSNP Annotation

Tags each variant with a dbSNP rsID if it appears in the population database. Two methods are tried:

1. Remote tabix query against NCBI (using correct RefSeq contig names `NC\\\_000012.11`, `NC\\\_000023.10`)
2. Ensembl REST API fallback

Variants with rsIDs are known in the general population and unlikely to be the sole cause of a rare de novo disease.

---

### Step 14 · Functional Annotation — snpEff

Predicts the consequence of each variant on every overlapping transcript:

|Impact|Examples|Effect|
|-|-|-|
|🔴 **HIGH**|Stop gained, frameshift, splice site|Likely loss of function|
|🟡 **MODERATE**|Missense, in-frame indel|Possible functional change|
|🟢 **LOW**|Synonymous|Unlikely to affect protein|
|⚪ **MODIFIER**|Intronic, UTR|No direct protein effect|

---

### Step 15 · Candidate Filtering

```bash
cat snpeff.vcf | vcfEffOnePerLine.pl \\\\
| SnpSift filter \\\\
    "((ANN\\\[\\\*].IMPACT has 'MODERATE') | (ANN\\\[\\\*].IMPACT has 'HIGH')) \\\& (na ID) \\\& (QUAL > 30)"
```

|Filter|Rationale|
|-|-|
|`IMPACT HIGH or MODERATE`|Only protein-altering variants|
|`na ID`|No rsID — novel, absent from population databases|
|`QUAL > 30`|≥ 99.9% confidence (Phred scale)|

---

### Step 16 · VEP Submission

Generates `results/vep\\\_submission.vcf` for cross-referencing:

|Database|Purpose|
|-|-|
|[Ensembl VEP](https://grch37.ensembl.org/Homo_sapiens/Tools/VEP)|SIFT, PolyPhen-2, ClinVar|
|[gnomAD](https://gnomad.broadinstitute.org)|Population allele frequency|
|[ExAC](http://exac.broadinstitute.org)|LoF intolerance (pLI score)|
|[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar)|Clinical significance|
|[OMIM](https://www.omim.org)|Disease–gene association|

---

## 📂 Output Files

```
results/patientA/
├── patientA.sorted.bam               ← coordinate-sorted alignment
├── patientA.markdups.bam             ← PCR duplicates flagged
├── patientA.markdups\\\_metrics.txt     ← duplicate rate report
├── patientA.realigned.bam            ← final BAM for variant calling
├── patientA.qc\\\_stats.txt             ← idxstats + flagstat
├── patientA.bcftools.vcf             ← bcftools calls
├── patientA.freebayes.vcf            ← freebayes calls
├── patientA.gatk.vcf                 ← GATK HaplotypeCaller calls
├── patientA.caller\\\_comparison.txt    ← Jaccard + intersection stats
├── patientA.fb.dbsnp.vcf             ← freebayes + dbSNP rsIDs
├── patientA.fb.snpeff.vcf            ← functionally annotated VCF
└── patientA.candidates.vcf           ← filtered candidates ✅

results/vep\\\_submission.vcf            ← paste into Ensembl VEP
```

---

## ⚙️ Configuration

Edit `config/pipeline.conf`:

```bash
THREADS=8                              # CPU threads (default: 4)
WORKDIR="/data/my-analysis"            # output directory

# Use a local dbSNP file (skips download)
# DBSNP\\\_VCF="/data/ref/dbsnp\\\_b151\\\_GRCh37.vcf.gz"

# Enable GATK3 indel realignment
# GATK3\\\_JAR="/opt/GATK3/GenomeAnalysisTK.jar"
```

---

## 🛠️ Troubleshooting

<details>
<summary><strong>SnpSift or snpEff fails with UnsupportedClassVersionError</strong></summary>

SnpSift 5.3+ and snpEff 5.1+ require Java 21. The pipeline auto-detects a compatible Java under `/usr/lib/jvm/`. If detection fails, install it:

```bash
sudo apt install openjdk-21-jre-headless
```

The pipeline will find it automatically on the next run.

</details>

<details>
<summary><strong>Remote BAM download hangs at Step 3</strong></summary>

Use the SRA fallback:

```bash
prefetch SRR622461
fastq-dump --split-files --gzip SRR622461
mv SRR622461\\\_1.fastq.gz \\\~/kabuki-wgs-output/data/patientA\\\_1.fq.gz
mv SRR622461\\\_2.fastq.gz \\\~/kabuki-wgs-output/data/patientA\\\_2.fq.gz
```

</details>

<details>
<summary><strong>dbSNP annotation produces a 0-byte or header-only file</strong></summary>

NCBI's dbSNP VCF uses RefSeq contig names (`NC\\\_000012.11`, `NC\\\_000023.10`), not plain `12`/`X`. Querying with bare chromosome numbers silently returns zero records. The pipeline now uses the correct contig names automatically. If both methods fail, delete the slice and re-run:

```bash
rm \\\~/kabuki-wgs-output/ref/dbsnp/dbsnp\\\_b151\\\_GRCh37.vcf.gz\\\*
bash pipeline.sh
```

</details>

<details>
<summary><strong>bedtools jaccard returns empty or errors</strong></summary>

bedtools jaccard requires coordinate-sorted input VCFs:

```bash
bcftools sort patientA.bcftools.vcf -o patientA.bcftools.sorted.vcf
```

</details>

<details>
<summary><strong>snpEff database not found</strong></summary>

```bash
snpEff download GRCh37.75
```

</details>

---

## 📚 References

|Tool|Citation|
|-|-|
|BWA|Li H \& Durbin R (2009). *Bioinformatics* 25(14):1754–1760|
|GATK|McKenna A et al. (2010). *Genome Research* 20:1297–1303|
|freebayes|Garrison E \& Marth G (2012). arXiv:1207.3907|
|snpEff/SnpSift|Cingolani P et al. (2012). *Fly* 6(2):80–92|
|bedtools|Quinlan AR \& Hall IM (2010). *Bioinformatics* 26(6):841–842|
|Kabuki syndrome|Ng SB et al. (2010). *Nature Genetics* 42:790–793|
|NA12878 / GIAB|Zook JM et al. (2014). *Nature Biotechnology* 32:246–251|

---

## 📄 License

MIT — see [LICENSE](LICENSE)

