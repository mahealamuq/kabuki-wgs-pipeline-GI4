
# kabuki-wgs-pipeline

A reproducible whole-genome sequencing (WGS) pipeline for detecting disease-causing variants in **KMT2D** and **KDM6A** — the two genes responsible for Kabuki syndrome.

The pipeline runs from raw FASTQ reads through alignment, duplicate marking, variant calling with three independent callers, functional annotation, and candidate filtering down to a shortlist of high-confidence novel variants.

---

## Table of Contents

- [Background](#background)
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## Background

**Kabuki syndrome** is a rare congenital disorder characterised by intellectual disability, post-natal growth deficiency, and distinct facial features. It is caused by mutations in one of two genes:

| | KMT2D | KDM6A |
|---|---|---|
| **Syndrome** | Kabuki syndrome type 1 | Kabuki syndrome type 2 |
| **Chromosome** | 12q13.12 | Xp11.3 |
| **OMIM** | [#147920](https://www.omim.org/entry/147920) | [#300867](https://www.omim.org/entry/300867) |
| **Protein** | H3K4 methyltransferase | H3K27 demethylase |
| **Inheritance** | Autosomal dominant | X-linked dominant |
| **Prevalence** | ~1 in 32,000 | Less common |

Both genes encode chromatin remodelling enzymes. Most causative variants are **de novo** (not inherited from either parent) and are absent from population databases such as gnomAD and ExAC — which is why filtering for novel variants is a powerful enrichment strategy.

---

## Pipeline Overview

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
 │ Picard           │  ← mark PCR duplicates
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
 ┌─────────────────────────────┐
 │  bedtools jaccard           │  ← caller concordance
 │  bedtools intersect         │
 └──────────────┬──────────────┘
                │
                ▼
 ┌──────────────────────┐
 │  SnpSift annotate    │  ← tag known variants with dbSNP rsIDs
 └──────────┬───────────┘
            │
            ▼
 ┌──────────────────────┐
 │  snpEff              │  ← predict functional consequence per transcript
 └──────────┬───────────┘
            │
            ▼
 ┌──────────────────────────────────────────────────┐
 │  SnpSift filter                                   │
 │  HIGH or MODERATE impact  +  no rsID  +  QUAL>30  │
 └──────────────────────────────────────────────────┘
            │
            ▼
    patientA.candidates.vcf
            │
            ▼
    vep_submission.vcf  →  Ensembl VEP / gnomAD / ClinVar
```

---

## Requirements

| Tool | Version | Role |
|---|---|---|
| [BWA](https://github.com/lh3/bwa) | ≥ 0.7.17 | Short-read alignment |
| [samtools](https://www.htslib.org) | ≥ 1.17 | BAM processing and QC |
| [bcftools](https://www.htslib.org) | ≥ 1.17 | Pileup-based variant calling |
| [freebayes](https://github.com/freebayes/freebayes) | ≥ 1.3.6 | Haplotype-based variant calling |
| [GATK](https://gatk.broadinstitute.org) | ≥ 4.4 | HaplotypeCaller variant calling |
| [Picard](https://broadinstitute.github.io/picard) | ≥ 3.0 | PCR duplicate marking |
| [bedtools](https://bedtools.readthedocs.io) | ≥ 2.31 | Variant set comparison |
| [snpEff / SnpSift](https://pcingola.github.io/SnpEff) | ≥ 5.1 | Functional annotation + filtering |
| [htslib](https://www.htslib.org) | ≥ 1.17 | tabix / bgzip |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/kabuki-wgs-pipeline.git
cd kabuki-wgs-pipeline
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate kabuki-wgs
```

### 3. Download the snpEff database (one-time, ~1.5 GB)

```bash
snpEff download GRCh37.75
```

---

## Usage

```bash
conda activate kabuki-wgs
bash pipeline.sh
```

**Options:**

```bash
bash pipeline.sh --threads 8                       # use 8 CPU threads
bash pipeline.sh --workdir /data/kabuki-analysis   # custom output directory
bash pipeline.sh --threads 8 --workdir /data/out   # both
```

**Bring your own reads:**

Place paired FASTQ files here before running and the automatic download is skipped:

```
kabuki-wgs-output/
└── data/
    ├── patientA_1.fq.gz
    ├── patientA_2.fq.gz
    ├── patientB_1.fq.gz
    └── patientB_2.fq.gz
```

---

## Pipeline Steps

### Step 0 — Dependency check

Verifies all required tools are on `PATH`. If `SnpSift` is missing but `snpEff` is installed, the script auto-detects and adds the snpEff bin directory to `PATH` (they ship together in the conda package).

---

### Step 1 — Directory structure

Creates the working directory layout:

```
kabuki-wgs-output/
├── data/              raw FASTQ reads
├── ref/
│   ├── genome/        reference FASTA + BWA/samtools indexes
│   ├── dbsnp/         dbSNP VCF (bgzipped + tabix indexed)
│   └── intervals/     target regions BED file
├── results/
│   ├── patientA/
│   └── patientB/
└── logs/
```

---

### Step 2 — Reference genome

Downloads chr12 and chrX from [UCSC hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/), merges them into a single FASTA, and builds all required indexes:

```bash
wget hg19/chr12.fa.gz
wget hg19/chrX.fa.gz
zcat chr12.fa.gz chrX.fa.gz > hg19_chr12_chrX.fasta

bwa index hg19_chr12_chrX.fasta       # BWA index (~4 min)
samtools faidx hg19_chr12_chrX.fasta  # FASTA index
picard CreateSequenceDictionary        # sequence dictionary for GATK
```

> **Why only chr12 and chrX?** *KMT2D* sits on chromosome 12 and *KDM6A* on the X chromosome. Restricting the reference to these two chromosomes keeps file sizes manageable (~550 MB vs ~3 GB for the full genome) while covering all target regions.

---

### Step 3 — Input reads

Streams NA12878 (HG001) reads from the [1000 Genomes Project](https://www.internationalgenome.org/) — only reads overlapping the two target gene regions are downloaded (~50–100 MB instead of the full ~100 GB BAM).

**Target regions (hg19):**

| Gene | Region | Syndrome |
|---|---|---|
| *KMT2D* | chr12:49,400,000–49,500,000 | Kabuki type 1 |
| *KDM6A* | chrX:44,900,000–45,100,000 | Kabuki type 2 |

The script writes the slice to a local temp BAM first, name-sorts it, then converts to paired FASTQ — this three-step approach avoids a known failure where `samtools sort -n` cannot seek back through a remote HTTP stream.

**Fallback if download fails:**

```bash
# Option A — SRA toolkit
prefetch SRR622461
fastq-dump --split-files --gzip SRR622461
mv SRR622461_1.fastq.gz data/patientA_1.fq.gz
mv SRR622461_2.fastq.gz data/patientA_2.fq.gz

# Option B — wget full BAM (~30 GB) then slice locally
wget https://ftp.sra.ebi.ac.uk/vol1/ERR/ERR194/ERR194147/ERR194147.bam
samtools view -b ERR194147.bam 12:49400000-49500000 X:44900000-45100000 -o slice.bam
```

---

### Step 4 — Alignment (BWA MEM)

Aligns paired-end reads to the reference using BWA MEM with a read group header required by downstream GATK tools:

```bash
bwa mem \
  -t 4 \
  -R "@RG\tID:patientA\tSM:patientA\tPL:ILLUMINA\tLB:lib1" \
  hg19_chr12_chrX.fasta \
  patientA_1.fq.gz \
  patientA_2.fq.gz \
  > patientA.sam
```

> **Why BWA MEM?** BWA MEM handles split reads and local alignments better than BWA backtrack and is the recommended aligner for Illumina reads ≥ 70 bp in modern germline sequencing workflows.

---

### Step 5 — SAM → sorted BAM

Converts the SAM alignment to a coordinate-sorted, indexed BAM. Coordinate sorting is required for duplicate marking and all downstream variant callers.

```bash
samtools view -bS patientA.sam \
| samtools sort -@ 4 -o patientA.sorted.bam

samtools index patientA.sorted.bam
```

---

### Step 6 — PCR duplicate marking (Picard MarkDuplicates)

Identifies and flags reads that are PCR duplicates — multiple reads originating from the same original DNA molecule rather than independent sequencing events. These are marked (not removed) so variant callers can exclude them from allele frequency calculations.

```bash
picard MarkDuplicates \
  I=patientA.sorted.bam \
  O=patientA.markdups.bam \
  M=patientA.markdups_metrics.txt
```

The metrics file reports the duplicate rate. A very high duplicate rate (> 50%) suggests the input library had low complexity.

---

### Step 7 — Indel realignment

**GATK4 (default):** Local indel realignment is performed internally by `HaplotypeCaller`, so no separate step is needed.

**GATK3 (optional):** Set `GATK3_JAR` in `config/pipeline.conf` and uncomment the GATK3 block in the script to run `RealignerTargetCreator` + `IndelRealigner` as a separate step before calling.

---

### Step 8 — Alignment QC

Runs two samtools QC tools and saves results to `patientA.qc_stats.txt`:

**`samtools idxstats`** — reads mapped per chromosome:

```
chr12   133851895   12543   0
chrX    155270560    8901   0
```

Columns: `chromosome | length | mapped reads | unmapped mates`

**`samtools flagstat`** — overall mapping statistics including total reads, duplicate count, mapping rate, and paired-end concordance.

---

### Step 9 — Variant calling: bcftools + freebayes

Two independent variant callers are run so their results can be compared for concordance.

**bcftools mpileup | call** — pileup-based, fast and conservative:

```bash
bcftools mpileup -f reference.fasta patientA.realigned.bam \
| bcftools call -vc \
  > patientA.bcftools.vcf
```

> **Ploidy note:** bcftools assumes diploid (2 alleles/site) by default. chrX in biological males is hemizygous (1 copy). Use `--ploidy GRCh37` or a ploidy file for accurate chrX calls.

**freebayes** — haplotype-based, more sensitive for indels and nearby variants:

```bash
freebayes \
  -F 0.2 \
  --min-repeat-entropy 0 \
  -f reference.fasta \
  patientA.realigned.bam \
  > patientA.freebayes.vcf
```

---

### Step 10 — Variant calling: GATK HaplotypeCaller

HaplotypeCaller performs local de novo assembly in active regions before genotyping — the most accurate germline variant caller for clinical WGS.

```bash
gatk HaplotypeCaller \
  -R reference.fasta \
  -I patientA.realigned.bam \
  -L ref/intervals/kabuki_regions.bed \
  -O patientA.gatk.vcf
```

The `-L` flag restricts calling to the KMT2D and KDM6A intervals, significantly reducing runtime.

---

### Step 11 — Caller comparison (bedtools)

Measures concordance between callers using the **Jaccard statistic**:

```
Jaccard = |intersection| / |union|    (0 = no overlap, 1 = perfect agreement)
```

```bash
bedtools jaccard \
  -a patientA.bcftools.vcf \
  -b patientA.freebayes.vcf
```

A three-way intersection table is also produced. Variants called by multiple independent callers are more likely to be true positives.

---

### Step 12 — dbSNP annotation (SnpSift)

Tags each variant with a dbSNP rsID if it appears in the population database. Variants with an rsID are seen in the general population and are unlikely to be the sole cause of a rare Mendelian disease.

```bash
SnpSift annotate \
  ref/dbsnp/dbsnp_b151_GRCh37.vcf.gz \
  patientA.freebayes.vcf \
  > patientA.fb.dbsnp.vcf
```

---

### Step 13 — Functional annotation (snpEff)

Predicts the functional consequence of each variant on every overlapping transcript using Ensembl GRCh37.75. Each variant-transcript pair receives an **impact** label:

| Impact | Consequence types | Likely effect |
|---|---|---|
| **HIGH** | Stop gained, frameshift, splice site | Loss of function |
| **MODERATE** | Missense, in-frame indel | Altered protein function |
| **LOW** | Synonymous, splice region | Unlikely to affect protein |
| **MODIFIER** | Intronic, UTR, intergenic | No direct protein change |

```bash
snpEff -v GRCh37.75 patientA.fb.dbsnp.vcf > patientA.fb.snpeff.vcf
```

> A single variant overlapping a gene with 10 transcripts will generate 10 `ANN` entries in the output — one per transcript. This is why the candidates VCF may have multiple lines per genomic position.

---

### Step 14 — Candidate variant filtering (SnpSift)

Applies three simultaneous filters to reduce thousands of variants to a shortlist of high-confidence candidates:

```bash
cat patientA.fb.snpeff.vcf \
| vcfEffOnePerLine.pl \
| SnpSift filter \
    "((ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH')) \
     & (na ID) \
     & (QUAL > 30)" \
  > patientA.candidates.vcf
```

| Filter | Rationale |
|---|---|
| `IMPACT has 'HIGH'` or `'MODERATE'` | Keeps only protein-altering variants |
| `na ID` | No dbSNP rsID — novel, not seen in population databases |
| `QUAL > 30` | ≥ 99.9% probability the variant is real (Phred scale) |

To count unique variant positions (not annotation lines):

```bash
grep -v '^#' patientA.candidates.vcf | cut -f1,2 | sort -u | wc -l
```

---

### Step 15 — VEP submission

Writes `vep_submission.vcf` containing candidate variants ready to paste into the [Ensembl VEP](https://grch37.ensembl.org/Homo_sapiens/Tools/VEP) web tool for SIFT, PolyPhen-2, and ClinVar annotations.

**Cross-reference databases:**

| Database | URL | Purpose |
|---|---|---|
| Ensembl VEP | https://grch37.ensembl.org/Homo_sapiens/Tools/VEP | Full consequence annotation |
| gnomAD | https://gnomad.broadinstitute.org | Population allele frequency |
| ExAC | http://exac.broadinstitute.org | Exome frequency + LoF intolerance (pLI) |
| ClinVar | https://www.ncbi.nlm.nih.gov/clinvar | Clinical significance |
| OMIM | https://www.omim.org | Disease–gene association |

---

## Output Files

| File | Description |
|---|---|
| `patientA.sorted.bam` | Coordinate-sorted alignment |
| `patientA.markdups.bam` | PCR duplicates flagged |
| `patientA.markdups_metrics.txt` | Duplicate rate statistics |
| `patientA.realigned.bam` | Final BAM for variant calling |
| `patientA.qc_stats.txt` | idxstats + flagstat QC report |
| `patientA.bcftools.vcf` | bcftools variant calls |
| `patientA.freebayes.vcf` | freebayes variant calls |
| `patientA.gatk.vcf` | GATK HaplotypeCaller calls |
| `patientA.caller_comparison.txt` | Jaccard + intersection statistics |
| `patientA.fb.dbsnp.vcf` | freebayes VCF + dbSNP rsIDs |
| `patientA.fb.snpeff.vcf` | Functionally annotated VCF |
| `patientA.candidates.vcf` | Final filtered candidate variants |
| `results/vep_submission.vcf` | Ready for Ensembl VEP web input |

---

## Configuration

Edit `config/pipeline.conf` to override defaults:

```bash
THREADS=8                           # CPU threads (default: 4)
WORKDIR="/data/my-analysis"         # output directory

# Use a local dbSNP VCF to skip the download
# DBSNP_VCF="/data/ref/dbsnp_b151_GRCh37.vcf.gz"

# Enable GATK3 indel realignment
# GATK3_JAR="/opt/GATK3/GenomeAnalysisTK.jar"
```

---

## Troubleshooting

**`samtools sort: failed to read header from "-"`**
Caused by piping a remote BAM directly into `samtools sort`. The pipeline writes to a local temp BAM first — ensure you have the latest version of the script.

**`SnpSift: command not found`**
SnpSift ships inside the snpEff conda package but is sometimes not symlinked on PATH. The script auto-detects this. If it still fails:
```bash
conda install -c bioconda snpsift
```

**`bedtools jaccard` fails or returns empty**
bedtools jaccard requires coordinate-sorted VCFs:
```bash
bcftools sort patientA.bcftools.vcf -o patientA.bcftools.sorted.vcf
```

**snpEff database error**
```bash
snpEff download GRCh37.75
```

**Remote BAM download hangs at Step 3**
The remote source may be temporarily unavailable. Use the SRA fallback:
```bash
prefetch SRR622461
fastq-dump --split-files --gzip SRR622461
mv SRR622461_1.fastq.gz ~/kabuki-wgs-output/data/patientA_1.fq.gz
mv SRR622461_2.fastq.gz ~/kabuki-wgs-output/data/patientA_2.fq.gz
```

---

## References

- Li H & Durbin R (2009). Fast and accurate short read alignment with Burrows-Wheeler Aligner. *Bioinformatics* 25(14):1754–1760
- McKenna A et al. (2010). The Genome Analysis Toolkit. *Genome Research* 20:1297–1303
- Garrison E & Marth G (2012). Haplotype-based variant detection from short-read sequencing. arXiv:1207.3907
- Cingolani P et al. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly* 6(2):80–92
- Quinlan AR & Hall IM (2010). BEDTools. *Bioinformatics* 26(6):841–842
- Ng SB et al. (2010). Exome sequencing identifies MLL2 mutations as a cause of Kabuki syndrome. *Nature Genetics* 42:790–793
- Zook JM et al. (2014). Integrating human sequence data sets provides a resource of benchmark SNP and indel genotype calls. *Nature Biotechnology* 32:246–251

---

## License

MIT — see [LICENSE](LICENSE)
