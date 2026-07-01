
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
git clone https://github.com/mahealamuq/kabuki-wgs-pipeline-GI4.git
cd kabuki-wgs-pipeline-GI4/
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
bash kabuki_wgs_pipeline.sh
```

**Options:**

```bash
bash kabuki_wgs_pipeline.sh --threads 8
bash kabuki_wgs_pipeline.sh --workdir /data/my-analysis
bash kabuki_wgs_pipeline.sh --threads 8 --workdir /data/out
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

**GATK4 (Default)**

GATK4 performs local reassembly within **HaplotypeCaller**, so a separate indel realignment step is **not required**.

**GATK3 (Optional)**

If you wish to use the legacy GATK3 indel realignment workflow, specify the location of the GATK3 JAR file in `config/pipeline.conf`:

```bash
GATK3_JAR="/path/to/GenomeAnalysisTK.jar"
```

---

### Step 8 · Alignment QC

`samtools idxstats` — mapped reads per chromosome  
`samtools flagstat` — total reads, duplicate rate, mapping rate, pairing statistics

---

### Steps 9–11 · Variant Calling

To improve confidence in detected variants, the pipeline performs variant calling using **three independent algorithms**. Comparing the results helps identify high-confidence variants that are consistently detected across multiple callers.

| Variant Caller | Method | Description |
|----------------|--------|-------------|
| **bcftools** | Pileup-based | Fast and conservative variant calling using sequence pileups and likelihood models. |
| **FreeBayes** | Haplotype-based | Detects variants using local haplotype reconstruction and offers improved sensitivity for small insertions and deletions (indels). |
| **GATK HaplotypeCaller** | Local de novo assembly | Performs local reassembly of reads around candidate variants and is widely regarded as the gold standard for germline variant discovery. |

> [!IMPORTANT]
> **Ploidy Considerations**
>
> By default, **bcftools** assumes diploid genotypes for all chromosomes. When analysing **chromosome X** in **male samples**, the correct ploidy should be specified because males are **hemizygous** for chromosome X.
>
> Use an appropriate ploidy definition (for example, `--ploidy GRCh37` or a custom ploidy file) to ensure accurate genotype calling.
---

### Step 12 · Variant Caller Concordance

To assess the consistency of variant detection, the pipeline compares the variant sets produced by **bcftools**, **FreeBayes**, and **GATK HaplotypeCaller**.

Agreement between callers is evaluated using the **Jaccard statistic**, which measures the similarity between two sets of variants.

```text
Jaccard Index = |Intersection| / |Union|

0 = No shared variants
1 = Perfect agreement between callers
```

The pipeline also generates a comprehensive three-way comparison of all variant callers, including:

- Total variants detected by each caller
- Shared variants between callers
- Unique variants identified by each caller
- Pairwise Jaccard similarity statistics

The comparison results are written to:

```text
results/patientA/patientA.caller_comparison.txt
results/patientB/patientB.caller_comparison.txt
```

These reports provide an overview of concordance between variant callers and help assess the reliability of detected variants before downstream annotation and filtering.

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

### Step 15 · Candidate Variant Filtering

Following functional annotation, the pipeline prioritises potential disease-causing variants using **SnpSift**. The filtering strategy retains only high-confidence variants that are predicted to have a functional impact and are absent from the dbSNP database.

```bash
cat snpeff.vcf \
| vcfEffOnePerLine.pl \
| SnpSift filter \
"((ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH')) & (na ID) & (QUAL > 30)"
```

The filtering criteria are summarised below:

| Filter | Description | Purpose |
|---------|-------------|---------|
| `ANN[*].IMPACT = HIGH or MODERATE` | Retains variants predicted to alter protein function | Prioritises potentially pathogenic variants |
| `na ID` | Keeps variants without a dbSNP rsID | Removes known common variants and enriches for novel candidates |
| `QUAL > 30` | Keeps high-confidence variant calls | Reduces false-positive variant calls |

The resulting VCF contains a shortlist of candidate variants suitable for downstream interpretation and validation.

**Output**

```text
patientA.candidates.vcf
patientB.candidates.vcf
```

These files represent the final set of prioritised variants for each sample and can be further investigated using databases such as **ClinVar**, **OMIM**, **gnomAD**, or the **Ensembl Variant Effect Predictor (VEP)**.

---

### Step 16 · VEP Submission

The final step prepares a **Variant Call Format (VCF)** file for downstream annotation using the **Ensembl Variant Effect Predictor (VEP)**.

The generated file:

```text
results/vep_submission.vcf
```

can be uploaded directly to the Ensembl VEP web interface for additional annotation and clinical interpretation.

The following resources can be explored through VEP:

| Database | Purpose |
|----------|---------|
| **Ensembl VEP** | Predicts variant consequences and integrates annotations from SIFT, PolyPhen-2, ClinVar, and other resources. |
| **gnomAD** | Provides population allele frequencies to identify rare and common variants. |
| **ExAC** | Offers population frequency data and gene constraint metrics, including pLI scores. |
| **ClinVar** | Contains clinically curated information on variant pathogenicity and disease associations. |
| **OMIM** | Provides detailed information on gene–disease relationships and inherited disorders. |

The VEP output complements the annotations generated by **snpEff** and provides additional evidence for prioritising candidate disease-causing variants.

---

## 📂 Output Files

### Patient A

| File | Description |
|------|-------------|
| `patientA.sorted.bam` | Coordinate-sorted alignment file |
| `patientA.markdups.bam` | BAM file with PCR duplicates marked |
| `patientA.markdups_metrics.txt` | Duplicate marking statistics |
| `patientA.realigned.bam` | Final BAM used for variant calling |
| `patientA.qc_stats.txt` | Alignment quality statistics (`idxstats` and `flagstat`) |
| `patientA.bcftools.vcf` | Variants identified by **bcftools** |
| `patientA.freebayes.vcf` | Variants identified by **FreeBayes** |
| `patientA.gatk.vcf` | Variants identified by **GATK HaplotypeCaller** |
| `patientA.caller_comparison.txt` | Concordance analysis between variant callers |
| `patientA.fb.dbsnp.vcf` | FreeBayes variants annotated with dbSNP identifiers |
| `patientA.fb.snpeff.vcf` | Functionally annotated variants generated by **snpEff** |
| `patientA.candidates.vcf` | Final filtered candidate variants |

### Patient B

| File | Description |
|------|-------------|
| `patientB.sorted.bam` | Coordinate-sorted alignment file |
| `patientB.markdups.bam` | BAM file with PCR duplicates marked |
| `patientB.markdups_metrics.txt` | Duplicate marking statistics |
| `patientB.realigned.bam` | Final BAM used for variant calling |
| `patientB.qc_stats.txt` | Alignment quality statistics (`idxstats` and `flagstat`) |
| `patientB.bcftools.vcf` | Variants identified by **bcftools** |
| `patientB.freebayes.vcf` | Variants identified by **FreeBayes** |
| `patientB.gatk.vcf` | Variants identified by **GATK HaplotypeCaller** |
| `patientB.caller_comparison.txt` | Concordance analysis between variant callers |
| `patientB.fb.dbsnp.vcf` | FreeBayes variants annotated with dbSNP identifiers |
| `patientB.fb.snpeff.vcf` | Functionally annotated variants generated by **snpEff** |
| `patientB.candidates.vcf` | Final filtered candidate variants |

### Shared Output

| File | Description |
|------|-------------|
| `results/vep_submission.vcf` | VCF file prepared for annotation with the **Ensembl Variant Effect Predictor (VEP)** |

---
## ⚙️ Configuration

The pipeline can be customised using command-line options.

### Available Options

| Option | Description | Example |
|--------|-------------|---------|
| `--threads` | Number of CPU threads to use during the analysis | `--threads 8` |
| `--workdir` | Directory where all output files will be stored | `--workdir ~/kabuki-analysis` |

### Examples

Run the pipeline using the default settings:

```bash
bash kabuki_wgs_pipeline.sh
```

Run the pipeline using 8 CPU threads:

```bash
bash kabuki_wgs_pipeline.sh --threads 8
```

Run the pipeline with a custom output directory:

```bash
bash kabuki_wgs_pipeline.sh --workdir ~/kabuki-analysis
```

Run the pipeline using both options:

```bash
bash kabuki_wgs_pipeline.sh \
    --threads 8 \
    --workdir ~/kabuki-analysis
```

## 📊 Pipeline Performance

The pipeline has been successfully tested on Ubuntu Linux using publicly available whole-genome sequencing data.

### Recommended System Requirements

| Component | Requirement |
|-----------|-------------|
| Operating System | Ubuntu 22.04 LTS (or compatible Linux distribution) |
| CPU | Quad-core processor or higher |
| Memory | Minimum 8 GB RAM (16 GB recommended) |
| Disk Space | Approximately 5 GB for analysis outputs (excluding downloaded datasets) |

Pipeline execution time depends on the available hardware, internet speed (during data download), and the number of CPU threads specified using the `--threads` option.

---
## 📌 Notes

- The pipeline is designed for the **GRCh37 (hg19)** human reference genome.
- Analysis focuses on the **KMT2D** and **KDM6A** genomic regions associated with Kabuki syndrome.
- Publicly available **NA12878/HG001** whole-genome sequencing data are used to demonstrate the complete workflow.
- All intermediate and final analysis files are generated automatically within the specified working directory.
- The pipeline is fully automated and reproducible using a single Bash script.

---

## 📚 References

The following publications describe the software tools and datasets used in this project.

| Tool / Resource | Citation |
|-----------------|----------|
| **BWA** | Li, H., & Durbin, R. (2009). *Fast and accurate short read alignment with Burrows-Wheeler Transform.* **Bioinformatics**, 25(14), 1754–1760. |
| **GATK** | McKenna, A., et al. (2010). *The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.* **Genome Research**, 20(9), 1297–1303. |
| **FreeBayes** | Garrison, E., & Marth, G. (2012). *Haplotype-based variant detection from short-read sequencing.* arXiv:1207.3907. |
| **snpEff / SnpSift** | Cingolani, P., et al. (2012). *A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff.* **Fly**, 6(2), 80–92. |
| **BEDTools** | Quinlan, A. R., & Hall, I. M. (2010). *BEDTools: A flexible suite of utilities for comparing genomic features.* **Bioinformatics**, 26(6), 841–842. |
| **Kabuki Syndrome** | Ng, S. B., et al. (2010). *Exome sequencing identifies MLL2 mutations as a cause of Kabuki syndrome.* **Nature Genetics**, 42(9), 790–793. |
| **Genome in a Bottle (NA12878/HG001)** | Zook, J. M., et al. (2014). *Integrating human sequence datasets provides a resource of benchmark SNP and indel genotype calls.* **Nature Biotechnology**, 32(3), 246–251. |

---

## 📄 License

This project is distributed under the **MIT License**.

See the [LICENSE](LICENSE) file for the complete license text.

---

## 👨‍💻 Author

**Mahe Alam**

Bioinformatics | Human Genomics | Variant Analysis | Next-Generation Sequencing (NGS)

GitHub: **https://github.com/mahealamuq**

---

