
# Kabuki Syndrome WGS Variant Analysis Pipeline

![Bioinformatics](https://img.shields.io/badge/Workflow-WGS%20Variant%20Analysis-blue)
![Genome Build](https://img.shields.io/badge/Genome-GRCh37%2Fhg19-green)
![Tools](https://img.shields.io/badge/Tools-BWA%20%7C%20samtools%20%7C%20GATK%20%7C%20snpEff-orange)
![License](https://img.shields.io/badge/Use-Educational-lightgrey)

## Overview

This repository contains a complete whole-genome sequencing variant analysis pipeline for investigating variants in **Kabuki syndrome-associated genes**.

The analysis focuses on two important genes:

| Gene | Associated Condition | Chromosome | Region Used in Pipeline |
|---|---|---|---|
| **KMT2D** | Kabuki syndrome type 1 | Chromosome 12 | `12:49400000-49500000` |
| **KDM6A** | Kabuki syndrome type 2 | Chromosome X | `X:44900000-45100000` |

The pipeline starts from public sequencing data, prepares the reference genome, aligns reads, processes BAM files, calls variants using multiple variant callers, annotates variants, filters candidate variants, and creates a VCF file for Ensembl VEP submission.

## Project Aim

The aim of this project is to demonstrate a reproducible human genomics workflow for identifying and annotating variants in disease-associated genes using WGS data.

This project shows how raw sequencing data can be converted into biologically meaningful variant results using standard command-line bioinformatics tools.

## Dataset

The pipeline uses publicly available **NA12878 / HG001 whole-genome sequencing data**.

Instead of downloading the entire genome dataset, the script extracts only reads that map to the KMT2D and KDM6A target regions. This keeps the project lightweight and suitable for a local Ubuntu machine.

| Data Type | Description |
|---|---|
| Sample | NA12878 / HG001 |
| Sequencing Type | Whole-genome sequencing |
| Reference Genome | GRCh37 / hg19 |
| Reference Chromosomes | chr12 and chrX |
| Target Genes | KMT2D and KDM6A |

## Pipeline Workflow

```text
1.  Create project directories
2.  Download hg19 chr12 and chrX reference genome
3.  Index the reference genome
4.  Download target-region NA12878 reads
5.  Convert extracted reads to paired FASTQ
6.  Align reads to reference genome using BWA MEM
7.  Convert SAM to BAM
8.  Sort and index BAM files
9.  Mark PCR duplicates using Picard
10. Prepare final BAM files for variant calling
11. Run alignment QC using samtools
12. Call variants using bcftools
13. Call variants using freebayes
14. Call variants using GATK HaplotypeCaller
15. Compare variant callers using bedtools
16. Annotate variants with dbSNP using SnpSift
17. Annotate functional effects using snpEff
18. Filter candidate disease-relevant variants
19. Create VEP submission file
```

## Workflow Diagram

```text
Public WGS BAM
     ↓
Target-region read extraction
     ↓
Paired FASTQ files
     ↓
BWA MEM alignment
     ↓
SAM file
     ↓
Sorted and indexed BAM
     ↓
Duplicate-marked BAM
     ↓
Alignment QC
     ↓
Variant calling
 ┌────────────┬────────────┬────────────┐
 │ bcftools   │ freebayes  │ GATK       │
 └────────────┴────────────┴────────────┘
     ↓
Caller comparison
     ↓
dbSNP annotation
     ↓
snpEff functional annotation
     ↓
Candidate variant filtering
     ↓
VEP-ready VCF file
```

## Tools Used

| Tool | Purpose |
|---|---|
| **BWA MEM** | Aligns sequencing reads to the reference genome |
| **samtools** | Converts, sorts, indexes, and checks BAM files |
| **Picard MarkDuplicates** | Marks PCR duplicate reads |
| **bcftools** | Calls variants from BAM files |
| **freebayes** | Calls variants using haplotype-based variant detection |
| **GATK HaplotypeCaller** | Calls variants using local haplotype assembly |
| **bedtools** | Compares variant calls between callers |
| **SnpSift** | Adds dbSNP annotation and filters variants |
| **snpEff** | Predicts functional effects of variants |
| **Ensembl VEP** | Provides web-based variant consequence annotation |

## Repository Structure

```text
kabuki-syndrome-wgs-variant-analysis/
│
├── README.md
├── kabuki_wgs_pipeline.sh
├── environment.yml
├── .gitignore
│
└── examples/
    └── vep_submission_example.vcf
```

Generated files are not stored in the repository. The pipeline creates them automatically when it runs.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/mahealamuq/kabuki-wgs-pipeline-GI4.git
cd kabuki-wgs-pipeline-GI4
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
```

### 3. Activate the environment

```bash
conda activate kabuki-wgs
```

### 4. Make the pipeline executable

```bash
chmod +x kabuki_wgs_pipeline.sh
```

## Running the Pipeline

Run the pipeline with default settings:

```bash
bash kabuki_wgs_pipeline.sh
```

Run with custom thread number:

```bash
bash kabuki_wgs_pipeline.sh --threads 8
```

Run with a custom output directory:

```bash
bash kabuki_wgs_pipeline.sh --threads 8 --workdir ~/kabuki-wgs-output
```

By default, the output directory is:

```text
~/kabuki-wgs-output
```

## Pipeline Steps Explained

### Step 1: Directory Setup

The script creates folders for input data, reference files, logs, and patient results.

Main folders:

```text
data/
ref/genome/
ref/dbsnp/
ref/intervals/
results/patientA/
results/patientB/
logs/
```

### Step 2: Reference Genome Preparation

The pipeline downloads chromosome 12 and chromosome X from the hg19 reference genome.

These chromosomes are selected because:

- **KMT2D** is located on chromosome 12
- **KDM6A** is located on chromosome X

The reference genome is indexed using BWA, samtools, and Picard so that alignment and variant calling tools can use it.

### Step 3: Input Read Preparation

The script downloads reads from public NA12878 WGS data and extracts only reads mapped to the KMT2D and KDM6A regions.

The extracted BAM data is converted into paired FASTQ files:

```text
patientA_1.fq.gz
patientA_2.fq.gz
patientB_1.fq.gz
patientB_2.fq.gz
```

### Step 4: Read Alignment

Reads are aligned to the hg19 reference genome using BWA MEM.

Output:

```text
patientA.sam
patientB.sam
```

The SAM file stores the alignment location of each sequencing read.

### Step 5: BAM Sorting and Indexing

SAM files are converted to BAM format, sorted by genomic coordinate, and indexed.

Output:

```text
patientA.sorted.bam
patientA.sorted.bam.bai
patientB.sorted.bam
patientB.sorted.bam.bai
```

Sorted and indexed BAM files are required for efficient downstream analysis.

### Step 6: Mark PCR Duplicates

Picard MarkDuplicates identifies duplicate reads that may have been produced during PCR amplification.

Output:

```text
patientA.markdups.bam
patientA.markdups_metrics.txt
patientB.markdups.bam
patientB.markdups_metrics.txt
```

Marking duplicates improves the reliability of variant calling.

### Step 7: Final BAM Preparation

The pipeline prepares final BAM files for variant calling.

GATK4 performs local realignment internally during variant calling, so the duplicate-marked BAM is used as the final analysis BAM.

Output:

```text
patientA.realigned.bam
patientB.realigned.bam
```

### Step 8: Alignment Quality Control

The pipeline uses samtools to assess alignment quality.

Commands used:

```bash
samtools idxstats
samtools flagstat
```

Output:

```text
patientA.qc_stats.txt
patientB.qc_stats.txt
```

These files report mapping statistics, chromosome-level read counts, properly paired reads, and duplicate read information.

### Step 9: Variant Calling with bcftools and freebayes

The pipeline calls variants using two variant callers:

```text
bcftools
freebayes
```

Output:

```text
patientA.bcftools.vcf
patientA.freebayes.vcf
patientB.bcftools.vcf
patientB.freebayes.vcf
```

Using more than one caller allows comparison of variant detection results.

### Step 10: Variant Calling with GATK HaplotypeCaller

The pipeline also calls variants using GATK HaplotypeCaller.

Output:

```text
patientA.gatk.vcf
patientB.gatk.vcf
```

This adds a third caller for comparison.

### Step 11: Caller Comparison

bedtools is used to compare variant calls between bcftools, freebayes, and GATK.

Output:

```text
patientA.caller_comparison.txt
patientB.caller_comparison.txt
```

This file includes variant counts, shared variants, and overlap statistics.

### Step 12: dbSNP Annotation

SnpSift annotates variants with dbSNP rsIDs.

Output:

```text
patientA.fb.dbsnp.vcf
patientB.fb.dbsnp.vcf
```

This helps distinguish known variants from novel or rare variants.

### Step 13: Functional Annotation with snpEff

snpEff predicts the biological consequence of each variant.

Output:

```text
patientA.fb.snpeff.vcf
patientB.fb.snpeff.vcf
```

snpEff can classify variants as:

```text
synonymous_variant
missense_variant
frameshift_variant
stop_gained
splice_region_variant
LOW impact
MODERATE impact
HIGH impact
```

### Step 14: Candidate Variant Filtering

The pipeline filters variants using the following criteria:

```text
HIGH or MODERATE impact
no dbSNP rsID
QUAL > 30
```

Output:

```text
patientA.candidates.vcf
patientB.candidates.vcf
```

These are candidate variants with potential functional importance.

### Step 15: VEP Submission File

The pipeline creates a VCF file for Ensembl VEP submission.

Output:

```text
results/vep_submission.vcf
```

Example:

```text
##fileformat=VCFv4.1
#CHROM  POS       ID  REF  ALT  QUAL     FILTER  INFO
12      49420214  .   G    A    322.788  .       .
X       44963994  .   C    T    322.000  .       .
```

Submit this file to Ensembl VEP GRCh37:

```text
https://grch37.ensembl.org/Homo_sapiens/Tools/VEP
```

## Main Output Files

| File | Description |
|---|---|
| `patientA.sorted.bam` | Sorted alignment file |
| `patientA.markdups.bam` | BAM file after duplicate marking |
| `patientA.realigned.bam` | Final BAM file used for variant calling |
| `patientA.qc_stats.txt` | Alignment QC summary |
| `patientA.bcftools.vcf` | Variants called by bcftools |
| `patientA.freebayes.vcf` | Variants called by freebayes |
| `patientA.gatk.vcf` | Variants called by GATK HaplotypeCaller |
| `patientA.caller_comparison.txt` | Comparison between variant callers |
| `patientA.fb.dbsnp.vcf` | VCF annotated with dbSNP IDs |
| `patientA.fb.snpeff.vcf` | VCF annotated with snpEff effects |
| `patientA.candidates.vcf` | Filtered candidate variants |
| `vep_submission.vcf` | VCF file for Ensembl VEP submission |

## Useful Result-Checking Commands

Count variants in a VCF file:

```bash
grep -v '^#' results/patientA/patientA.freebayes.vcf | wc -l
```

View candidate variants:

```bash
grep -v '^#' results/patientA/patientA.candidates.vcf | head
```

Count unique candidate positions:

```bash
grep -v '^#' results/patientA/patientA.candidates.vcf | cut -f1,2 | sort -u | wc -l
```

Check alignment statistics:

```bash
samtools flagstat results/patientA/patientA.realigned.bam
```

Check reads mapped to chr12 and chrX:

```bash
samtools idxstats results/patientA/patientA.realigned.bam
```

## Biological Interpretation

The final candidate variants are filtered to prioritize variants that may affect gene function.

A candidate disease-relevant variant is expected to have:

- a moderate or high predicted functional impact
- good variant quality
- no dbSNP rsID, suggesting it may be rare or novel

For Kabuki syndrome analysis, special attention is given to variants in:

```text
KMT2D
KDM6A
```

## Example Known Variants Used for VEP Practice

The pipeline writes known Kabuki-related example variants into the VEP submission file:

| Gene | Variant | Notes |
|---|---|---|
| KMT2D | `12:49420214 G>A` | Example missense variant, p.Arg5179Cys |
| KDM6A | `X:44963994 C>T` | Example variant for Kabuki syndrome type 2 region |

These variants can be submitted to VEP for consequence annotation and population frequency checking.

## Files to Upload to GitHub

Upload these files:

```text
README.md
kabuki_wgs_pipeline.sh
environment.yml
.gitignore
```

Do not upload generated sequencing output files because they are large and can be regenerated by the pipeline.

## Files Generated by the Pipeline

The following files and folders are created automatically when the pipeline runs:

```text
data/
ref/
results/
logs/
*.fastq.gz
*.sam
*.bam
*.bai
*.vcf
*.vcf.gz
```

## Conclusion

This project demonstrates a complete WGS variant analysis workflow for Kabuki syndrome-associated genes. It shows how public sequencing data can be processed through alignment, BAM processing, variant calling, variant comparison, annotation, filtering, and VEP-ready reporting.

The repository is suitable for demonstrating practical skills in human genomics, mutation detection, Linux-based bioinformatics, and reproducible pipeline development.
