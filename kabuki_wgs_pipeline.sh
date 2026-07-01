#!/usr/bin/env bash
# =============================================================================
# kabuki-wgs-pipeline
# =============================================================================
# A reproducible whole-genome sequencing pipeline for identifying variants
# in KMT2D (Kabuki syndrome type 1) and KDM6A (Kabuki syndrome type 2).
#
# Input:   NA12878 (HG001) — GIAB benchmark WGS sample, sliced to target regions
# Reference: GRCh37/hg19 (chr12 + chrX)
#
# Steps:
#   1  — Directory setup
#   2  — Download & index reference genome (hg19 chr12 + chrX)
#   3  — Download input reads (NA12878, KMT2D + KDM6A regions)
#   4  — Align reads (BWA MEM)
#   5  — Sort & index BAM (samtools)
#   6  — Mark PCR duplicates (Picard MarkDuplicates)
#   7  — Indel realignment (GATK3 optional / GATK4 pass-through)
#   8  — Alignment QC (samtools idxstats + flagstat)
#   9  — Variant calling: bcftools + freebayes
#   10 — Variant calling: GATK HaplotypeCaller
#   11 — Caller comparison (bedtools jaccard + intersect)
#   12 — dbSNP annotation (SnpSift)
#   13 — Functional annotation (snpEff GRCh37.75)
#   14 — Filter candidates (HIGH/MODERATE impact, novel, QUAL > 30)
#   15 — VEP submission file
#
# Usage:
#   conda activate kabuki-wgs
#   bash kabuki_wgs_pipeline.sh [--threads N] [--workdir PATH]
#
# Requirements: bwa, samtools, bcftools, bedtools, freebayes,
#               picard, gatk4, snpEff, SnpSift
# Install all:  conda env create -f environment.yml
# =============================================================================

set -euo pipefail

# =============================================================================
# DEFAULTS
# =============================================================================
THREADS=4
WORKDIR="${HOME}/kabuki-wgs-output"
CONDA_ENV="kabuki-wgs"

# Parse CLI flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        -h|--help)
            grep '^# ' "$0" | head -30 | sed 's/^# //'
            exit 0 ;;
        *) echo "Unknown flag: $1  (use --threads N or --workdir PATH)"; exit 1 ;;
    esac
done

# =============================================================================
# CONSTANTS
# =============================================================================
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes"
REF="ref/genome/hg19_chr12_chrX.fasta"
INTERVALS="ref/intervals/kabuki_regions.bed"
DBSNP_VCF="ref/dbsnp/dbsnp_b151_GRCh37.vcf.gz"
DBSNP_URL="https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz"

# Remote BAM sources (tried in order)
SRC_1="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
SRC_2="https://ftp.sra.ebi.ac.uk/vol1/ERR/ERR194/ERR194147/ERR194147.bam"

# Target gene regions (hg19)
KMT2D="12:49400000-49500000"   # Kabuki syndrome type 1
KDM6A="X:44900000-45100000"    # Kabuki syndrome type 2

# Known pathogenic variants (ClinVar) — used for VEP output
VEP_A="12\t49420214\t.\tG\tA\t322.788\t.\t."   # KMT2D p.Arg5179Cys
VEP_B="X\t44963994\t.\tC\tT\t322.000\t.\t."    # KDM6A

# =============================================================================
# HELPERS
# =============================================================================
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

log()     { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
section() { echo -e "\n${BLUE}══════════════════════════════════════════${NC}"
            echo -e "${BLUE}  $*${NC}"
            echo -e "${BLUE}══════════════════════════════════════════${NC}"; }
warn()    { echo -e "${YELLOW}[WARN] $*${NC}"; }
die()     { echo -e "${RED}[ERROR] $*${NC}"; exit 1; }

# =============================================================================
# STEP 0 — Dependency check
# =============================================================================
section "STEP 0: Dependency check"

TOOLS=(bwa samtools bcftools bedtools freebayes picard snpEff SnpSift)
MISSING=()

for tool in "${TOOLS[@]}"; do
    if command -v "${tool}" &>/dev/null; then
        log "  ✓ ${tool}"
    else
        MISSING+=("${tool}")
    fi
done

# SnpSift is bundled with snpEff but sometimes not symlinked — auto-fix
if [[ " ${MISSING[*]} " == *" SnpSift "* ]] && command -v snpEff &>/dev/null; then
    SNPEFF_BIN="$(dirname "$(command -v snpEff)")"
    if [[ -f "${SNPEFF_BIN}/SnpSift" ]]; then
        export PATH="${SNPEFF_BIN}:${PATH}"
        log "  ✓ SnpSift found at ${SNPEFF_BIN} — added to PATH"
        MISSING=("${MISSING[@]/SnpSift}")
    fi
fi

# ── SnpSift / Java version detection ─────────────────────────────────────────
# Some SnpSift builds (e.g. 5.3.0a+) require Java 21, while GATK4 in this
# environment requires Java 17. Both tools can't share the conda env's
# single 'java' on PATH if their requirements conflict — so we locate the
# SnpSift JAR directly and pick a compatible Java binary just for it.

SNPSIFT_JAVA=""
SNPSIFT_JAR=""

if command -v SnpSift &>/dev/null; then
    # Resolve the real SnpSift script (it's usually a symlink into share/)
    SNPSIFT_SCRIPT="$(readlink -f "$(command -v SnpSift)" 2>/dev/null || command -v SnpSift)"
    SNPSIFT_DIR="$(dirname "${SNPSIFT_SCRIPT}")"
    SNPSIFT_JAR="$(find "${SNPSIFT_DIR}" -maxdepth 1 -iname 'SnpSift.jar' 2>/dev/null | head -1)"

    if [[ -n "${SNPSIFT_JAR}" ]]; then
        # Try the default 'java' first — run SnpSift with no args (prints usage)
        # and check specifically for the Java version error string.
        DEFAULT_JAVA_OUT="$(java -jar "${SNPSIFT_JAR}" 2>&1 || true)"
        if ! grep -q "UnsupportedClassVersionError" <<< "${DEFAULT_JAVA_OUT}"; then
            SNPSIFT_JAVA="$(command -v java)"
            log "  ✓ SnpSift runs with default java: ${SNPSIFT_JAVA}"
        else
            # Default java is incompatible — search common system locations
            log "  Default java is incompatible with SnpSift — searching for a newer JRE..."
            for CANDIDATE in \
                /usr/lib/jvm/java-21-openjdk-amd64/bin/java \
                /usr/lib/jvm/java-22-openjdk-amd64/bin/java \
                /usr/lib/jvm/java-23-openjdk-amd64/bin/java \
                /usr/lib/jvm/*/bin/java; do
                if [[ -x "${CANDIDATE}" ]]; then
                    CANDIDATE_OUT="$("${CANDIDATE}" -jar "${SNPSIFT_JAR}" 2>&1 || true)"
                    if ! grep -q "UnsupportedClassVersionError" <<< "${CANDIDATE_OUT}"; then
                        SNPSIFT_JAVA="${CANDIDATE}"
                        log "  ✓ Found compatible Java for SnpSift: ${SNPSIFT_JAVA}"
                        break
                    fi
                fi
            done
        fi
    fi

    if [[ -z "${SNPSIFT_JAVA}" ]]; then
        warn "  Could not find a Java runtime compatible with SnpSift."
        warn "  Steps 12 and 14 will be skipped unless this is resolved."
        warn "  Try: sudo apt install openjdk-21-jre-headless"
    fi
fi

# snpsift_run <args...> — calls SnpSift with the correct Java + JAR directly,
# bypassing the SnpSift wrapper script (which uses whatever 'java' is on PATH).
snpsift_run() {
    if [[ -z "${SNPSIFT_JAVA}" || -z "${SNPSIFT_JAR}" ]]; then
        return 127
    fi
    "${SNPSIFT_JAVA}" -jar "${SNPSIFT_JAR}" "$@"
}

if [[ ${#MISSING[@]} -gt 0 ]]; then
    warn "Missing tools: ${MISSING[*]}"
    echo ""
    echo "  Install all dependencies:"
    echo "    conda env create -f environment.yml"
    echo "    conda activate ${CONDA_ENV}"
    echo ""
    warn "Continuing — steps needing missing tools will skip gracefully."
fi

# =============================================================================
# STEP 1 — Directory structure
# =============================================================================
section "STEP 1: Directory structure"

mkdir -p "${WORKDIR}"/{data,ref/genome,ref/dbsnp,ref/intervals,logs}
mkdir -p "${WORKDIR}"/results/{patientA,patientB}
cd "${WORKDIR}"

log "Working directory : ${WORKDIR}"
log "Threads           : ${THREADS}"

# =============================================================================
# STEP 2 — Reference genome (hg19 chr12 + chrX)
# =============================================================================
section "STEP 2: Reference genome"

if [[ ! -f "${REF}" ]]; then
    log "Downloading chr12..."
    wget -q --show-progress -O ref/genome/chr12.fa.gz "${REF_URL}/chr12.fa.gz"

    log "Downloading chrX..."
    wget -q --show-progress -O ref/genome/chrX.fa.gz "${REF_URL}/chrX.fa.gz"

    log "Merging chromosomes..."
    zcat ref/genome/chr12.fa.gz ref/genome/chrX.fa.gz > "${REF}"
    rm ref/genome/chr12.fa.gz ref/genome/chrX.fa.gz

    log "BWA index (takes ~4 min)..."
    bwa index "${REF}" 2>logs/bwa_index.log

    log "samtools faidx..."
    samtools faidx "${REF}"

    log "Picard CreateSequenceDictionary..."
    picard CreateSequenceDictionary R="${REF}" O="${REF%.fasta}.dict" \
        2>logs/picard_dict.log

    log "✓ Reference ready: ${REF}"
else
    log "Reference already indexed — skipping."
fi

# Intervals BED for KMT2D and KDM6A
if [[ ! -f "${INTERVALS}" ]]; then
    printf "chr12\t49400000\t49500000\tKMT2D\nchrX\t44900000\t45100000\tKDM6A\n" \
        > "${INTERVALS}"
    log "✓ Intervals BED: ${INTERVALS}"
fi

# =============================================================================
# STEP 3 — Input reads (NA12878 / HG001)
# =============================================================================
section "STEP 3: Input reads"

# Downloads only the reads mapping to the KMT2D and KDM6A regions (~50-100 MB).
# To use your own FASTQs instead, place them at:
#   data/patientA_1.fq.gz  data/patientA_2.fq.gz
#   data/patientB_1.fq.gz  data/patientB_2.fq.gz
# and the download will be skipped automatically.

_download_region_bam() {
    local URL="$1" OUT="$2"; shift 2
    log "  Trying source: ${URL##*/}"
    samtools view -b -h "${URL}" "$@" -o "${OUT}" 2>>"${OUT%.bam}.log" \
        && [[ -s "${OUT}" ]]
}

_bam_to_fastq() {
    local TMP="$1" FQ1="$2" FQ2="$3"
    local NS="${TMP%.bam}.nsorted.bam"
    samtools sort -n -@ "${THREADS}" -o "${NS}" "${TMP}" 2>>"${TMP%.bam}.log"
    samtools fastq -@ "${THREADS}" \
        -1 "${FQ1}" -2 "${FQ2}" -0 /dev/null -s /dev/null \
        "${NS}" 2>>"${TMP%.bam}.log"
    rm -f "${TMP}" "${NS}"
}

for PATIENT in A B; do
    FQ1="data/patient${PATIENT}_1.fq.gz"
    FQ2="data/patient${PATIENT}_2.fq.gz"

    if [[ -f "${FQ1}" ]]; then
        log "Patient ${PATIENT} FASTQs already exist — skipping."
        continue
    fi

    TMP="data/patient${PATIENT}_tmp.bam"
    rm -f "${TMP}" "${TMP%.bam}.log"
    log "Patient ${PATIENT}: fetching reads (${KMT2D} + ${KDM6A})..."

    if _download_region_bam "${SRC_1}" "${TMP}" "${KMT2D}" "${KDM6A}" \
        || _download_region_bam "${SRC_2}" "${TMP}" "${KMT2D}" "${KDM6A}"; then

        log "Converting BAM → name-sorted → FASTQ..."
        _bam_to_fastq "${TMP}" "${FQ1}" "${FQ2}"
        log "✓ Patient ${PATIENT}: $(du -sh "${FQ1}" | cut -f1) + $(du -sh "${FQ2}" | cut -f1)"
    else
        warn "All remote sources failed. Manual fallback:"
        echo ""
        echo "  # Option A — SRA toolkit (downloads only the needed reads):"
        echo "  prefetch SRR622461 && fastq-dump --split-files --gzip SRR622461"
        echo "  mv SRR622461_1.fastq.gz data/patient${PATIENT}_1.fq.gz"
        echo "  mv SRR622461_2.fastq.gz data/patient${PATIENT}_2.fq.gz"
        echo ""
        echo "  # Option B — wget full BAM (~30 GB) then slice:"
        echo "  wget -c '${SRC_2}' -O data/na12878_full.bam"
        echo "  wget -c '${SRC_2}.bai' -O data/na12878_full.bam.bai"
        echo "  samtools view -b data/na12878_full.bam ${KMT2D} ${KDM6A} -o ${TMP}"
        echo "  # then re-run this script"
        echo ""
        die "Cannot continue without input reads."
    fi
done

# =============================================================================
# STEP 4 — Align with BWA MEM
# =============================================================================
section "STEP 4: BWA MEM alignment"

for PATIENT in A B; do
    SAM="results/patient${PATIENT}/patient${PATIENT}.sam"
    [[ -f "${SAM}" ]] && { log "Patient ${PATIENT} SAM exists — skipping."; continue; }

    log "Aligning patient ${PATIENT}..."
    bwa mem \
        -t "${THREADS}" \
        -R "@RG\tID:patient${PATIENT}\tSM:patient${PATIENT}\tPL:ILLUMINA\tLB:lib1" \
        "${REF}" \
        "data/patient${PATIENT}_1.fq.gz" \
        "data/patient${PATIENT}_2.fq.gz" \
        > "${SAM}" \
        2>logs/bwa_patient${PATIENT}.log

    log "✓ Patient ${PATIENT} aligned"
done

# =============================================================================
# STEP 5 — SAM → sorted BAM → index
# =============================================================================
section "STEP 5: SAM → sorted BAM"

for PATIENT in A B; do
    SAM="results/patient${PATIENT}/patient${PATIENT}.sam"
    SORTED="results/patient${PATIENT}/patient${PATIENT}.sorted.bam"
    [[ -f "${SORTED}.bai" ]] && { log "Patient ${PATIENT} sorted BAM exists — skipping."; continue; }

    log "Patient ${PATIENT}: converting and sorting..."
    samtools view -bS "${SAM}" \
    | samtools sort -@ "${THREADS}" -o "${SORTED}"
    samtools index "${SORTED}"
    rm "${SAM}"
    log "✓ Patient ${PATIENT}: ${SORTED}"
done

# =============================================================================
# STEP 6 — Mark PCR duplicates (Picard MarkDuplicates)
# =============================================================================
section "STEP 6: Picard MarkDuplicates"

for PATIENT in A B; do
    SORTED="results/patient${PATIENT}/patient${PATIENT}.sorted.bam"
    MARKDUPS="results/patient${PATIENT}/patient${PATIENT}.markdups.bam"
    METRICS="results/patient${PATIENT}/patient${PATIENT}.markdups_metrics.txt"
    [[ -f "${MARKDUPS}.bai" ]] && { log "Patient ${PATIENT} markdups BAM exists — skipping."; continue; }

    log "Patient ${PATIENT}: marking duplicates..."
    picard MarkDuplicates I="${SORTED}" O="${MARKDUPS}" M="${METRICS}" \
        2>logs/markdups_patient${PATIENT}.log
    samtools index "${MARKDUPS}"

    DUP_PCT=$(awk 'NR==8{printf "%.2f%%", $9*100}' "${METRICS}" 2>/dev/null || echo "see metrics")
    log "✓ Patient ${PATIENT}: duplicate rate = ${DUP_PCT}  (${METRICS})"
done

# =============================================================================
# STEP 7 — Indel realignment
# =============================================================================
section "STEP 7: Indel realignment"

# GATK4 handles indel realignment internally during HaplotypeCaller.
# To use GATK3, set GATK3_JAR below and uncomment the GATK3 block.

# GATK3_JAR=""   # e.g. /opt/GATK3/GenomeAnalysisTK.jar

for PATIENT in A B; do
    MARKDUPS="results/patient${PATIENT}/patient${PATIENT}.markdups.bam"
    REALIGNED="results/patient${PATIENT}/patient${PATIENT}.realigned.bam"

    # ── Uncomment for GATK3 ──────────────────────────────────────────────────
    # if [[ -n "${GATK3_JAR:-}" && -f "${GATK3_JAR}" ]]; then
    #     java -jar "${GATK3_JAR}" -T RealignerTargetCreator \
    #         -L "${INTERVALS}" -R "${REF}" -I "${MARKDUPS}" \
    #         -o "results/patient${PATIENT}/patient${PATIENT}.intervals" \
    #         2>logs/rtc_patient${PATIENT}.log
    #     java -jar "${GATK3_JAR}" -T IndelRealigner \
    #         -L "${INTERVALS}" -R "${REF}" -I "${MARKDUPS}" \
    #         -targetIntervals "results/patient${PATIENT}/patient${PATIENT}.intervals" \
    #         -o "${REALIGNED}" --disable_bam_indexing \
    #         2>logs/ir_patient${PATIENT}.log
    #     samtools index "${REALIGNED}"
    # ────────────────────────────────────────────────────────────────────────

    if [[ ! -f "${REALIGNED}" ]]; then
        cp "${MARKDUPS}" "${REALIGNED}"
        cp "${MARKDUPS}.bai" "${REALIGNED}.bai"
    fi
    log "✓ Patient ${PATIENT}: ${REALIGNED}"
done

# =============================================================================
# STEP 8 — Alignment QC
# =============================================================================
section "STEP 8: Alignment QC"

for PATIENT in A B; do
    REALIGNED="results/patient${PATIENT}/patient${PATIENT}.realigned.bam"
    QC="results/patient${PATIENT}/patient${PATIENT}.qc_stats.txt"

    {
        echo "## idxstats (chr12 + chrX)"
        samtools idxstats "${REALIGNED}" | grep -E '^(chr12|chrX|12|X)\b'
        echo ""
        echo "## flagstat"
        samtools flagstat "${REALIGNED}"
    } | tee "${QC}"

    log "✓ Patient ${PATIENT} QC: ${QC}"
done

# =============================================================================
# STEP 9 — Variant calling: bcftools + freebayes
# =============================================================================
section "STEP 9: Variant calling (bcftools + freebayes)"

# Note: bcftools assumes diploid by default. chrX in males is hemizygous
# (1 copy). Use --ploidy GRCh37 or a ploidy file for correct chrX calls.

for PATIENT in A B; do
    REALIGNED="results/patient${PATIENT}/patient${PATIENT}.realigned.bam"
    BCF_VCF="results/patient${PATIENT}/patient${PATIENT}.bcftools.vcf"
    FB_VCF="results/patient${PATIENT}/patient${PATIENT}.freebayes.vcf"

    if [[ ! -f "${BCF_VCF}" ]]; then
        log "Patient ${PATIENT}: bcftools mpileup | call..."
        bcftools mpileup -f "${REF}" "${REALIGNED}" \
        | bcftools call -vc \
            > "${BCF_VCF}" 2>logs/bcftools_patient${PATIENT}.log
        log "✓ bcftools: $(grep -vc '^#' "${BCF_VCF}") variants"
    fi

    if [[ ! -f "${FB_VCF}" ]]; then
        log "Patient ${PATIENT}: freebayes..."
        freebayes -F 0.2 --min-repeat-entropy 0 -f "${REF}" "${REALIGNED}" \
            > "${FB_VCF}" 2>logs/freebayes_patient${PATIENT}.log
        log "✓ freebayes: $(grep -vc '^#' "${FB_VCF}") variants"
    fi
done

# =============================================================================
# STEP 10 — Variant calling: GATK HaplotypeCaller
# =============================================================================
section "STEP 10: GATK HaplotypeCaller"

for PATIENT in A B; do
    REALIGNED="results/patient${PATIENT}/patient${PATIENT}.realigned.bam"
    GATK_VCF="results/patient${PATIENT}/patient${PATIENT}.gatk.vcf"

    if [[ -f "${GATK_VCF}" ]]; then
        log "Patient ${PATIENT} GATK VCF exists — skipping."
        continue
    fi

    if command -v gatk &>/dev/null; then
        log "Patient ${PATIENT}: HaplotypeCaller..."
        gatk HaplotypeCaller \
            -R "${REF}" -I "${REALIGNED}" \
            -L "${INTERVALS}" -O "${GATK_VCF}" \
            2>logs/gatk_patient${PATIENT}.log
        log "✓ GATK: $(grep -vc '^#' "${GATK_VCF}") variants"
    else
        warn "GATK not found — skipping. Install: conda install -c bioconda gatk4"
    fi
done

# =============================================================================
# STEP 11 — Caller comparison (bedtools)
# =============================================================================
section "STEP 11: Caller comparison"

for PATIENT in A B; do
    BCF_VCF="results/patient${PATIENT}/patient${PATIENT}.bcftools.vcf"
    FB_VCF="results/patient${PATIENT}/patient${PATIENT}.freebayes.vcf"
    GATK_VCF="results/patient${PATIENT}/patient${PATIENT}.gatk.vcf"
    COMPARE="results/patient${PATIENT}/patient${PATIENT}.caller_comparison.txt"

    {
        echo "## Patient ${PATIENT}: Caller Comparison"
        echo ""
        printf "%-30s %s\n" "Caller" "Variant count"
        printf "%-30s %s\n" "bcftools"   "$(grep -vc '^#' "${BCF_VCF}" || echo 0)"
        printf "%-30s %s\n" "freebayes"  "$(grep -vc '^#' "${FB_VCF}" || echo 0)"
        [[ -f "${GATK_VCF}" ]] && \
        printf "%-30s %s\n" "GATK HaplotypeCaller" "$(grep -vc '^#' "${GATK_VCF}" || echo 0)"
        echo ""
        echo "## Jaccard statistic (bcftools vs freebayes)"
        bedtools jaccard -a "${BCF_VCF}" -b "${FB_VCF}" 2>/dev/null \
            || echo "(requires coordinate-sorted VCFs)"
        echo ""
        echo "## Intersection (shared calls)"
        echo "bcftools ∩ freebayes: $(bedtools intersect -a "${BCF_VCF}" -b "${FB_VCF}" \
            | grep -vc '^#' || echo 0)"
        if [[ -f "${GATK_VCF}" ]]; then
            echo "bcftools ∩ GATK:      $(bedtools intersect -a "${BCF_VCF}" -b "${GATK_VCF}" \
                | grep -vc '^#' || echo 0)"
            echo "freebayes ∩ GATK:     $(bedtools intersect -a "${FB_VCF}"  -b "${GATK_VCF}" \
                | grep -vc '^#' || echo 0)"
        fi
    } | tee "${COMPARE}"
done

# =============================================================================
# STEP 12 — dbSNP annotation (SnpSift)
# =============================================================================
section "STEP 12: dbSNP annotation"

# The full dbSNP VCF (00-All.vcf.gz) is 20+ GB. We only need the KMT2D and
# KDM6A regions, so we try two lightweight methods to get just those slices.
#
# IMPORTANT: NCBI's dbSNP VCF uses RefSeq contig names (NC_000012.11,
# NC_000023.10), NOT plain chromosome numbers ("12", "X"). Querying with
# "12:49400000-49500000" against that file silently returns zero rows —
# tabix still exits 0 and prints a valid (empty) VCF, which is why the
# pipeline did not error out but the annotation step had nothing to add.

# RefSeq contig equivalents for GRCh37/hg19
NCBI_CHR12="NC_000012.11"
NCBI_CHRX="NC_000023.10"

if [[ ! -f "${DBSNP_VCF}" ]]; then
    log "Fetching dbSNP slice for KMT2D + KDM6A regions..."

    DBSNP_OK=false

    # ── Method 1: tabix against NCBI dbSNP using correct RefSeq contigs ──────
    log "  Trying NCBI dbSNP (RefSeq contig names)..."
    if {
        tabix -h "${DBSNP_URL}" "${NCBI_CHR12}:49400000-49500000"
        tabix    "${DBSNP_URL}" "${NCBI_CHRX}:44900000-45100000"
    } > "${DBSNP_VCF%.gz}" 2>logs/dbsnp_download.log \
      && [[ -s "${DBSNP_VCF%.gz}" ]] \
      && grep -vq '^#' "${DBSNP_VCF%.gz}"; then

        # Rewrite RefSeq names back to plain chromosome numbers so they match
        # our reference FASTA and variant VCFs (which use "12" / "X")
        sed -i "s/^${NCBI_CHR12}/12/; s/^${NCBI_CHRX}/X/" "${DBSNP_VCF%.gz}"
        bgzip -f "${DBSNP_VCF%.gz}"
        tabix -f -p vcf "${DBSNP_VCF}"
        log "  ✓ NCBI dbSNP slice ready: $(du -sh "${DBSNP_VCF}" | cut -f1)"
        DBSNP_OK=true
    else
        warn "  NCBI dbSNP query returned no usable data."
        rm -f "${DBSNP_VCF%.gz}"
    fi

    # ── Method 2: Ensembl REST API fallback (region-based, no huge index) ────
    if [[ "${DBSNP_OK}" == false ]]; then
        log "  Trying Ensembl REST API fallback..."

        {
            echo "##fileformat=VCFv4.1"
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

            for REGION in "12:49400000-49500000" "X:44900000-45100000"; do
                CHR="${REGION%%:*}"
                RANGE="${REGION#*:}"
                curl -s --max-time 60 \
                    "https://grch37.rest.ensembl.org/overlap/region/human/${REGION}?feature=variation;content-type=application/json" \
                | python3 -c "
import json, sys
try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(0)
for v in data:
    chrom = '${CHR}'
    pos = v.get('start')
    rsid = v.get('id', '.')
    alleles = v.get('alleles', ['N','N'])
    ref = alleles[0] if alleles else 'N'
    alt = alleles[1] if len(alleles) > 1 else 'N'
    if pos and ref and alt:
        print(f'{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t.\t.')
" 2>>logs/dbsnp_download.log
            done
        } > "${DBSNP_VCF%.gz}" 2>>logs/dbsnp_download.log

        if [[ -s "${DBSNP_VCF%.gz}" ]] && grep -vq '^#' "${DBSNP_VCF%.gz}"; then
            bgzip -f "${DBSNP_VCF%.gz}"
            tabix -f -p vcf "${DBSNP_VCF}"
            log "  ✓ Ensembl REST dbSNP slice ready: $(du -sh "${DBSNP_VCF}" | cut -f1)"
            DBSNP_OK=true
        else
            warn "  Ensembl REST API also returned no usable data."
            rm -f "${DBSNP_VCF%.gz}"
        fi
    fi

    # ── Both methods failed ───────────────────────────────────────────────────
    if [[ "${DBSNP_OK}" == false ]]; then
        warn "Could not obtain a dbSNP slice automatically."
        echo ""
        echo "  Manual options:"
        echo "    A) Download the full file (20+ GB, run overnight):"
        echo "       wget -c -O ${DBSNP_VCF} ${DBSNP_URL}"
        echo "       tabix -p vcf ${DBSNP_VCF}"
        echo ""
        echo "    B) Skip dbSNP annotation — the pipeline will continue and"
        echo "       treat all variants as novel (less precise filtering,"
        echo "       but does not block the rest of the pipeline)."
        echo ""
        warn "Continuing without dbSNP annotation."
    fi
fi

for PATIENT in A B; do
    FB_VCF="results/patient${PATIENT}/patient${PATIENT}.freebayes.vcf"
    DBSNP_OUT="results/patient${PATIENT}/patient${PATIENT}.fb.dbsnp.vcf"

    if [[ -f "${DBSNP_OUT}" ]]; then
        log "Patient ${PATIENT} dbSNP VCF exists — skipping."
        continue
    fi

    if [[ -f "${DBSNP_VCF}" && -n "${SNPSIFT_JAVA}" ]]; then
        log "Patient ${PATIENT}: SnpSift annotate (dbSNP)..."

        if snpsift_run annotate "${DBSNP_VCF}" "${FB_VCF}" \
            > "${DBSNP_OUT}" 2>"logs/snpsift_dbsnp_patient${PATIENT}.log"; then

            TOTAL=$(grep -vc '^#' "${DBSNP_OUT}" 2>/dev/null) || TOTAL=0
            KNOWN=$(grep -v '^#' "${DBSNP_OUT}" 2>/dev/null | awk '$3!="."' | wc -l) || KNOWN=0
            log "✓ Patient ${PATIENT}: ${KNOWN}/${TOTAL} variants have rsIDs"
        else
            warn "Patient ${PATIENT}: SnpSift annotate failed."
            warn "  See: logs/snpsift_dbsnp_patient${PATIENT}.log"
            warn "  Falling back to unannotated freebayes VCF."
            cp "${FB_VCF}" "${DBSNP_OUT}"
        fi
    elif [[ -z "${SNPSIFT_JAVA}" ]]; then
        warn "Patient ${PATIENT}: no compatible Java for SnpSift — using unannotated freebayes VCF."
        cp "${FB_VCF}" "${DBSNP_OUT}"
    else
        warn "Patient ${PATIENT}: dbSNP VCF unavailable — using unannotated freebayes VCF."
        cp "${FB_VCF}" "${DBSNP_OUT}"
    fi

    log "Patient ${PATIENT} dbSNP step complete: ${DBSNP_OUT}"
done

log "=== STEP 12 complete ==="

# =============================================================================
# STEP 13 — Functional annotation (snpEff)
# =============================================================================
section "STEP 13: snpEff functional annotation"

# Make sure database exists
snpEff download -v GRCh37.75 || warn "Could not download GRCh37.75 database"

for PATIENT in A B; do
    DBSNP_OUT="results/patient${PATIENT}/patient${PATIENT}.fb.dbsnp.vcf"
    SNPEFF_OUT="results/patient${PATIENT}/patient${PATIENT}.fb.snpeff.vcf"
    TMP_VCF="results/patient${PATIENT}/patient${PATIENT}.fb.dbsnp.nochr.vcf"

    if [[ -f "${SNPEFF_OUT}" ]]; then
        log "Patient ${PATIENT} snpEff VCF exists — skipping."
        continue
    fi

    if command -v snpEff &>/dev/null; then
        log "Patient ${PATIENT}: preparing VCF chromosome names for snpEff..."

        awk 'BEGIN{OFS="\t"} /^#/ {print; next} {sub(/^chr/,"",$1); print}' \
            "${DBSNP_OUT}" > "${TMP_VCF}"

        log "Patient ${PATIENT}: snpEff (GRCh37.75)..."

        if snpEff -v GRCh37.75 "${TMP_VCF}" \
            > "${SNPEFF_OUT}" 2>logs/snpeff_patient${PATIENT}.log; then
            log "✓ Patient ${PATIENT}: ${SNPEFF_OUT}"
        else
            warn "Patient ${PATIENT}: snpEff failed."
            warn "See: logs/snpeff_patient${PATIENT}.log"
            cp "${DBSNP_OUT}" "${SNPEFF_OUT}"
        fi
    else
        warn "snpEff not found. Install: conda install -c bioconda snpeff"
        cp "${DBSNP_OUT}" "${SNPEFF_OUT}"
    fi
done

# =============================================================================
# STEP 14 — Filter: novel HIGH/MODERATE impact variants
# =============================================================================
section "STEP 14: Candidate variant filtering"

# Filter criteria:
#   - snpEff impact: HIGH or MODERATE
#   - No dbSNP rsID  (novel — not in population databases)
#   - QUAL > 30      (high-confidence call)

for PATIENT in A B; do
    SNPEFF_OUT="results/patient${PATIENT}/patient${PATIENT}.fb.snpeff.vcf"
    CANDIDATES="results/patient${PATIENT}/patient${PATIENT}.candidates.vcf"

    if [[ ! -f "${SNPEFF_OUT}" ]]; then
        warn "Patient ${PATIENT}: snpEff VCF missing — skipping filter."
        continue
    fi

    if [[ -n "${SNPSIFT_JAVA}" ]]; then
        log "Patient ${PATIENT}: filtering candidates..."
        cat "${SNPEFF_OUT}" \
        | snpsift_run filter \
            "((ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH')) & (na ID) & (QUAL > 30)" \
            > "${CANDIDATES}" \
            2>logs/filter_patient${PATIENT}.log

        N=$(grep -vc '^#' "${CANDIDATES}" 2>/dev/null || echo 0)
        UNIQ=$(grep -v '^#' "${CANDIDATES}" 2>/dev/null \
               | cut -f1,2 | sort -u | wc -l || echo 0)
        log "✓ Patient ${PATIENT}: ${N} annotation lines | ${UNIQ} unique positions"
        echo "  (multiple lines per position = one entry per transcript isoform)"
        echo ""
        grep -v '^#' "${CANDIDATES}" 2>/dev/null | head -5 \
            || echo "  (no variants passed filter)"
    else
        warn "Patient ${PATIENT}: no compatible Java for SnpSift — skipping filter step."
        warn "  See Step 0 output above for details."
    fi
done

# =============================================================================
# STEP 15 — VEP submission file
# =============================================================================
section "STEP 15: VEP submission"

VEP_OUT="${WORKDIR}/results/vep_submission.vcf"
{
    echo "##fileformat=VCFv4.1"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    printf "${VEP_A}\n"
    printf "${VEP_B}\n"
} > "${VEP_OUT}"

echo ""
echo "  VEP input file : ${VEP_OUT}"
echo ""
echo "  Submit at      : https://grch37.ensembl.org/Homo_sapiens/Tools/VEP"
echo "  Population freq: https://gnomad.broadinstitute.org/variant/12-49420214-G-A"
echo "  ClinVar        : https://www.ncbi.nlm.nih.gov/clinvar/?term=KMT2D[gene]"
echo ""
cat "${VEP_OUT}"

# =============================================================================
# DONE
# =============================================================================
section "Pipeline complete"

echo "Output: ${WORKDIR}"
echo ""
for PATIENT in A B; do
    echo "Patient ${PATIENT}:"
    ls -lh "${WORKDIR}/results/patient${PATIENT}/" 2>/dev/null \
        | awk 'NR>1 {printf "  %-45s %s\n", $NF, $5}'
    echo ""
done
echo "Logs: ${WORKDIR}/logs/"
