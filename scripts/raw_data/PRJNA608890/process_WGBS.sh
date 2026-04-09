#!/usr/bin/env bash
set -euo pipefail

##############################
# User-defined directories
##############################
# Path to the reference genome FASTA file (for building the Bismark index)
GENOME_DIR="/lustre/groups/crna01/projects/collabs/shared/annots/GENECODE/hsapiens/release_44"
GENOME_FASTA="$GENOME_DIR/GRCh38.primary_assembly.genome.fa"
# GENOME_DIR is the directory where the Bismark genome index will be created.
# this folder should contain the reference FASTA file.

# WDI is the directory containing folders with the raw paired-end FASTQ files.
# Files are assumed to be named: replicate1_1.fastq.gz, replicate1_2.fastq.gz, etc.
WD="/lustre/groups/crna01/projects/collabs/gabi/ipsc/data_raw/PRJNA608890/WGBS"
TMP_DIR="$WD/tmp"
if [[ ! -d $TMP_DIR ]]; then
    mkdir -p "$TMP_DIR"
fi
##############################
# Step 1: Prepare the Genome Index
##############################
# If the genome index has not been prepared, copy the genome file and build the index.
## Filter the FASTA file to keep only autosomes and sex chromosomes to speed up bismark
### scaffolds and contigs consume a lot of memory and significatly slow down bismark downstream
FILTERED_GENOME_DIR="${GENOME_DIR}/autosexosomes"
FILTERED_FASTA="${FILTERED_GENOME_DIR}/GRCh38.primary_assembly.autosexosomes.fa"
if [[ ! -d "$FILTERED_GENOME_DIR" ]]; then
    mkdir -p "$FILTERED_GENOME_DIR"
elif [[ ! -f "$FILTERED_FASTA" ]]; then
    echo ""
    echo "Filtering genome to autosomes and sex chromosomes..."
    echo "-----------------------------------"
    awk '/^>/{f=($1 ~ /^>chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/)} f' "$GENOME_FASTA" > "$FILTERED_FASTA"
    samtools faidx "$FILTERED_FASTA"
elif [[ ! -d "${FILTERED_GENOME_DIR}/Bisulfite_Genome" ]]; then
    echo ""
    echo "Preparing Bismark genome index (autosomes + sex chromosomes only)..."
    echo "-----------------------------------"
    bismark_genome_preparation --parallel 4 "$FILTERED_GENOME_DIR"
else 
    echo "Filtered Bismark genome index already exists in ${FILTERED_GENOME_DIR}"
fi

##############################
# Define sample names
##############################
samples=("SRR11185407")

##############################
# Process each sample
##############################
for sample in "${samples[@]}"; do
    echo ""
    echo "Processing sample: ${sample}"
    echo "-----------------------------------"
    SD="$WD/$sample"
    TRIM_R1="${SD}/${sample}_1_val_1.fq.gz"
    TRIM_R2="${SD}/${sample}_2_val_2.fq.gz"
    if [[ ! -f "$TRIM_R1" ]]; then
        echo ""
        echo "Trimming reads...."
        echo "-----------------------------------"
        trim_galore --paired --fastqc \
            --clip_R1 10 --clip_R2 10 \
            --quality 20 --length 20 \
            -o "${SD}" \
            "${WD}/${sample}_1.fastq.gz" "${WD}/${sample}_2.fastq.gz" > "${SD}/${sample}_trimgalore.log"
    else
        echo ""
        echo "Found trimmed reads:"
        echo "$TRIM_R1"
        echo "$TRIM_R2"
        echo "-----------------------------------"
    fi
    
    
    BAM_FILE="${SD}/${sample}_1_val_1_bismark_bt2_pe.bam"
    echo ""
    if [[ -f "$BAM_FILE" ]]; then
        echo "Found aligned BAM file: $BAM_FILE"
        echo "-----------------------------------"
    else
        echo "Aligning reads with Bismark..."
        echo "-----------------------------------"
        bismark \
            -p 8 --parallel 2 \
            -o "${SD}" --temp_dir "$TMP_DIR" \
            --genome_folder "${FILTERED_GENOME_DIR}" \
            -1 "${TRIM_R1}" -2 "${TRIM_R2}"
    fi
    
    
    DEDUP_BAM="${BAM_FILE%.bam}.deduplicated.bam"
    echo ""
    if [[ -f "$DEDUP_BAM" ]]; then
        echo "Found deduplicated BAM file: $DEDUP_BAM"
        echo "-----------------------------------"
    else
        echo "Deduplicating BAM file..."
        echo "-----------------------------------"
        deduplicate_bismark --bam --paired --output_dir "${SD}" "${BAM_FILE}" 
    fi
    
    
    echo ""
    echo "Extracting methylation calls..."
    echo "-----------------------------------"
    bismark_methylation_extractor --parallel 4 --genome_folder "${FILTERED_GENOME_DIR}" --gzip --paired-end --comprehensive --bedGraph --report --output_dir "${SD}" "${DEDUP_BAM}" 
    
    echo "Finished processing sample: ${sample}"
    echo "-----------------------------------"
done
echo "All samples processed successfully."
