#!/usr/bin/env bash

# Directory containing your raw paired-end FASTQ files.
# Files are assumed to be named: replicate1_R1.fastq.gz, replicate1_R2.fastq.gz, etc.
WD="/lustre/groups/crna01/projects/collabs/gabi/ipsc/data_raw/PRJNA608890/WGBS"
NCORES=4
samples=("SRR11185407")

for sample in "${samples[@]}"; do
    echo ""
    echo "Processing sample: ${sample}"
    SD="$WD/$sample"
    mkdir -p "$SD"
    # Step 2: Trim reads with Trim Galore (paired-end mode)
    # Trim Galore will output files with _val_1.fq.gz and _val_2.fq.gz suffixes.
    echo "Trimming reads...."
    trim_galore --paired --fastqc --cores "$NCORES" \
        --clip_R1 10 --clip_R2 10 \
        --quality 20 --length 20 \
        -o "${SD}" \
        "${WD}/${sample}_1.fastq.gz" "${WD}/${sample}_2.fastq.gz" > "${SD}/${sample}_trimgalore.log"
done
