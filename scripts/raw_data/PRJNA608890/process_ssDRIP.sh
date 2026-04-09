#!/usr/bin/env bash

WD="./ssDRIP"
NCORES=8
GIDX="/lustre/groups/crna01/projects/collabs/shared/annots/GENECODE/hsapiens/release_44/STAR_index"
#
declare -A SAMPLE_LIST
while IFS=, read -r SRR SAMPLE; do
    SAMPLE_LIST["$SRR"]=$SAMPLE
done < "$WD/sample_list.csv"

echo ""
echo "processing samples:"
echo "${SAMPLE_LIST[@]}"
for SRR in "${!SAMPLE_LIST[@]}"; do
    #
    SAMPLE="${SAMPLE_LIST[$SRR]}"
    #
    echo ""
    echo "Processing $SAMPLE..."
    #
    SD="$WD/$SAMPLE"
    if [[ ! -d  "$SD" ]]; then
        mkdir -p "$SD"
    fi
    #
    if [[ -f "$WD/$SRR"_1.fastq.gz ]]; then
        mv "$WD/$SRR"_1.fastq.gz "$SD"/"$SAMPLE"_1.fastq.gz
    fi
    if [[ -f "$WD/$SRR"_2.fastq.gz ]]; then
        mv "$WD/$SRR"_2.fastq.gz "$SD"/"$SAMPLE"_2.fastq.gz
    fi
    # trim adapters & clip tails according to protocol reccomendations:
    # https://www.genetargetsolutions.com.au/wp-content/uploads/2015/05/Accel-NGS%E2%84%A2-1S-Plus-and-Methyl-Seq-Tail-Triming-for-Better-Data.pdf
    if [[ ! -f "$SD"/"$SAMPLE"_1_val_1.fastq.gz || ! -f "$SD"/"$SAMPLE"_2_val_2.fastq.gz ]]; then
        echo "Trimming $SAMPLE..."
        mkdir "$SD/fastqc"
        trim_galore \
            --paired --gzip --cores "$NCORES" \
            --fastqc_args "--outdir $SD/fastqc" \
            --stringency 3 \
            --clip_R1 10 --clip_R2 10 \
            --output_dir "$SD" \
            "$SD"/"$SAMPLE"_1.fastq.gz "$SD"/"$SAMPLE"_2.fastq.gz > "$SD"/"$SAMPLE"_trimgalore.Log
        wait
        # rename fq to fastq for consitency
        if compgen -G "$SD"/*.fq.gz > /dev/null; then
            for FASTA in "$SD"/*.fq.gz; do
                mv -- "$FASTA" "$SD"/"${FASTA%.fq.gz}.fastq.gz"
            done
        fi
    else
        echo "Already trimmed: $SAMPLE"
    fi
    # ------------- {MAP}
    if [[ ! -f "$SD"/"$SAMPLE"_Log.final.out ]]; then
        echo "Mapping $SAMPLE..."
        STAR \
            --runThreadN "$NCORES" \
            --genomeDir "$GIDX" \
            --readFilesCommand zcat \
            --readFilesIn "$SD"/"$SAMPLE"_1_val_1.fastq.gz "$SD"/"$SAMPLE"_2_val_2.fastq.gz \
            --outFileNamePrefix "$SD"/"$SAMPLE"_ \
            --outSAMtype BAM SortedByCoordinate
        wait
    else
        echo "Already mapped: $SAMPLE"
    fi
    if [[ ! -f "$SD"/"$SAMPLE"_Aligned.sortedByCoord.out_duplRm.bam ]]; then
        echo "Removing duplicates from $SAMPLE..."
        ##
        samtools index "$SD"/"$SAMPLE"_Aligned.sortedByCoord.out.bam
        ## {DEDUP}
        samtools sort -n -o "$SD"/tmp_nameSorted.bam "$SD"/"$SAMPLE"_Aligned.sortedByCoord.out.bam  # Sort by name
        samtools fixmate -m "$SD"/tmp_nameSorted.bam "$SD"/tmp_nameSorted_fixedMate.bam  # Add mate information
        samtools sort -o "$SD"/tmp_coordSorted_fixedMate.bam "$SD"/tmp_nameSorted_fixedMate.bam  # Sort by position
        samtools markdup -r "$SD"/tmp_coordSorted_fixedMate.bam "$SD"/"$SAMPLE"_Aligned.sortedByCoord.out_duplRm.bam  # Remove duplicates
        ##
        samtools index "$SD"/"$SAMPLE"_Aligned.sortedByCoord.out_duplRm.bam
        rm "$SD"/tmp_*
    else
        echo "Already deduplicated: $SAMPLE"
    fi
done
