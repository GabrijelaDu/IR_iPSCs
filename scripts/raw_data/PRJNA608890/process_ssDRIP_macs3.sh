#!/usr/bin/env bash

# ----------- Call Strand-Specific R-loops ----------
WD="./ssDRIP"
NCORES=8
PD="$WD/peaks_macs"
TMP="$WD/tmp"

mkdir -p "$PD" "$TMP"

for BATCH in batch1 batch2; do
    echo -e "\nProcessing $BATCH..."

    BAMS=($(find "$WD" -type f -name "*_${BATCH}_rep*_Aligned.sortedByCoord.out_duplRm.bam"))
    if [[ ${#BAMS[@]} -eq 0 ]]; then
        echo "No BAM files found for batch $BATCH, skipping..."
        continue
    fi

    # Create arrays for forward and reverse strand BAMs
    FWD_BAMS=()
    REV_BAMS=()

    for BAM in "${BAMS[@]}"; do
        NAME=$(basename "$BAM" .bam)

        # Generate forward strand BAM (-F 16: NOT reverse)
        FWD="$TMP/${NAME}_fwd.bam"
        samtools view -@ "$NCORES" -b -F 16 "$BAM" -o "$FWD"
        samtools index "$FWD"
        FWD_BAMS+=("$FWD")

        # Generate reverse strand BAM (-f 16: only reverse)
        REV="$TMP/${NAME}_rev.bam"
        samtools view -@ "$NCORES" -b -f 16 "$BAM" -o "$REV"
        samtools index "$REV"
        REV_BAMS+=("$REV")
    done

    # Call peaks for forward strand
    echo "Calling forward strand peaks for $BATCH..."
    macs3 callpeak \
        -t "${FWD_BAMS[@]}" -f BAMPE \
        --keep-dup all --nomodel \
        -g hs -n "${BATCH}_fwd" \
        --tempdir "$TMP" --outdir "$PD"

    # Call peaks for reverse strand
    echo "Calling reverse strand peaks for $BATCH..."
    macs3 callpeak \
        -t "${REV_BAMS[@]}" -f BAMPE \
        --keep-dup all --nomodel \
        -g hs -n "${BATCH}_rev" \
        --tempdir "$TMP" --outdir "$PD"

    wait
done
