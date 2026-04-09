#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 accessions.txt"
    exit 1
fi

ACCESSION_FILE=$1
NCORES=4

while IFS= read -r ACCESSION || [[ -n "$ACCESSION" ]]; do
    if [ -z "$ACCESSION" ]; then
        continue
    fi
    echo "Downloading $ACCESSION..."
    prefetch -X 60G  "$ACCESSION"
    if [ $? -eq 0 ]; then
        echo "Converting $ACCESSION to FASTQ format..."
        fasterq-dump -3 --skip-technical "$ACCESSION"
        wait
        echo "Converting $ACCESSION to FASTQ format complete."
        #
        echo "Deleting $ACCESSION SRA files..."
        rm -r "$ACCESSION"
        echo "Compressing $ACCESSION..."
        pigz -p "$NCORES" "${ACCESSION}*fastq"
        wait
    else
        echo "Error downloading $ACCESSION. Skipping..."
    fi
done < "$ACCESSION_FILE"

echo "All downloads complete."
