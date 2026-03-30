#!/bin/bash
# blongpi - Step 2: Functional Search
# Requirements: hmmer

# Get the absolute path of the blongpi directory
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
HMM_DB="$BASE_DIR/data/hmms/markers.hmm"
FAA_FILES=$(ls $BASE_DIR/results/annotations/*/*.faa)
OUT_DIR="$BASE_DIR/results/functions"

mkdir -p $OUT_DIR

if [ ! -f "$HMM_DB" ]; then
    echo "Error: HMM database not found at $HMM_DB"
    echo "Please place your combined .hmm file there."
    exit 1
fi

echo ">>> Searching for markers using HMMER..."
for f in $FAA_FILES; do
    id=$(basename $f .faa)
    echo "Processing $id..."
    hmmsearch --cpu 4 --tblout $OUT_DIR/${id}_hits.txt $HMM_DB $f
done

# Summarize hits
echo ">>> Collating results into combined_markers.tsv..."
# Header
echo -e "genome_id\ttarget_name\tquery_name\tfull_sequence_e_value" > $OUT_DIR/combined_markers.tsv
# Extract hits and include genome ID from filename
for f in $OUT_DIR/*_hits.txt; do
    id=$(basename $f _hits.txt)
    grep -v "^#" "$f" | awk -v gid="$id" 'BEGIN {OFS="\t"} {print gid,$1,$3,$5}' >> $OUT_DIR/combined_markers.tsv
done
