#!/bin/bash
# blongpi - Step 6: Evolutionary Context (Cell 2026 Reference Fetcher)
# Logic: Identify representative clades from the Global Atlas.

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
REF_DIR="$BASE_DIR/data/references"
mkdir -p $REF_DIR

echo ">>> INSTRUCTIONS FOR CONTEXTUAL ANALYSIS:"
echo "1. Download Supplementary Table S1 from https://doi.org/10.1016/j.cell.2026.01.007"
echo "2. Filter for: 'B. infantis clades' and 'Shared strains'."
echo "3. Save the Accession numbers to data/ref_accessions.txt"
echo ""
echo ">>> This script will then use NCBI datasets to fetch genomes."

# Check if accession list exists
if [ -f "$BASE_DIR/data/ref_accessions.txt" ]; then
    echo ">>> Fetching Reference Genomes from NCBI..."
    # Requirements: ncbi-datasets CLI
    datasets download genome accession --inputfile $BASE_DIR/data/ref_accessions.txt --include genome --filename $REF_DIR/cell_refs.zip
    unzip $REF_DIR/cell_refs.zip -d $REF_DIR/
    echo ">>> References ready. You can now re-run Step 1 including these genomes."
else
    echo ">>> SKIP: data/ref_accessions.txt not found. Please provide accession list to proceed."
fi
