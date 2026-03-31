#!/bin/bash
# blongpi - Step 4: Pangenome GWAS (Scoary)
# Requirements: scoary

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
PA_FILE="$BASE_DIR/results/pangenome/gene_presence_absence.csv"
TRAIT_FILE="$BASE_DIR/data/metadata_traits.csv"
OUT_DIR="$BASE_DIR/results/pgwas"

mkdir -p $OUT_DIR

echo ">>> Preparing metadata for Scoary..."
# Scoary needs a simple CSV with: Name,Trait1,Trait2
# We extract this from your detailed analysis table
METADATA="/mnt/store/users/zhujie/HOLA/assay/results/tables/b_longum_detailed_analysis.tsv"

# Create Scoary trait file: Column 1 is ID, Column 2 is Group (Binary 0/1)
echo "f.id,is_exposure" > $TRAIT_FILE
tail -n +2 $METADATA | awk -F'\t' '{
    id=basename($1); 
    gsub(".fa.gz|.fa.fna|.fa|.fasta","",id); 
    group=($11=="Exposure"?1:0); 
    print id","group
}' >> $TRAIT_FILE

echo ">>> Running Scoary P-GWAS..."
# -g: pangenome matrix, -t: traits, -c: correction method
scoary -g $PA_FILE -t $TRAIT_FILE -o $OUT_DIR --p_value_cutoff 0.05 --correction bh

echo ">>> P-GWAS complete. Check results in $OUT_DIR"
