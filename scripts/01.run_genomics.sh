#!/bin/bash
# blongpi - Step 1: Foundation
# Requirements: prokka, panaroo, iqtree

# Get the absolute path of the blongpi directory
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
MAG_DIR="$BASE_DIR/data/mags"
RES_DIR="$BASE_DIR/results"

mkdir -p $RES_DIR/annotations $RES_DIR/pangenome $RES_DIR/phylogeny

# 1. Annotation
echo ">>> Running Prokka..."
for f in $MAG_DIR/*.fasta; do
    id=$(basename $f .fasta)
    prokka --outdir $RES_DIR/annotations/$id --prefix $id --genus Bifidobacterium --species longum --cpus 8 $f
done

# 2. Pangenome (Panaroo)
echo ">>> Running Panaroo..."
# Collect all GFF files
GFFS=$(ls $RES_DIR/annotations/*/*.gff)
panaroo -i $GFFS -o $RES_DIR/pangenome --clean-mode strict -a core --threads 8

# 3. Phylogeny
echo ">>> Running IQ-TREE..."
iqtree -s $RES_DIR/pangenome/core_gene_alignment.aln -m GTR+G -nt AUTO -bb 1000 -pre $RES_DIR/phylogeny/blongum_core
