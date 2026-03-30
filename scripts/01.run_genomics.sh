#!/bin/bash
# blongpi - Step 1: Foundation
# Requirements: bakta, panaroo, iqtree

# Get the absolute path of the blongpi directory
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
MAG_DIR="$BASE_DIR/data/mags"
RES_DIR="$BASE_DIR/results"

mkdir -p $RES_DIR/annotations $RES_DIR/pangenome $RES_DIR/phylogeny

# 1. Annotation
echo ">>> Running Bakta..."
for f in $MAG_DIR/*.fa; do
    id=$(basename $f .fa)
    prokka --outdir $RES_DIR/annotations/$id --prefix $id --genus Bifidobacterium --species longum --cpus 32 $f
done

# 2. Pangenome (Panaroo)
echo ">>> Running Panaroo..."
# Collect all GFF files
GFFS=$(ls $RES_DIR/annotations/*/*.gff)
# Suppress Biopython warnings and adjust core threshold for MAGs
export PYTHONWARNINGS="ignore::BiopythonDeprecationWarning"
panaroo -i $GFFS -o $RES_DIR/pangenome --clean-mode strict -a core --core_threshold 0.5 --threads 32

# 3. Phylogeny
echo ">>> Running IQ-TREE..."
iqtree -s $RES_DIR/pangenome/core_gene_alignment.aln -m GTR+G -nt AUTO -bb 1000 -pre $RES_DIR/phylogeny/blongum_core
