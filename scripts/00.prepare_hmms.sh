#!/bin/bash
# blongpi - Step 0: HMM Preparation
# This script downloads targeted Pfam HMMs for B. longum functional analysis.

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
HMM_DIR="$BASE_DIR/data/hmms"
HMM_OUT="$HMM_DIR/markers.hmm"

mkdir -p $HMM_DIR

# Define PFAM IDs and their biological roles
# PF00437: TadA (Adhesion)
# PF07811: TadE (Adhesion)
# PF04964: Flp (Adhesion)
# PF04203: Sortase (Adhesion)
# PF02275: BSH (Stress)
# PF01120: Fucosidase GH29 (HMO)
# PF22124: Fucosidase GH95 (HMO)
# PF02973: Sialidase GH33 (HMO)
# PF17385: LNB-phosphorylase (HMO)
# PF00006: F1F0-ATPase alpha/beta (Stress)

IDs=("PF00437" "PF07811" "PF04964" "PF04203" "PF02275" "PF01120" "PF22124" "PF02973" "PF17385" "PF00006")

echo ">>> Cleaning old HMMs..."
rm -f $HMM_OUT $HMM_OUT.h3*

for id in "${IDs[@]}"; do
    echo "Processing $id..."
    # Download from InterPro API (using annotation=hmm parameter)
    curl -s -L "https://www.ebi.ac.uk/interpro/api/entry/pfam/$id?annotation=hmm" -o "$HMM_DIR/$id.hmm.gz"
    
    if file "$HMM_DIR/$id.hmm.gz" | grep -q "gzip compressed data"; then
        gunzip -c "$HMM_DIR/$id.hmm.gz" >> $HMM_OUT
        echo -e "\n" >> $HMM_OUT
        rm "$HMM_DIR/$id.hmm.gz"
    else
        echo "Error: Failed to download or invalid file for $id"
    fi
done

echo ">>> Pressing HMM database..."
hmmpress $HMM_OUT

echo ">>> HMM database prepared at $HMM_OUT"
