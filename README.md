# blongpi - Bifidobacterium longum Comparative Genomics Pipeline

This project compares 28 B. longum MAGs between Control and Exposure infants.

## Workflow:
1. Place MAG fasta files in `data/mags/`.
2. Run `scripts/01.run_genomics.sh` for annotation and pangenome.
3. Run `scripts/02.search_functions.sh` for targeted HMM searches.
4. Run `scripts/03.analysis_viz.R` for statistical enrichment and visualization.

## Tasks Covered:
- Task 1: Pangenome enrichment (Exposure vs Control).
- Task 2: Phylogenomics & Clustering.
- Task 3: Metabolic potential (HMO/Adhesion).
- Task 4: Mobile Genetic Elements.
- Task 5: Correlation with ASV sharing metrics.
