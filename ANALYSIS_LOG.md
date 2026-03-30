# Comprehensive Analysis Log & Technical Record: blongpi

This document serves as both a technical report and a pedagogical record of the troubleshooting process used to analyze 28 *Bifidobacterium longum* MAGs.

---

## 1. Executive Summary of Biological Principles

### 1.1 The Research Question
**Objective**: Identify genomic drivers of "Body-Site Sharing" (Nose ↔ Gut translocation) in infants.
**Meaningfulness**: *B. longum* is a keystone species. Understanding why some strains transition between body sites better than others informs our knowledge of infant health and probiotic colonization.

### 1.2 The Dual-Track Study Design
1.  **Discovery (Pangenome-wide)**: Using **Panaroo** to identify any gene enriched in the "Exposure" group. This is unbiased and can find novel genes.
2.  **Validation (Targeted HMM)**: Using **HMMER** to look for specific "smoking gun" domains related to adhesion, HMO metabolism, and stress survival.

### 1.3 Tool Principles
- **Panaroo**: Uses a graph-based approach to "repair" fragmented MAGs. It looks at the flanking genes of an assembly break to decide if a gene is actually missing or just failed to assemble.
- **HMMER**: Uses Profile Hidden Markov Models. It is far more sensitive than BLAST because it searches for the "evolutionary signature" of a protein family rather than just exact sequence matches.
- **IQ-TREE**: Uses Maximum Likelihood to construct a tree from core genes, allowing us to see if "sharing" strains are closely related (clonal) or if they independently acquired sharing genes (convergence).

---

## 2. Troubleshooting & Session Narrative (Learning Record)

### Phase 1: The Pangenome Failure (Step 1)
*   **Problem**: Panaroo reported "No gene clusters present above the core frequency threshold."
*   **Investigation**: We checked MAG sizes. One MAG was only 928KB. Standard *B. longum* is ~2.4MB.
*   **Lesson**: Default thresholds (95%) assume high-quality isolate genomes. With Metagenome-Assembled Genomes (MAGs), fragmentation is expected.
*   **Solution**: Lowered `--core_threshold` to `0.5`. This ensures that genes present in at least 50% of the genomes (14/28) are captured for the phylogenomic tree construction.
*   **How the tree was built**: By lowering the threshold to 50%, Panaroo could generate a `core_gene_alignment.aln`. For MAGs missing a specific gene, Panaroo fills the alignment with gaps (`---`). IQ-TREE then uses the concatenated sequence of these "pseudo-core" genes to calculate the phylogeny, allowing us to build a tree even with fragmented data.
*   **Fix**: Modified `scripts/01.run_genomics.sh`.
*   **Note on Summary Statistics**: Even after this fix, `summary_statistics.txt` will still show **0 Core genes**. This is because Panaroo's summary file uses fixed reporting bins (99% and 95%). Our `--core_threshold 0.5` parameter affects the **processing logic** (which genes are included in the alignment) rather than the **reporting labels**.
*   **Why aligned genes (1550) < shell genes (2104)?**: The "Shell" bin in the summary includes every gene present in 15% to 95% of genomes. However, our tree only uses genes present in >= 50% of genomes. Thus, many "Shell" genes that were too rare (e.g., present in only 20% of genomes) were correctly excluded from the phylogenetic alignment to ensure tree accuracy.

### Phase 2: The HMM Preparation (Step 2)
*   **Problem**: `markers.hmm` was missing, and initial strategy documents had incorrect IDs.
*   **Investigation**: `PF01391` was listed as TadA but returned "Collagen."
*   **Lesson**: Bioinformatic IDs change or can be mislabeled in strategy drafts. Always verify against the latest InterPro/Pfam API (2026).
*   **Solution**: Systematically verified every ID. Corrected TadA to `PF00437`, BSH to `PF02275`, etc.
*   **Fix**: Created `scripts/00.prepare_hmms.sh` to automate the download and indexing.

### Phase 3: The Functional Search & Collation
*   **Problem**: The initial search script didn't associate hits with specific genomes in the final table.
*   **Investigation**: `awk` was grabbing columns but the genome ID was buried in the filename.
*   **Lesson**: Data collation is the "bridge" between raw tools and R analysis. If the IDs don't match, the R script fails.
*   **Fix**: Updated `scripts/02.search_functions.sh` to loop through files and prepend the `genome_id`.

### Phase 4: The R Visualization Crash (Step 3)
*   **Problem**: `scripts/03.analysis_viz.R` failed with "at least two distinct break values" and missing group info.
*   **Investigation**: 
    1. The "Efficiency" metadata column was 80% `NA`s, crashing the heatmap scale.
    2. `ggtree` failed to color tips because the metadata ID column wasn't the first column.
    3. Paths were hardcoded to a different directory structure.
*   **Lesson**: Data hygiene is critical. Standardize IDs (strip `.fa.gz`, `.fna`) across all inputs (Pangenome, HMM, Metadata).
*   **Fix**: 
    - Removed the continuous "Efficiency" scale.
    - Reordered metadata columns (`select(mag_id_clean, everything())`).
    - Standardized string cleaning for IDs.

---

## 3. Final Technical Configuration

### Core Gene Definition
- **Threshold**: 0.5 (adjusted for MAGs).
- **Tool**: Panaroo (strict mode).

### Targeted Functional Markers (The HMM Suite)
| Marker | Pfam ID | Hypothesis |
| :--- | :--- | :--- |
| TadA | PF00437 | Adhesion (Pili Motor) |
| TadE | PF07811 | Adhesion (Minor Pilin) |
| Flp | PF04964 | Adhesion (Major Pilin) |
| Sortase | PF04203 | Adhesion (Cell Wall Anchor) |
| BSH | PF02275 | Stress (Bile Salt Survival) |
| GH29 | PF01120 | Metabolism (Fucosidase) |
| GH95 | PF22124 | Metabolism (Fucosidase) |
| GH33 | PF02973 | Metabolism (Sialidase) |
| lnpA | PF17385 | Metabolism (LNB-phosphorylase) |
| atpD | PF00006 | Stress (Acid resistance) |

---

## 4. The Relationship Between Discovery (Panaroo) and Validation (HMMER)

One might ask: *Were the targeted HMMs identified FROM the Panaroo results?*

The answer is **No**. In this pipeline, they are **Complementary but Independent** tracks. This "Dual-Track" approach is used to prevent bias:

| Feature | Panaroo (Bottom-Up) | HMMER (Top-Down) |
| :--- | :--- | :--- |
| **Logic** | "Show me every gene that is different." | "Look specifically for Tad Pili and HMO genes." |
| **Bias** | **Unbiased**: Can discover novel genes with no known name. | **Biased**: Only finds what we tell it to look for. |
| **Sensitivity** | **Lower**: Relies on general annotation (Bakta/Prokka). | **Higher**: Uses Profile HMMs to find "hidden" signatures. |
| **Role** | **Discovery**: Finds the "What." | **Validation**: Confirms the "Why" (Mechanism). |

### The "Convergence" Goal
The R script (`03.analysis_viz.R`) brings these together. We perform a pangenome-wide Fisher's test to find *statistically enriched clusters*. We then compare these clusters to our *Targeted HMM hits*. 
- If a Panaroo cluster is enriched AND it hits an HMM marker, we have **high-confidence evidence** of a biological driver.
- If an HMM marker is present but Panaroo doesn't show it as "enriched," it means the gene is likely a "core" trait of *B. longum* and not the reason why some strains "share" better than others.

---

## 5. Phylogenomic Interpretation: Clonal vs. Convergent

We build a tree using >50% genes to establish the "Null Model" of the species evolution. 

1.  **Scenario A (Clustering)**: If Exposure MAGs cluster together, the "Sharing" phenotype is **Clonal**. A specific lineage of *B. longum* is responsible.
2.  **Scenario B (Mixing)**: If Exposure and Control MAGs are mixed, but only Exposure MAGs have specific functional genes (from HMMER), the phenotype is **Convergent**. Different strains are independently adapting to the "Exposure" environment, potentially via Horizontal Gene Transfer (HGT).

Building the tree allows us to distinguish between these two fundamental biological processes.

---

## 6. Future Improvements (Comprehensive Analysis)
1.  **Operon Context**: Verify if the *entire* `tad` locus is present, not just single genes.
2.  **CAZyme Profiling**: Run **dbCAN3** to get a full map of sugar utilization.
3.  **MGE Detection**: Use **IslandPath** to see if sharing genes are on genomic islands.

---
**Workflow Execution Order**:
1. `bash scripts/00.prepare_hmms.sh`
2. `bash scripts/01.run_genomics.sh`
3. `bash scripts/02.search_functions.sh`
4. `Rscript scripts/03.analysis_viz.R`
