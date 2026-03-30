#!/usr/bin/env Rscript
# blongpi - Step 3: Statistical Analysis & Visualization

library(tidyverse)
library(ComplexHeatmap)
library(ggtree)
library(patchwork)
library(here)

# Paths
base_dir <- "/mnt/store/users/zhujie/HOLA/assay/blongpi"
meta_path <- "/mnt/store/users/zhujie/HOLA/assay/results/tables/b_longum_detailed_analysis.tsv"
pa_path <- file.path(base_dir, "results/pangenome/gene_presence_absence.csv")
tree_path <- file.path(base_dir, "results/phylogeny/blongum_core.treefile")
out_dir <- file.path(base_dir, "results")

# 1. Load and Clean Data
# We use your existing analysis table as metadata for MAG IDs
metadata_full <- read_tsv(meta_path) %>%
  # Extract group from user_genome if possible, or join with infant metadata
  distinct(user_genome, .keep_all = TRUE)

# Ensure user_genome matches column names in Panaroo (strip extensions if needed)
metadata_full <- metadata_full %>%
  mutate(mag_id_clean = gsub(".fa.gz|.fa|.fasta", "", basename(user_genome)))

if (!file.exists(pa_path)) {
  stop(paste("Pangenome file not found at:", pa_path, ". Run Step 1 first."))
}

pa_matrix <- read_csv(pa_path)
tree <- if(file.exists(tree_path)) read.tree(tree_path) else NULL

# 2. Enrichment Analysis (Exposure vs Control)
# Extract the binary presence/absence matrix
gene_data <- pa_matrix %>%
  select(Gene, matches("HOLA|SRR")) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, 1)))

# Clean column names of gene_data to match metadata
colnames(gene_data) <- gsub(".fa.gz|.fa.fna|.fa|.fasta", "", basename(colnames(gene_data)))

# Align metadata with columns
mag_cols <- colnames(gene_data)
# Filter metadata to match only present MAGs and ensure mag_id_clean is the FIRST column for ggtree
metadata <- metadata_full %>% 
  filter(mag_id_clean %in% mag_cols) %>%
  select(mag_id_clean, everything()) %>%
  # Match the order of columns in gene_data
  arrange(match(mag_id_clean, mag_cols))

# Subset gene_data to match metadata order (only those present in metadata)
gene_data_sub <- gene_data[, metadata$mag_id_clean]

test_gene <- function(gene_row, groups) {
  tab <- table(factor(gene_row, levels=c(0,1)), groups)
  if(ncol(tab) < 2) return(NA)
  fisher.test(tab)$p.value
}

p_vals <- apply(gene_data_sub, 1, test_gene, groups = metadata$group)
enriched <- data.frame(Gene = names(p_vals), p_val = p_vals) %>%
  filter(!is.na(p_val)) %>%
  mutate(adj_p = p.adjust(p_val, method = "BH")) %>%
  filter(p_val < 0.05) %>%
  arrange(p_val)

# 3. Visualization
pdf(file.path(out_dir, "blongum_final_report.pdf"), width=14, height=10)

# Tree (if exists)
if (!is.null(tree)) {
  # Clean tree labels to match metadata
  tree$tip.label <- gsub(".fa.gz|.fa.fna|.fa|.fasta", "", basename(tree$tip.label))
  
  p_tree <- ggtree(tree) %<+% metadata +
    geom_tippoint(aes(color = group), size=3) +
    geom_tiplab(size = 2, offset = 0.001) +
    scale_color_manual(values = c("Control" = "#3288BD", "Exposure" = "#D53E4F")) +
    theme(legend.position = "right") +
    ggtitle("Phylogenomic Clustering of B. longum MAGs")
  print(p_tree)
}

# Targeted Markers from HMMER
hmm_path <- file.path(out_dir, "functions/combined_markers.tsv")
if (file.exists(hmm_path)) {
  hmm_data <- read_tsv(hmm_path) %>%
    mutate(genome_id_clean = gsub(".fa.gz|.fa.fna|.fa|.fasta", "", genome_id)) %>%
    # Filter for significant hits (e.g. E-value < 1e-5)
    filter(full_sequence_e_value < 1e-5) %>%
    # Map HMM names to biological markers
    mutate(marker = case_when(
      query_name == "T2SSE" ~ "tadA (Adhesion)",
      query_name == "TadE" ~ "tadE (Adhesion)",
      query_name == "Flp_Fap" ~ "flp (Adhesion)",
      query_name == "Sortase" ~ "sortase (Adhesion)",
      query_name == "CBAH" ~ "bsh (Stress)",
      query_name == "Alpha_L_fucos" ~ "fucosidase_GH29 (HMO)",
      query_name == "Glyco_hydro_95_cat" ~ "fucosidase_GH95 (HMO)",
      query_name == "Sialidase" ~ "sialidase_GH33 (HMO)",
      query_name == "LBP_M" ~ "lnpA (HMO)",
      query_name == "ATP-synt_ab" ~ "atpD (Stress)",
      TRUE ~ query_name
    )) %>%
    distinct(genome_id_clean, marker) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = genome_id_clean, values_from = present, values_fill = 0) %>%
    column_to_rownames("marker")

  # Align with metadata
  hmm_viz_mat <- as.matrix(hmm_data[, metadata$mag_id_clean, drop=FALSE])

  ha = HeatmapAnnotation(
    Group = metadata$group,
    col = list(Group = c("Control" = "#3288BD", "Exposure" = "#D53E4F"))
  )

  h2 <- Heatmap(hmm_viz_mat, name = "Targeted Marker", 
                top_annotation = ha,
                col = c("0" = "white", "1" = "darkgreen"),
                row_names_side = "left",
                column_title = "B. longum: Targeted Functional Markers (HMMER)")
  draw(h2)
}

# Heatmap of key traits from Panaroo (if needed)
# ...

dev.off()

write_tsv(enriched, file.path(out_dir, "enriched_genes.tsv"))
cat("Analysis complete. Results in 'blongpi/results/'.\n")
