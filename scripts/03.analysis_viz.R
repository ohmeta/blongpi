#!/usr/bin/env Rscript
# blongpi - Step 3: Statistical Analysis & Visualization

library(tidyverse)
library(ComplexHeatmap)
library(ggtree)
library(patchwork)
library(here)

# Paths
base_dir <- here("blongpi")
meta_path <- here("analysis/results/tables/asv_mag/b_longum_detailed_analysis.tsv")
pa_path <- file.path(base_dir, "results/pangenome/gene_presence_absence.csv")
tree_path <- file.path(base_dir, "results/phylogeny/blongum_core.treefile")
out_dir <- file.path(base_dir, "results")

# 1. Load and Clean Data
# We use your existing analysis table as metadata for MAG IDs
metadata_full <- read_tsv(meta_path) %>%
  # Extract group from mag_id if possible, or join with infant metadata
  # Here we take unique MAGs and assume the 'group' info is present or needs matching
  distinct(mag_id, .keep_all = TRUE)

# Ensure mag_id matches column names in Panaroo (strip extensions if needed)
metadata_full <- metadata_full %>%
  mutate(mag_id_clean = gsub(".fa.gz|.fa|.fasta", "", basename(mag_id)))

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

# Align metadata with columns
mag_cols <- colnames(gene_data)
# Filter metadata to match only present MAGs
metadata <- metadata_full %>% 
  filter(mag_id_clean %in% mag_cols)

# Define groups (Mock logic: if group is not in TSV, assign based on Sample name pattern)
if (!"group" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(group = ifelse(grepl("C$", isolate), "Control", "Exposure"))
}

# Subset gene_data to match metadata order
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
  p_tree <- ggtree(tree) %<+% metadata +
    geom_tippoint(aes(color = group), size=3) +
    geom_tiplab(size = 2, offset = 0.001) +
    scale_color_manual(values = c("Control" = "#3288BD", "Exposure" = "#D53E4F")) +
    theme(legend.position = "right") +
    ggtitle("Phylogenomic Clustering of B. longum MAGs")
  print(p_tree)
}

# Heatmap of key traits
target_genes <- c("tadA", "bsh", "fnb", "aga", "lacZ", "lnpA", "bln") # Common B. longum markers
viz_mat <- as.matrix(gene_data_sub[rownames(gene_data_sub) %in% target_genes, ])

if (nrow(viz_mat) > 0) {
  ha = HeatmapAnnotation(
    Group = metadata$group,
    # Use efficiency if available
    Efficiency = if("efficiency_exposure" %in% colnames(metadata)) metadata$efficiency_exposure else NULL,
    col = list(Group = c("Control" = "#3288BD", "Exposure" = "#D53E4F"))
  )

  h1 <- Heatmap(viz_mat, name = "Gene Presence", 
                top_annotation = ha,
                col = c("0" = "white", "1" = "steelblue"),
                row_names_side = "left",
                column_title = "B. longum: Targeted Functional Profiles")
  draw(h1)
}

dev.off()

write_tsv(enriched, file.path(out_dir, "enriched_genes.tsv"))
cat("Analysis complete. Results in 'blongpi/results/'.\n")
