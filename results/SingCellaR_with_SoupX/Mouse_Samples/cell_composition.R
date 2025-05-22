# Load required libraries
library(SingCellaR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(Seurat)
library(monocle3)
library(fgsea)
library(biomaRt)


# Set working directory
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR_with_SoupX/Mouse_Samples/")

# Load the processed data
mouse <- readRDS("./all_mouse_samples.SoupX.SingCellaR.Harmony.v1.rds")

###############################
# 1. Cross-Species Analysis
###############################
# Convert human genes to mouse homologs
gene_map <- read.table("~/Downloads/HOM_MouseHumanSequence.rpt.txt", sep="\t", header=TRUE)

# Separate Mouse and Human Genes
mouse_genes <- gene_map[gene_map$Common.Organism.Name == "mouse, laboratory", c("DB.Class.Key", "Symbol")]
human_genes <- gene_map[gene_map$Common.Organism.Name == "human", c("DB.Class.Key", "Symbol")]

# Merge on Homologous Cluster ID
gene_mapping <- merge(mouse_genes, human_genes, by = "DB.Class.Key", suffixes = c("_Mouse", "_Human"))

# Select Relevant Columns
gene_mapping <- gene_mapping[, c("Symbol_Mouse", "Symbol_Human")]

# Remove Duplicates
gene_mapping <- unique(gene_mapping)

# View First Few Mappings
head(gene_mapping)

  # Map human genes to mouse
m2 = mouse
mouse_gene_expression <- mouse@meta.data
rownames(mouse_gene_expression) <- gene_map$mgi_symbol[match(rownames(mouse_gene_expression), gene_map$hgnc_symbol)]

# Compare cell type distributions across species
p1 <- plot_umap_label_by_clusters(mouse) + ggtitle("Mouse Tumor Cells")
p2 <- plot_umap_label_by_clusters(human) + ggtitle("Human Tumor Cells")

ggsave("cross_species_umap.pdf", p1 + p2, width = 12, height = 6)

###############################
# 2. Cell Proportion Shifts Across Passages
###############################
cell_counts <- table(mouse@sc.clusters$Condition, mouse@sc.clusters$Cell_Type)
cell_proportions <- prop.table(cell_counts, margin = 2)

# Plot Sankey diagram
library(networkD3)
sankey_data <- as.data.frame(as.table(cell_proportions))
names(sankey_data) <- c("Condition", "Cell type", "Proportion")
write.table(sankey_data, "cell_proportions_sankey.tsv", sep = "\t")

###############################
# 3. Tumor-Associated Macrophage (TAM) Analysis
###############################
# Define M1 and M2 markers
M1_markers <- c("Tnf", "Il6", "Nos2")
M2_markers <- c("Cd163", "Mrc1", "Il10")

plot_violin_for_marker_genes(mouse,gene_list = M1_markers)
plot_violin_for_marker_genes(mouse,gene_list = M2_markers)
plot_umap_label_by_genes(mouse,gene_list = M1_markers)

# Compute M1/M2 ratio per sample
immune_cells <- subset(mouse@meta.data, rownames(mouse) %in% c(M1_markers, M2_markers))
M1_scores <- colMeans(immune_cells[M1_markers, ], na.rm = TRUE)
M2_scores <- colMeans(immune_cells[M2_markers, ], na.rm = TRUE)
M1_M2_ratio <- M1_scores / (M2_scores + 1e-6)

# Visualize M1 vs. M2 shift
m1m2_df <- data.frame(Sample = colnames(immune_cells), M1_M2_ratio)
ggplot(m1m2_df, aes(x = Sample, y = M1_M2_ratio)) +
  geom_bar(stat = "identity") +
  ggtitle("M1/M2 Ratio Across Mouse Models")
ggsave("M1_M2_ratio.pdf", width = 8, height = 6)

###############################
# 4. MES1 vs MES2 Tumor States
###############################
MES1_markers <- c("CD44", "VIM", "TGFBI")
MES2_markers <- c("HIF1A", "VEGFA", "CA9")

MES1_scores <- colMeans(mouse@normalized.umi[MES1_markers, ], na.rm = TRUE)
MES2_scores <- colMeans(mouse@normalized.umi[MES2_markers, ], na.rm = TRUE)

# Compare MES1 vs MES2 proportions
mes_df <- data.frame(MES1 = MES1_scores, MES2 = MES2_scores)
ggplot(mes_df, aes(x = rownames(mes_df), y = MES1)) +
  geom_bar(stat = "identity", fill = "blue") +
  ggtitle("MES1 vs MES2 Proportions")
ggsave("MES1_vs_MES2.pdf", width = 8, height = 6)

###############################
# 5. Immune Composition in CD1 vs NSG Models
###############################
immune_composition <- table(mouse@sc.clusters$louvain_cluster, mouse@sc.clusters$mouse_strain)
immune_proportions <- prop.table(immune_composition, margin = 2)

immune_df <- as.data.frame(as.table(immune_proportions))
names(immune_df) <- c("Cluster", "Strain", "Proportion")
ggplot(immune_df, aes(x = Strain, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Immune Cell Differences in CD1 vs NSG Mice")
ggsave("immune_composition_CD1_vs_NSG.pdf", width = 10, height = 6)
