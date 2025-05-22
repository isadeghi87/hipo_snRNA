library(SingCellaR)
library(harmony)
library(infercnv)

## set working dir
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR/Mouse/")
mouse <- readRDS("./mouse_all_samples_SingCellaR.Harmony.v1.Rds")


# Step 1: Prepare the input data for inferCNV
# Extract raw counts and metadata from SingCellaR object
raw_counts <- mouse@raw.data  # Extract raw counts
cell_annotations <- get_cells_annotation(mouse)  # Extract cell annotations

# Ensure cell barcodes are consistent
rownames(cell_annotations) <- cell_annotations$Cell

# Create annotation file for inferCNV
# This file should map each cell to a cell type or state (e.g., tumor vs. normal)
annotation_file <- cell_annotations[, c("Cell", "SingCellaR_anno_layer1")]  # Replace "SingCellaR_anno_layer1" with the appropriate column for annotation
colnames(annotation_file) <- c("Cell", "Annotation")
write.table(annotation_file, "cell_annotations.txt", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)

# Create a gene ordering file for inferCNV
# Map genes to their respective chromosomes and positions
# If not already available, download or generate this information (e.g., using a GTF file or reference genome)
# Example format: GeneID, Chromosome, Start, End
gene_order_file <- "mouse_gene_order.txt"  # Provide the path to your gene ordering file

# Step 2: Create the inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,  # Raw expression matrix
  annotations_file = "cell_annotations.txt",  # Annotation file
  delim = "\t",  # Tab-separated annotation file
  gene_order_file = gene_order_file,  # Gene order file
  ref_group_names = c("Normal")  # Define reference group (e.g., normal cells)
)

# Step 3: Run inferCNV
infercnv_results <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,  # Minimum average expression threshold (adjust based on data)
  out_dir = "./inferCNV_results",  # Output directory
  cluster_by_groups = TRUE,  # Cluster cells by groups
  denoise = TRUE,  # Enable denoising
  HMM = TRUE,  # Enable Hidden Markov Model for CNV calling
  analysis_mode = "samples"  # Use sample-level analysis
)

# Step 4: Visualize inferCNV results
# The results folder will contain heatmaps and CNV predictions
# You can explore and visualize the CNV data with the provided outputs
