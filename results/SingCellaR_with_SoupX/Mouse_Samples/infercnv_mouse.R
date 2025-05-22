library(infercnv)
library(SingCellaR)
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR_with_SoupX/Mouse_Samples/")

mouse <- readRDS("all_mouse_samples.SoupX.SingCellaR.Harmony.v1.rds")
# Step 1: Extract raw counts and annotation from SingCellaR object
expr_matrix <- get_normalized_umi(mouse)  # Replace with raw UMI counts if needed
annotation <- mouse@sc.clusters

# Step 2: Create the annotation file for inferCNV
# Assuming 'SingCellaR_anno_layer2' contains cell type labels
annotation_file <- annotation[, c("Cell",  "louvain_cluster")]
colnames(annotation_file) <- c("Cell", "Group")
write.table(annotation_file, "./cell_annotations.txt",  sep="\t", quote=FALSE, row.names=FALSE, col.names=F)

# Step 3: Gene order file for mouse genome
# Download a gene order file for mouse from inferCNV resources or create one.
# Ensure it contains columns: "Gene", "Chromosome", "Start", "End"
gene_order_file <- "~/Downloads/mouse_gene_position.txt"  # Replace with the path to your gene order file

refs = c('cl5','cl14','cl17','cl19','cl20','cl25','cl27')

# Step 4: Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expr_matrix,
  annotations_file = "./cell_annotations.txt",
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = refs  # Replace with the names of normal cells in your data
)

# Step 5: Run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  out_dir = "inferCNV_output",   # Output directory
  cutoff = 0.1,                 # Minimum expression threshold
  cluster_by_groups = TRUE,     # Cluster cells by annotation groups
  denoise = TRUE,               # Enable denoising
  HMM = TRUE, # Enable HMM to infer CNV events
  output_format = 'pdf',
  resume_mode = F                   
)

# Step 6: Visualize Results
# inferCNV generates heatmaps and other outputs in the specified output directory.
