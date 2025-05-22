# Set working directory
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR_with_SoupX/Mouse_Samples/")

# Load the processed data
mouse <- readRDS("./all_mouse_samples.SoupX.SingCellaR.Harmony.annotated.Rds")


# Run differential expression analysis
tumor_vs_normal_DE <- FindMarkers(
  mouse,
  ident.1 = "Tumor",
  ident.2 = "Normal",
  group.by = "Condition",
  test.use = "wilcox"
)

# Save results
write.csv(tumor_vs_normal_DE, "DEGs_Tumor_vs_Normal.csv")

# View top differentially expressed genes
head(tumor_vs_normal_DE)