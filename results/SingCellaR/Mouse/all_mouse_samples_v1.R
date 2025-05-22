#############################################
library(SingCellaR)
library(harmony)

## set working dir
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR")

mouse <- readRDS("mouse_SingCellaR_no_batch_correction.rds")

runPCA(mouse,use.components=40,use.regressout.data = T)

plot_PCA_Elbowplot(mouse)

################Harmony2##################
SingCellaR::runHarmony(mouse,covariates = c("sampleID","tenx_version","sequencer"),
                       n.dims.use = 40,harmony.max.iter = 10,
                       n.seed = 1,harmony.theta = c(3,2,2))

SingCellaR::runUMAP(mouse,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",n.dims.use = 40,
                    n.neighbors = 40,uwot.metric = "euclidean")

SingCellaR::runTSNE(mouse,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",n.dims.use = 40)

plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
plot_tsne_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "tenx_version",point.size = 1)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_type",point.size = 1)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "SingCellaR_anno_layer1",point.size = 1)

plot_umap_label_by_a_feature_of_interest(mouse,feature = "SingCellaR_anno_layer2",point.size = 1)
plot_tsne_label_by_a_feature_of_interest(mouse,feature = "SingCellaR_anno_layer2",point.size = 1)

SingCellaR::identifyClusters(mouse,useIntegrativeEmbeddings = T,
                             integrative_method = "harmony", 
                             n.dims.use = 40,n.neighbors = 30, 
                             knn.metric = "euclidean")

plot_tsne_label_by_clusters(mouse,show_method = "louvain")
plot_umap_label_by_clusters(mouse,show_method = "louvain")

plot_jaccard_similarity_among_clusters(mouse,cluster.type = "louvain")

findMarkerGenes(mouse,cluster.type = "louvain")

export_marker_genes_to_table(mouse,cluster.type = "louvain",
                             n.TopGenes=40,min.log2FC=0.25,min.expFraction=0.2,
                             write.to.file="All_mouse_Top40_marker_genes.txt")

plot_heatmap_for_marker_genes(mouse,cluster.type = "louvain",
                              n.TopGenes = 10,rowFont.size = 5)

######################GSEA analysis##########################
mouse_pre_rankedGenes<-identifyGSEAPrerankedGenes_for_all_clusters(mouse,
                                                                   cluster.type = "louvain")

fgsea_Results <- Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = mouse_pre_rankedGenes, 
                                                  eps = 0,nPermSimple = 10000,
                                                  gmt.file = "../../../datasets/mouse_all_gene_set.gmt")

# plot_heatmap_for_fGSEA_all_clusters(fGSEA_results.data.frame = fgsea_Results,
#                                     isApplyCutoff = FALSE,
#                                     use_pvalues_for_clustering=T,
#                                     show_NES_score = T,fontsize_row = 10,
#                                     adjusted_pval = 0.30,
#                                     show_only_NES_positive_score = T,
#                                     format.digits = 3,
#                                     clustering_method = "ward.D",
#                                     clustering_distance_rows = "euclidean",
#                                     clustering_distance_cols = "euclidean",show_text_for_ns = F)


# Step 1: Filter data (if needed)
filtered_data <- as.data.frame(fgsea_Results)# Example: Filter by adjusted p-value threshold

library(tibble)
library(tidyverse)

# Step 2: Create a matrix for the heatmap
heatmap_matrix <- filtered_data %>%
  select(pathway, cluster, NES) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES,
    values_fill = NA,
    values_fn = mean  # Summarize duplicates by their mean NES values
  ) %>%
  column_to_rownames("pathway")

library(ComplexHeatmap)
# Replace NULL with NA
heatmap_matrix[heatmap_matrix == "NULL"] <- NA


# Convert to matrix
heatmap_matrix[is.na(heatmap_matrix)] = 0
heatmap_matrix <- as.matrix(heatmap_matrix)

# Verify the matrix
print(heatmap_matrix)

# Plot the heatmap
library(ComplexHeatmap)
Heatmap(
  heatmap_matrix,
  name = "NES",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D",
  clustering_distance_columns = "euclidean",
  clustering_distance_rows = "euclidean",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  #heatmap_legend_param = list(title = "Normalized Enrichment Score"),
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

saveRDS(mouse,file="mouse_all_samples_SingCellaR.Harmony.v1.rds")
mouse = readRDS(file = "mouse_all_samples_SingCellaR.Harmony.v1.rds")

write.table(filtered_data,'mouse_annotated_cell_typs.tsv',sep = '\t')

################CCA######################
##########################################
#meta.data<-get_cells_annotation(mouse)
#rownames(meta.data)<-meta.data$Cell
#
#runSeuratIntegration(mouse,Seurat.metadata=meta.data,n.dims.use = 20,
#                     Seurat.split.by = "data_set",use.SingCellaR.varGenes = FALSE)
#
#SingCellaR::runUMAP(mouse,useIntegrativeEmbeddings = T, 
#                    integrative_method = "seurat",n.dims.use = 20,
#                    n.neighbors = 30,uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)

#################rpca###################
########################################
#library(Seurat)
#
#meta.data<-get_cells_annotation(mouse)
#rownames(meta.data)<-meta.data$Cell
#
#runSeuratIntegration_with_rpca(mouse,Seurat.metadata=meta.data,n.dims.use = 20,
#                     Seurat.split.by = "sample_name",use.SingCellaR.varGenes = FALSE)
#
#SingCellaR::runUMAP(mouse,useIntegrativeEmbeddings = T, 
#                    integrative_method = "seurat",n.dims.use = 20,
#                    n.neighbors = 30,uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)
#
#
####################combat##############
#library(sva)
#
#runCombat(mouse,use.reduced_dim = F,batch_identifier = c("sampleID"))
#
#runPCA(mouse,use.regressout.data = T)
#
#SingCellaR::runUMAP(mouse,dim_reduction_method = "pca",
#                    n.dims.use = 30, n.neighbors = 30,
#                    uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)
#
####################Limma###############
#remove_unwanted_confounders(mouse,residualModelFormulaStr="~UMI_count+percent_mito+sampleID")
#runPCA(mouse,use.regressout.data = T)
#SingCellaR::runUMAP(mouse,dim_reduction_method = "pca",
#                    n.dims.use = 20, n.neighbors = 20,
#                    uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)

