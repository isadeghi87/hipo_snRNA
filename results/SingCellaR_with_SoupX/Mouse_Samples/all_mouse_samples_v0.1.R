###################################
library(SingCellaR)
###################################
###################################
setwd("/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/hipo_temp/results/SingCellaR_with_SoupX/Mouse_Samples/")
memory.limit(size = 32000)
mouse <- readRDS(file="all_Mouse_samples.SoupX.SingCellaR_no_batch_correction.rds")

plot_umap_label_by_qc(mouse)

normalize_UMIs(mouse,use.scaled.factor = T)
remove_unwanted_confounders(mouse,
                            residualModelFormulaStr="~UMI_count+percent_mito")

runPCA(mouse,use.components=30,use.regressout.data = T)

plot_PCA_Elbowplot(mouse)

SrunHarmony <- function(object,n.dims.use=30,covariates=c("data_set"),harmony.theta = NULL,harmony.lambda = NULL,harmony.sigma = 0.1,harmony.nclust = NULL,
                        harmony.max.iter = 10,n.seed=1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  ######################################################
  orig.pca<-get_pca.result(object)$x[,1:n.dims.use]
  cell.meta<-get_cells_annotation(object)
  cell.meta<-subset(cell.meta,IsPassed==TRUE)
  ######################################################
  ###Run Harmony########################################
  set.seed(seed=n.seed)
  ######################################################
  harmony_embeddings <- harmony::RunHarmony(orig.pca, 
                                            cell.meta,  
                                            covariates,
                                            theta = harmony.theta,
                                            lambda = harmony.lambda,
                                            sigma = harmony.sigma,
                                            nclust = harmony.nclust,
                                            max_iter = harmony.max.iter,
                                            verbose=TRUE)
  ######################################################
  object@Harmony.embeddings<-harmony_embeddings
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Harmony analysis is done!.")
}

################Harmony2##################
SrunHarmony(mouse,covariates = c("sampleID","tenx_version","sequencer"),
            n.dims.use = 30,harmony.max.iter = 10,
            n.seed = 1,harmony.theta = c(3,2,2))

SingCellaR::runUMAP(mouse,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",n.dims.use = 30,
                    n.neighbors = 30,uwot.metric = "euclidean")

p1 = plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 1)
p2 = plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 1)
p3 = plot_umap_label_by_a_feature_of_interest(mouse,feature = "tenx_version",point.size = 1)
p4 = plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_type",point.size = 1)

pp = (p1+p2)/(p3+p4)
ggsave(filename = "./figures/umap_samples_soupx.pdf",plot = pp,width = 10,height = 8,units = 'in')

SingCellaR::identifyClusters(mouse,useIntegrativeEmbeddings = T,
                             integrative_method = "harmony", 
                             n.dims.use = 30,n.neighbors = 50, 
                             knn.metric = "euclidean")

plot_umap_label_by_clusters(mouse,show_method = "louvain")+Seurat::NoLegend()

plot_jaccard_similarity_among_clusters(mouse,cluster.type = "louvain")

findMarkerGenes(mouse,cluster.type = "louvain")

clusters = unique(mouse@sc.clusters$louvain_cluster)
markerdata = mouse@marker.genes$louvain
for (cl in clusters){
  subdat = markerdata[[cl]]
  genes = subdat %>% arrange(desc(log2FC)) %>% dplyr::select(Gene) %>% as.data.frame() 
  topmarkers = genes$Gene[1:20]
  
  print(paste0('cluster:',cl,"-markers: "))
  print(topmarkers)
}

#export_marker_genes_to_table(human,cluster.type = "louvain",
#                             n.TopGenes=40,min.log2FC=0.25,min.expFraction=0.2,
#                             write.to.file="All_human_Top40_marker_genes.txt")

plot_heatmap_for_marker_genes(mouse,cluster.type = "louvain",
                              n.TopGenes = 8,rowFont.size = 5)

saveRDS(mouse,file="all_mouse_samples.SoupX.SingCellaR.Harmony.v1.rds")
#############################################################
######################GSEA analysis##########################
mouse = readr::read_rds("all_mouse_samples.SoupX.SingCellaR.Harmony.v1.rds")
human_pre_rankedGenes<-identifyGSEAPrerankedGenes_for_all_clusters(human,
                                                                   cluster.type = "louvain")

fgsea_Results<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = human_pre_rankedGenes, 
                                                  eps = 0,nPermSimple = 10000,
                                                  gmt.file = "../../../GeneSets/Kati_signature_genes_human.gmt")

plot_heatmap_for_fGSEA_all_clusters(fgsea_Results,isApplyCutoff = T,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row = 10,
                                    adjusted_pval = 0.30,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)

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
pdf("../../../results/SingCellaR_with_SoupX/Mouse_Samples/figures/cell_annotation_heatmap_soupx.pdf",width = 12,
    height = 8)
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

dev.off()

write.table(filtered_data,'mouse_annotated_cell_typs.tsv',sep = '\t')

