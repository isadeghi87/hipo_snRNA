#############################################
library(SingCellaR)
library(harmony)

human<-readRDS(file="Human_SingCellaR_no_batch_correction.rds")

normalize_UMIs(human,use.scaled.factor = T)
remove_unwanted_confounders(human,
                            residualModelFormulaStr="~UMI_count")

runPCA(human,use.components=40,use.regressout.data = T)

plot_PCA_Elbowplot(human)

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
SrunHarmony(human,covariates = c("sampleID","tenx_version","sequencer"),
                       n.dims.use = 40,harmony.max.iter = 10,
                       n.seed = 1,harmony.theta = c(3,2,2))

SingCellaR::runUMAP(human,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",n.dims.use = 40,
                    n.neighbors = 30,uwot.metric = "euclidean")

SingCellaR::runTSNE(human,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",
                    n.dims.use = 40)

#runFA2_ForceDirectedGraph(human,n.dims.use = 40,
#                          n.neighbors = 5,n.seed = 1,
#                          fa2_n_iter = 1000)
#SingCellaR::runKNN_Graph(human,n.dims.use = 40,n.neighbors = 10)
#
#library(threejs)
#plot_3D_knn_graph_label_by_clusters(human,show_method = "louvain",vertext.size = 0.5)
plot_umap_label_by_qc(human)

plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
plot_tsne_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "tenx_version",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "sample_type",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "SingCellaR_anno_layer1",point.size = 1)

plot_umap_label_by_a_feature_of_interest(human,feature = "SingCellaR_anno_layer2",point.size = 1)
plot_tsne_label_by_a_feature_of_interest(human,feature = "SingCellaR_anno_layer2",point.size = 1)

SingCellaR::identifyClusters(human,useIntegrativeEmbeddings = T,
                             integrative_method = "harmony", 
                             n.dims.use = 40,n.neighbors = 40, 
                             knn.metric = "euclidean")

#plot_tsne_label_by_clusters(human,show_method = "louvain")
plot_umap_label_by_clusters(human,show_method = "louvain")

plot_jaccard_similarity_among_clusters(human,cluster.type = "louvain")

findMarkerGenes(human,cluster.type = "louvain")

#export_marker_genes_to_table(human,cluster.type = "louvain",
#                             n.TopGenes=40,min.log2FC=0.25,min.expFraction=0.2,
#                             write.to.file="All_human_Top40_marker_genes.txt")

plot_heatmap_for_marker_genes(human,cluster.type = "louvain",
                              n.TopGenes = 10,rowFont.size = 5)

######################GSEA analysis##########################
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

saveRDS(human,file="Human_all_samples_SingCellaR.Harmony.v1.rds")
################CCA######################
##########################################
#meta.data<-get_cells_annotation(human)
#rownames(meta.data)<-meta.data$Cell
#
#runSeuratIntegration(human,Seurat.metadata=meta.data,n.dims.use = 20,
#                     Seurat.split.by = "data_set",use.SingCellaR.varGenes = FALSE)
#
#SingCellaR::runUMAP(human,useIntegrativeEmbeddings = T, 
#                    integrative_method = "seurat",n.dims.use = 20,
#                    n.neighbors = 30,uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)

#################rpca###################
########################################
#library(Seurat)
#
#meta.data<-get_cells_annotation(human)
#rownames(meta.data)<-meta.data$Cell
#
#runSeuratIntegration_with_rpca(human,Seurat.metadata=meta.data,n.dims.use = 20,
#                     Seurat.split.by = "sample_name",use.SingCellaR.varGenes = FALSE)
#
#SingCellaR::runUMAP(human,useIntegrativeEmbeddings = T, 
#                    integrative_method = "seurat",n.dims.use = 20,
#                    n.neighbors = 30,uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)
#
#
####################combat##############
#library(sva)
#
#runCombat(human,use.reduced_dim = F,batch_identifier = c("sampleID"))
#
#runPCA(human,use.regressout.data = T)
#
#SingCellaR::runUMAP(human,dim_reduction_method = "pca",
#                    n.dims.use = 30, n.neighbors = 30,
#                    uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)
#
####################Limma###############
#remove_unwanted_confounders(human,residualModelFormulaStr="~UMI_count+percent_mito+sampleID")
#runPCA(human,use.regressout.data = T)
#SingCellaR::runUMAP(human,dim_reduction_method = "pca",
#                    n.dims.use = 20, n.neighbors = 20,
#                    uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)

