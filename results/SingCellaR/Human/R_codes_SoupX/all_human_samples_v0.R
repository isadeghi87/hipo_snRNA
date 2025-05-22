###################################
library(SingCellaR)
###################################
K27M<-readRDS("../R_codes/K27M.SingCellaR.v3.rds")
MET<-readRDS("../R_codes/MET.louvain.v2.rds")
###################################
x.meta<-rbind(K27M@meta.data,MET@meta.data)
rownames(x.meta)<-x.meta$Cell
###################################
human <- new("SingCellaR")
human@dir_path_SingCellR_object_files<-"../R_codes_SoupX/"
human@SingCellR_object_files=c("ICGC_GBM27.SingCellaR.rdata",
                               "ICGC_GBM60.SingCellaR.rdata",
                               "ICGC_GBM96.SingCellaR.rdata",
                               "ICGC_GBM15.SingCellaR.rdata",
                               "ICGC_GBM71.SingCellaR.rdata")

preprocess_integration(human)
human

filter_cells_and_genes(human,
                       min_UMIs=400,
                       max_UMIs=100000,
                       min_detected_genes=300,
                       max_detected_genes=15000,
                       max_percent_mito=3,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = FALSE)

##############################################################
meta.data<-human@meta.data

meta.data$sample_name<-""
meta.data$sample_name[meta.data$sampleID=="1_ICGC_GBM27"]<-"K27M_p1_GBM27"
meta.data$sample_name[meta.data$sampleID=="1_ICGC_GBM60"]<-"K27M_p2_GBM60"
meta.data$sample_name[meta.data$sampleID=="1_ICGC_GBM96"]<-"K27M_p3_GBM96"
meta.data$sample_name[meta.data$sampleID=="1_ICGC_GBM15"]<-"MET_p2_GBM15"
meta.data$sample_name[meta.data$sampleID=="1_ICGC_GBM71"]<-"MET_p2_GBM71"

meta.data$sample_type<-""
meta.data$sample_type[meta.data$sampleID=="1_ICGC_GBM27"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_ICGC_GBM60"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_ICGC_GBM96"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_ICGC_GBM15"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_ICGC_GBM71"]<-"MET"
##############################################################
meta.data$sequencer<-""
meta.data$sequencer[meta.data$sampleID=="1_ICGC_GBM27"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_ICGC_GBM60"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_ICGC_GBM96"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_ICGC_GBM15"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_ICGC_GBM71"]<-"novaseq6000"
##############################################################
meta.data$tenx_version<-""

meta.data$tenx_version[meta.data$sampleID=="1_ICGC_GBM27"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_ICGC_GBM60"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_ICGC_GBM96"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_ICGC_GBM15"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_ICGC_GBM71"]<-"v2"
##############################################################
x2.meta<-x.meta[as.character(meta.data$Cell),]
meta.data$SingCellaR_anno_layer1<-x2.meta$SingCellaR_anno_layer1
meta.data$SingCellaR_anno_layer2<-x2.meta$SingCellaR_anno_layer2

human@meta.data<-meta.data
##############################################################
normalize_UMIs(human,use.scaled.factor = T)
remove_unwanted_confounders(human,residualModelFormulaStr="~UMI_count")

get_variable_genes_for_integrative_data_by_fitting_GLM_model(human,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(human,gmt.file = "../../../GeneSets/human.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(human)

runPCA(human,use.components=50,use.regressout.data = T)

plot_PCA_Elbowplot(human)


SingCellaR::runUMAP(human,n.dims.use = 40,n.neighbors = 30,uwot.metric = "euclidean")

plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human,feature = "tenx_version",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human,feature = "SingCellaR_anno_layer1",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human,feature = "SingCellaR_anno_layer2",point.size = 0.5)

#SingCellaR::identifyClusters(human,n.dims.use = 40,n.neighbors = 30, 
#                             knn.metric = "euclidean")
#
#plot_umap_label_by_clusters(human,show_method = "louvain")

saveRDS(human,file="Human_SingCellaR_no_batch_correction.rds")
##################NNLM###################################
#library(NNLM)

#SingCellaR::runNNMF(human)
#
#SingCellaR::runUMAP(human,dim_reduction_method = "nnmf",n.dims.use = 30,n.neighbors = 30,
#                    uwot.metric = "euclidean")
#plot_umap_label_by_a_feature_of_interest(human,feature = "sampleID",point.size = 0.5)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_type",point.size = 0.5)
#plot_umap_label_by_a_feature_of_interest(human,feature = "tenx_version",point.size = 0.2)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 0.2)

##################HARMONY################################
#runHarmony(human,covariates = c("sampleID"),n.dims.use = 50,
#           harmony.max.iter = 20,n.seed = 1,harmony.theta = 3)
#
#SingCellaR::runUMAP(human,useIntegrativeEmbeddings = T, integrative_method = "harmony",
#                    n.dims.use = 50,n.neighbors = 30,uwot.metric = "euclidean")

#plot_umap_label_by_a_feature_of_interest(human,feature = "sampleID",point.size = 0.5)
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_type",point.size = 0.5)

####################Seurat################################
#library(Seurat)

#meta.data<-get_cells_annotation(human)
#rownames(meta.data)<-meta.data$Cell

#runSeuratIntegration_with_rpca(human,Seurat.metadata=meta.data,n.dims.use = 20,
#                     Seurat.split.by = "data_set",use.SingCellaR.varGenes = F)

#SingCellaR::runUMAP(human,n.dims.use = 20,n.neighbors = 30,uwot.metric = "euclidean")
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 2)
#
#plot_umap_label_by_a_feature_of_interest(human,feature = "UMI_count",
#                                         point.size = 0.1,mark.feature = F)
#
#identifyClusters(human,n.dims.use = 20,n.neighbors = 30,knn.metric = "euclidean")
#
#plot_umap_label_by_clusters(human,show_method = "louvain",point.size = 1.5)

#findMarkerGenes(human,cluster.type = "louvain")
#
#plot_heatmap_for_marker_genes(human,cluster.type = "louvain",
#                              n.TopGenes = 8,rowFont.size = 5)
#
#pre_rankedGenes_for_GSEA<-identifyGSEAPrerankedGenes_for_all_clusters(human,
#                                                                      cluster.type = "louvain")
#
#save(pre_rankedGenes_for_GSEA,file="all_human_preRankedGenes_for_GSEA.rdata")

#load("all_human_preRankedGenes_for_GSEA.rdata")
#fgsea_Results<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA, eps = 0,
#                                                  gmt.file = "../../../GeneSets/Kati_signature_genes_human.gmt")
#
#plot_heatmap_for_fGSEA_all_clusters(fgsea_Results,isApplyCutoff = TRUE,
#                                    use_pvalues_for_clustering=T,
#                                    show_NES_score = F,fontsize_row = 10,
#                                    adjusted_pval = 0.1,
#                                    show_only_NES_positive_score = T,format.digits = 3,
#                                    clustering_method = "ward.D",
#                                    clustering_distance_rows = "euclidean",
#                                    clustering_distance_cols = "euclidean",
#                                    show_text_for_ns = F)



#save(human,file="Human_samples_singCellaR.rdata")
