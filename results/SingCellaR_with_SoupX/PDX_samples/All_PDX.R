###################################
library(SingCellaR)
###################################
###################################
PDX <- new("SingCellaR")
PDX@dir_path_SingCellR_object_files<-"../R_codes_SoupX/"
PDX@SingCellR_object_files=c("SJ_DIPG_X7_XP1.SoupX.SingCellaR.rdata",
                             "SU_DIPG_XIII.SoupX.SingCellaR.rdata",
                             "Astro110fh.SoupX.SingCellaR.rdata",
                             "RCXXB47_XP1_pellet.SoupX.SingCellaR.rdata",
                             "Tr041.SoupX.SingCellaR.rdata",
                             "Tr068.SoupX.SingCellaR.rdata",
                             "Tr082.SoupX.SingCellaR.rdata",
                             "Tr083.SoupX.SingCellaR.rdata")

preprocess_integration(PDX,mito_genes_start_with = "^MT-")
PDX

filter_cells_and_genes(PDX,
                       min_UMIs=500,
                       max_UMIs=200000,
                       min_detected_genes=200,
                       max_detected_genes=15000,
                       max_percent_mito=20,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = FALSE)

##############################################################
meta.data<-PDX@meta.data

meta.data$sample_name<-""
meta.data$sample_name[meta.data$sampleID=="1_SJ_DIPG_X7_XP1"]<-"K27M_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr041"]<-"K27M_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_SU_DIPG_XIII"]<-"K27M_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Astro110fh"]<-"K27M_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_RCXXB47_XP1_pellet"]<-"MET_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr068"]<-"MET_xenograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr082"]<-"MET_xenograft_P3"
meta.data$sample_name[meta.data$sampleID=="1_Tr083"]<-"MET_xenograft_P3"

meta.data$sample_type<-""
meta.data$sample_type[meta.data$sampleID=="1_SJ_DIPG_X7_XP1"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_Tr041"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_SU_DIPG_XIII"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_Astro110fh"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_RCXXB47_XP1_pellet"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr068"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr082"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr083"]<-"MET"
##############################################################
meta.data$sequencer<-""
meta.data$sequencer[meta.data$sampleID=="1_SJ_DIPG_X7_XP1"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr041"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_SU_DIPG_XIII"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Astro110fh"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_RCXXB47_XP1_pellet"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr068"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr082"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr083"]<-"novaseq6000"
##############################################################

meta.data$tenx_version<-""

meta.data$tenx_version[meta.data$sampleID=="1_SJ_DIPG_X7_XP1"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr041"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_SU_DIPG_XIII"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_Astro110fh"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_RCXXB47_XP1_pellet"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr068"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr082"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr083"]<-"v3"
##############################################################
PDX@meta.data<-meta.data
##############################################################
normalize_UMIs(PDX,use.scaled.factor = T)
remove_unwanted_confounders(PDX,residualModelFormulaStr="~UMI_count+percent_mito")

get_variable_genes_for_integrative_data_by_fitting_GLM_model(PDX,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(PDX,gmt.file = "../../../../KATI_et_al_project/GeneSets/human.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(PDX)

runPCA(PDX,use.components=50,use.regressout.data = T)

plot_PCA_Elbowplot(PDX)

SingCellaR::runUMAP(PDX,n.dims.use = 30,n.neighbors = 30,uwot.metric = "euclidean")

plot_umap_label_by_a_feature_of_interest(PDX,feature = "sampleID",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(PDX,feature = "sample_name",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(PDX,feature = "sequencer",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(PDX,feature = "tenx_version",point.size = 0.5)

SingCellaR::identifyClusters(PDX,n.dims.use = 30,
                             n.neighbors = 40, knn.metric = "euclidean")

plot_umap_label_by_clusters(PDX,show_method = "louvain")

saveRDS(PDX,file="PDX_SoupX.SingCellaR.No_batch_correction.rds")
##################NNLM##################################