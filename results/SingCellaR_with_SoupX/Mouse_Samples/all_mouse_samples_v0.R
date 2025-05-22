###################################
library(SingCellaR)
###################################
###################################
mouse <- new("SingCellaR")
mouse@dir_path_SingCellR_object_files<-"../Rcodes_SoupX/"
mouse@SingCellR_object_files=c("9T.SoupX.SingCellaR.rdata",
                               "Tr001.SoupX.SingCellaR.rdata",
                               "Tr015.SoupX.SingCellaR.rdata",
                               "Tr089.SoupX.SingCellaR.rdata",
                               "Tr091.SoupX.SingCellaR.rdata",
                               "Tr110.SoupX.SingCellaR.rdata",
                               "Tr074.SoupX.SingCellaR.rdata",
                               "Tr075.SoupX.SingCellaR.rdata",
                               "Tr085.SoupX.SingCellaR.rdata",
                               "Tr086.SoupX.SingCellaR.rdata",
                               "NSG_pons_cntrl.SoupX.SingCellaR.rdata",
                               "NSG_cortex_cntrl.SoupX.SingCellaR.rdata",
                               "CD1_pons_cntrl.SoupX.SingCellaR.rdata",
                               "CD1_cortex_cntrl.SoupX.SingCellaR.rdata")


preprocess_integration(mouse,mito_genes_start_with = "mt-")
mouse

filter_cells_and_genes(mouse,
                       min_UMIs=300,
                       max_UMIs=200000,
                       min_detected_genes=200,
                       max_detected_genes=15000,
                       max_percent_mito=20,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = FALSE)

##############################################################
meta.data<-mouse@meta.data
meta.data$sample_name<-""
meta.data$sample_name[meta.data$sampleID=="1_9T"]<-"K27M_ep_mouse"
meta.data$sample_name[meta.data$sampleID=="1_Tr001"]<-"K27M_CD1_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr015"]<-"K27M_CD1_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr089"]<-"MET_CD1_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr091"]<-"MET_CD1_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr110"]<-"MET_CD1_allograft_P3"
meta.data$sample_name[meta.data$sampleID=="1_Tr074"]<-"MET_NSG_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr075"]<-"MET_NSG_allograft_P1"
meta.data$sample_name[meta.data$sampleID=="1_Tr085"]<-"MET_NSG_allograft_P3"
meta.data$sample_name[meta.data$sampleID=="1_Tr086"]<-"MET_NSG_allograft_P3"
meta.data$sample_name[meta.data$sampleID=="1_NSG_pons_cntrl"]<-"NSG_pons_cntrl"
meta.data$sample_name[meta.data$sampleID=="1_NSG_cortex_cntrl"]<-"NSG_cortex_cntrl"
meta.data$sample_name[meta.data$sampleID=="1_CD1_pons_cntrl"]<-"CD1_pons_cntrl"
meta.data$sample_name[meta.data$sampleID=="1_CD1_cortex_cntrl"]<-"CD1_cortex_cntrl"


meta.data$sample_type<-""
meta.data$sample_type[meta.data$sampleID=="1_9T"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_Tr001"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_Tr015"]<-"K27M"
meta.data$sample_type[meta.data$sampleID=="1_Tr089"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr091"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr110"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr074"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr075"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr085"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_Tr086"]<-"MET"
meta.data$sample_type[meta.data$sampleID=="1_NSG_pons_cntrl"]<-"control"
meta.data$sample_type[meta.data$sampleID=="1_NSG_cortex_cntrl"]<-"control"
meta.data$sample_type[meta.data$sampleID=="1_CD1_pons_cntrl"]<-"control"
meta.data$sample_type[meta.data$sampleID=="1_CD1_cortex_cntrl"]<-"control"

##############################################################
meta.data$sequencer<-""
meta.data$sequencer[meta.data$sampleID=="1_9T"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr001"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr015"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr089"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr091"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr110"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr074"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr075"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_Tr085"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_Tr086"]<-"novaseq6000"
meta.data$sequencer[meta.data$sampleID=="1_NSG_pons_cntrl"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_NSG_cortex_cntrl"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_CD1_pons_cntrl"]<-"hiseq4000"
meta.data$sequencer[meta.data$sampleID=="1_CD1_cortex_cntrl"]<-"hiseq4000"
##############################################################
meta.data$tenx_version<-""

meta.data$tenx_version[meta.data$sampleID=="1_9T"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr001"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr015"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr089"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_Tr091"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr110"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_Tr074"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr075"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_Tr085"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_Tr086"]<-"v3"
meta.data$tenx_version[meta.data$sampleID=="1_NSG_pons_cntrl"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_NSG_cortex_cntrl"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_CD1_pons_cntrl"]<-"v2"
meta.data$tenx_version[meta.data$sampleID=="1_CD1_cortex_cntrl"]<-"v2"
##############################################################
mouse@meta.data<-meta.data
##############################################################
normalize_UMIs(mouse,use.scaled.factor = T)
remove_unwanted_confounders(mouse,residualModelFormulaStr="~UMI_count+percent_mito")

get_variable_genes_for_integrative_data_by_fitting_GLM_model(mouse,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(mouse,gmt.file = "../../../GeneSets/mouse.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(mouse)

runPCA(mouse,use.components=50,use.regressout.data = T)

plot_PCA_Elbowplot(mouse)


SingCellaR::runUMAP(mouse,n.dims.use = 30,n.neighbors = 30,uwot.metric = "euclidean")

plot_umap_label_by_a_feature_of_interest(mouse,feature = "sampleID",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_name",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "sample_type",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "sequencer",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(mouse,feature = "tenx_version",point.size = 0.5)

#SingCellaR::identifyClusters(mouse,n.dims.use = 30,n.neighbors = 30, 
#                             knn.metric = "euclidean")

#plot_umap_label_by_clusters(mouse,show_method = "louvain")

saveRDS(mouse,file="all_Mouse_samples.SoupX.SingCellaR_no_batch_correction.rds")
