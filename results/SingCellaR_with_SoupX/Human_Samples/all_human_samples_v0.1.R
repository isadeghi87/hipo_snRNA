###################################
library(SingCellaR)
###################################
human<-readRDS(file="Human_SingCellaR_no_batch_correction.rds")
##############################################################
normalize_UMIs(human,use.scaled.factor = T)
remove_unwanted_confounders(human,residualModelFormulaStr="~UMI_count+percent_mito")

get_variable_genes_for_integrative_data_by_fitting_GLM_model(human,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(human,gmt.file = "../../../GeneSets/human.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(human)

runPCA(human,use.components=50,use.regressout.data = T)

plot_PCA_Elbowplot(human)


SingCellaR::runUMAP(human,n.dims.use = 40,n.neighbors = 30,uwot.metric = "euclidean")

plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",
                                         point.size = 0.3,mark.feature = F)
