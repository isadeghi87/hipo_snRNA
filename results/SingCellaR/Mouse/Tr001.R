library(SingCellaR)

data_matrices_dir<-"../../../cellranger_results/mouse/K35R-3Q7AGR/filtered_feature_bc_matrix/"
mouse<-new("SingCellaR")
mouse@dir_path_10x_matrix<-data_matrices_dir
mouse@sample_uniq_id<-"Tr001"

load_matrices_from_cellranger(mouse,cellranger.version = 7)
process_cells_annotation(mouse,mito_genes_start_with ="mt-")

plot_cells_annotation(mouse,type="histogram")
plot_cells_annotation(mouse,type="boxplot")
dev.off()
plot_UMIs_vs_Detected_genes(mouse)

DoubletDetection_with_scrublet(mouse)

filter_cells_and_genes(mouse,
                       min_UMIs=1000,
                       max_UMIs=100000,
                       min_detected_genes=500,
                       max_detected_genes=15000,
                       max_percent_mito=5,
                       genes_with_expressing_cells = 5,
                       isRemovedDoublets = TRUE)

normalize_UMIs(mouse,use.scaled.factor = T)
remove_unwanted_confounders(mouse,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(mouse,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(mouse,gmt.file = "../../../GeneSets/mouse.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(mouse)

runPCA(mouse,use.components=30,use.regressout.data = T)
plot_PCA_Elbowplot(mouse)

runUMAP(mouse,dim_reduction_method = "pca",n.dims.use = 20,n.neighbors = 30,
        uwot.metric = "euclidean")

identifyClusters(mouse,n.dims.use = 20,n.neighbors = 30,knn.metric = "euclidean")

plot_umap_label_by_clusters(mouse,show_method = "louvain",point.size = 1.5)

findMarkerGenes(mouse,cluster.type = "louvain")

save(mouse,file="Tr001.SingCellaR.rdata")
