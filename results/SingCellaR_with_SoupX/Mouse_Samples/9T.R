library(SingCellaR)
library(SingleCellExperiment)

my.umi<-readRDS("../../../cellranger_results/mouse/K35R-D2BCNW/9T.SoupX.adjusted.rds")
colnames(my.umi)<-paste(colnames(my.umi),"9T",sep="_")

sce<-SingleCellExperiment(assays = list(counts = my.umi))
mouse <-new("SingCellaR",sce)
mouse@sample_uniq_id<-"9T"

process_cells_annotation(mouse,mito_genes_start_with ="mt-")

plot_cells_annotation(mouse,type="histogram")
plot_cells_annotation(mouse,type="boxplot")
dev.off()
plot_UMIs_vs_Detected_genes(mouse)

DoubletDetection_with_scrublet(mouse)

filter_cells_and_genes(mouse,
                       min_UMIs=1000,
                       max_UMIs=50000,
                       min_detected_genes=800,
                       max_detected_genes=10000,
                       max_percent_mito=5,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = TRUE)

normalize_UMIs(mouse,use.scaled.factor = T)
remove_unwanted_confounders(mouse,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(mouse,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(mouse,gmt.file = "../../../GeneSets/mouse.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(mouse)

runPCA(mouse,use.components=30,use.regressout.data = T)
plot_PCA_Elbowplot(mouse)

runUMAP(mouse,dim_reduction_method = "pca",n.dims.use = 30,n.neighbors = 30,
        uwot.metric = "euclidean")

identifyClusters(mouse,n.dims.use = 30,n.neighbors = 30,knn.metric = "euclidean")

plot_umap_label_by_clusters(mouse,show_method = "louvain",point.size = 1.5)

#findMarkerGenes(mouse,cluster.type = "louvain")

save(mouse,file="9T.SoupX.SingCellaR.rdata")
