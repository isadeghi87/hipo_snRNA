library(SingCellaR)
library(SingleCellExperiment)

my.umi<-readRDS("../../../cellranger_results/human/K35R-7P6VFX/ICGC_GBM27.soupX.adjusted.rds")
colnames(my.umi)<-paste(colnames(my.umi),"ICGC_GBM27",sep="_")

sce<-SingleCellExperiment(assays = list(counts = my.umi))
human <-new("SingCellaR",sce)
human@sample_uniq_id<-"ICGC_GBM27"

process_cells_annotation(human,mito_genes_start_with ="MT-")

plot_cells_annotation(human,type="histogram")
plot_cells_annotation(human,type="boxplot")
dev.off()
plot_UMIs_vs_Detected_genes(human)

DoubletDetection_with_scrublet(human)

filter_cells_and_genes(human,
                       min_UMIs=400,
                       max_UMIs=20000,
                       min_detected_genes=300,
                       max_detected_genes=10000,
                       max_percent_mito=2,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = TRUE)

normalize_UMIs(human,use.scaled.factor = T)
remove_unwanted_confounders(human,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(human,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(human,gmt.file = "../../../GeneSets/human.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(human)

runPCA(human,use.components=30,use.regressout.data = T)
plot_PCA_Elbowplot(human)

runUMAP(human,dim_reduction_method = "pca",n.dims.use = 30,n.neighbors = 30,
        uwot.metric = "euclidean")

plot_umap_label_by_a_feature_of_interest(human,feature = "UMI_count",
                                         point.size = 0.1,mark.feature = F)

identifyClusters(human,n.dims.use = 30,n.neighbors = 40,knn.metric = "euclidean")

plot_umap_label_by_clusters(human,show_method = "louvain",point.size = 1.5)

findMarkerGenes(human,cluster.type = "louvain")


save(human,file="ICGC_GBM27.SingCellaR.rdata")
