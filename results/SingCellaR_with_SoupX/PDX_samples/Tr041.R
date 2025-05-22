library(SingCellaR)
library(SingleCellExperiment)

# Load the data set
human<-readRDS("../../../cellranger_results/PDX/K35R-Tr041/Tr041.soupX.adjusted.rds")
human.umi<-human[grep("^GRCh38",rownames(human)),]
rownames(human.umi)<-gsub("GRCh38_","",rownames(human.umi))

colnames(human.umi)<-paste(colnames(human.umi),"Tr041",sep="_")

sce<-SingleCellExperiment(assays = list(counts = human.umi))

obj<-new("SingCellaR", sce)
obj@sample_uniq_id<-"Tr041"

process_cells_annotation(obj,mito_genes_start_with ="MT-")

plot_cells_annotation(obj,type="histogram")
plot_cells_annotation(obj,type="boxplot")
dev.off()

plot_UMIs_vs_Detected_genes(obj)

DoubletDetection_with_scrublet(obj)

filter_cells_and_genes(obj,
                       min_UMIs=500,
                       max_UMIs=80000,
                       min_detected_genes=300,
                       max_detected_genes=10000,
                       max_percent_mito=20,
                       genes_with_expressing_cells = 10,
                       isRemovedDoublets = TRUE)

normalize_UMIs(obj,use.scaled.factor = T)
remove_unwanted_confounders(obj,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(obj,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.05)

remove_unwanted_genes_from_variable_gene_set(obj,gmt.file = "../../../GeneSets/human.ribosomal-mitocondrial.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
plot_variable_genes(obj)

runPCA(obj,use.components=30,use.regressout.data = T)
plot_PCA_Elbowplot(obj)

runUMAP(obj,dim_reduction_method = "pca",n.dims.use = 20,n.neighbors = 20,
        uwot.metric = "euclidean")

identifyClusters(obj,n.dims.use = 20,n.neighbors = 40,knn.metric = "euclidean")

plot_umap_label_by_clusters(obj,show_method = "louvain",point.size = 1.5)

findMarkerGenes(obj,cluster.type = "louvain")

save(obj,file="Tr041.SoupX.SingCellaR.rdata")

