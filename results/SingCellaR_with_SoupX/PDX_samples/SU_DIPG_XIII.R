library(SingCellaR)
library(SingleCellExperiment)

human<-readRDS("../../../cellranger_results/PDX/K35R-sn_SU_DIPG_13/SU_DIPG_13.soupX.adjusted.rds")
human.umi<-human[grep("^GRCh38",rownames(human)),]
rownames(human.umi)<-gsub("GRCh38_","",rownames(human.umi))

colnames(human.umi)<-paste(colnames(human.umi),"SU_DIPG_XIII",sep="_")

sce<-SingleCellExperiment(assays = list(counts = human.umi))

obj<-new("SingCellaR", sce)
obj@sample_uniq_id<-"SU_DIPG_XIII"

process_cells_annotation(obj,mito_genes_start_with ="MT-")

plot_cells_annotation(obj,type="histogram")
plot_cells_annotation(obj,type="boxplot")
dev.off()

plot_UMIs_vs_Detected_genes(obj)

DoubletDetection_with_scrublet(obj)

filter_cells_and_genes(obj,
                       min_UMIs=1000,
                       max_UMIs=150000,
                       min_detected_genes=1000,
                       max_detected_genes=12000,
                       max_percent_mito=5,
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

identifyClusters(obj,n.dims.use = 20,n.neighbors = 30,knn.metric = "euclidean")

plot_umap_label_by_clusters(obj,show_method = "louvain",point.size = 1.5)

#findMarkerGenes(obj,cluster.type = "louvain")

save(obj,file="SU_DIPG_XIII.SoupX.SingCellaR.rdata")

