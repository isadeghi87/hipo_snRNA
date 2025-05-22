#############################################
library(SingCellaR)
library(SingleCellExperiment)
library(sceasy)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

#############################################
merge_clusters(human,cluster.type = "louvain",merge_cluster_ids = c('cl3:cl4:cl5'))

findMarkerGenes(human,cluster.type = "merged_louvain")

plot_umap_label_by_clusters(human,show_method = "merged_louvain")

human@meta.data$merged_louvain_cluster<-human@sc.clusters$merged_louvain

saveRDS(human,file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation_merged_clusters.rds")
#############################################
export_marker_genes_to_table(human,cluster.type = "merged_louvain",
                             n.TopGenes=10,min.log2FC=0.5,min.expFraction=0.3,
                             write.to.file="All_human_Top10_marker_genes_merged_louvain.txt")
