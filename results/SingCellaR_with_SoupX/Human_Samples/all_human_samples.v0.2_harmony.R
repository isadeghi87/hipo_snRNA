#############################################
library(SingCellaR)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

plot_jaccard_similarity_among_clusters(human,cluster.type = "louvain")

plot_umap_label_by_clusters(human,show_method = "louvain",point.size = 0.3)

findMarkerGenes(human,cluster.type = "louvain")

export_marker_genes_to_table(human,cluster.type = "louvain",
                             n.TopGenes=40,min.log2FC=0.5,min.expFraction=0.3,
                             write.to.file="All_human_Top40_marker_genes.txt")

#pdf("top10_marker_genes_heatmap.pdf", width = 6, height = 8)
plot_heatmap_for_marker_genes(human,cluster.type = "louvain",
                              n.TopGenes = 8,rowFont.size = 5,use_raster = F)
#dev.off()
######################GSEA analysis##########################
human_pre_rankedGenes<-identifyGSEAPrerankedGenes_for_all_clusters(human,
                                                cluster.type = "louvain")

fgsea_Results<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = human_pre_rankedGenes, 
                                                  eps = 0,nPermSimple = 10000,
                                gmt.file = "../../../GeneSets/Kati_signature_genes_human.gmt")

plot_heatmap_for_fGSEA_all_clusters(fgsea_Results,isApplyCutoff = T,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row = 10,
                                    adjusted_pval = 0.15,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
