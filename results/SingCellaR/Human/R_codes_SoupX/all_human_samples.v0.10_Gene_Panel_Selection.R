#############################################
library(SingCellaR)
library(AUCell)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")


my_markers<-read.table("All_human_Top30_marker_genes.txt",header=T)

cl6<-subset(my_markers,cluster=="cl6")

plot_umap_label_by_genes(human,gene_list = cl6$Gene,point.size = 0.1)
