#############################################
library(SingCellaR)
library(SingleCellExperiment)
library(sceasy)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

#export_marker_genes_to_table(human,cluster.type = "louvain",
#                             n.TopGenes=5,min.log2FC=0.5,min.expFraction=0.3,
#                             write.to.file="All_human_Top10_marker_genes.txt")
#
#my_markers<-read.table("All_human_Top10_marker_genes.txt",header=T)


cl1.genes<-c("GFAP","TNC")
cl2.genes<-c("OPCML")
mes.genes<-c("VEGFA","AKAP12","PRKCB")
cl6.genes<-c("MT3","CLU","S100A6","DBI")
cl7.genes<-c("RRM2","AURKB")
cl8.genes<-c("CCL3","CSF2RA","SRGN")
cl9.genes<-c("RBFOX3","MYT1L")
cl10.genes<-c("IGHG1","MZB1")
cl11.genes<-c("MBP","PLP1")
cl12.genes<-c("CLDN5","VWF")
cl13.genes<-c("COL3A1","COL1A1")

my.genes<-c(cl1.genes,cl2.genes,mes.genes,
            cl6.genes,cl7.genes,cl8.genes,
            cl9.genes,cl10.genes,cl11.genes,
            cl12.genes,cl13.genes)

plot_bubble_for_genes_per_cluster(human,cluster.type = "louvain",IsApplyClustering = T,
                                  IsClusterByRow = T,IsClusterByColumn = T,
                                  gene_list=my.genes,show.percent = T,buble.scale = 8,point.color1 = "white",
                                  point.color2 = "red")