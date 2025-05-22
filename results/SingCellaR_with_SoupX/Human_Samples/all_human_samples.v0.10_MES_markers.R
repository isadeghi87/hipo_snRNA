#############################################
library(SingCellaR)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")



mes.cells<-subset(human@meta.data,assigned_cell_anno_layer3=="MES_like")
rest.cells<-subset(human@meta.data,assigned_cell_anno_layer3!="MES_like")

DE.genes<-identifyDifferentialGenes(objectA = human, objectB = human, cellsA = mes.cells$Cell, cellsB = rest.cells$Cell)

x.genes<-subset(DE.genes,log2FC > 1.5 & ExpFractionA > 0.3)

write.table(x.genes,file="MES-like.top_marker_genes.txt",row.names = F,col.names = T,sep="\t")
