#############################################
library(SingCellaR)
library(SingleCellExperiment)
library(sceasy)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

#############################################
my.umi<-get_umi_count(human)
#############################################
cell.info<-human@meta.data
cell.info<-cell.info[,-c(3,4,5,6,7)]
rownames(cell.info)<-cell.info$Cell
cell.info<-cell.info[,-c(1)]

s.umap<-get_umap.result(human)
rownames(s.umap)<-s.umap$Cell
s.umap<-s.umap[rownames(cell.info),]
s.umap<-s.umap[,c("UMAP1","UMAP2")]

sce <- SingleCellExperiment(assays=list(counts=my.umi),colData=cell.info)

reducedDim(sce, "UMAP") <- s.umap

sceasy::convertFormat(sce, from="sce", to="anndata",
                      outFile='all_human_samples.h5ad')
