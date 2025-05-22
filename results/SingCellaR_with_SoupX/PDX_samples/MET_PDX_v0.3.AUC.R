###################################
library(SingCellaR)
library(AUCell)
library(gridExtra)

human<-readRDS(file="MET_PDX_samples.SoupX.SingCellaR.Harmony_v1.rds")
Build_AUCell_Rankings(human,AUCell_buildRankings.file="MET_PDX_rankings.AUCells.rdata")

set.seed(1)
AUC.scores<-Run_AUCell(human,"MET_PDX_rankings.AUCells.rdata","../../../GeneSets/human.K27M_MET.gmt")

p <- list()

p[[1]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Astrocyte-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T,showLegend = T)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.3,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[6]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Microglia",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[7]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Neuron",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[8]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Plasma",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[9]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Oligodendrocyte",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[10]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Endothelial",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[11]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Pericyte",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:13,sep="")),
                                        selected.sampleID=c("1_RCXXB47_XP1_pellet","1_Tr068"),
                                        IsDownsample = T,downsample.size = 13262,
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)


pdf(file="AUCell-PDX-MET-P1_Cancer_Nornmal.pdf", width=18, height=8, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 4, nrow=3, as.table = FALSE)
dev.off()
