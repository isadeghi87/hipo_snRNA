###################################
library(SingCellaR)
library(AUCell)
library(gridExtra)

human<-readRDS(file="K27M_PDX_samples.SoupX.SingCellaR.Harmony_v1.rds")
Build_AUCell_Rankings(human,AUCell_buildRankings.file="K27M_PDX_rankings.AUCells.rdata")

set.seed(1)
AUC.scores<-Run_AUCell(human,"K27M_PDX_rankings.AUCells.rdata","../../../GeneSets/human.K27M_MET.gmt")

p <- list()


p[[1]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Astrocyte-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,point.size = 0.1,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T,IsDownsample = F,downsample.size = 19646,
                                        ,showLegend = T)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.2,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.3,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.2,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[6]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Microglia",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[7]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Neuron",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[8]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Plasma",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[9]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Oligodendrocyte",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[10]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Endothelial",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

p[[11]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Pericyte",AUCell_cutoff = 0.2,IsLimitedAUCscoreByClusters = T,
                                        selected.limited.clusters = c(paste('cl',1:15,sep="")),
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T
                                        ,showLegend = T)

pdf(file="AUCell-PDX-K27M-Cancer_Normal_cell_number.pdf", width=18, height=8, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 4, nrow=3, as.table = FALSE)
dev.off()