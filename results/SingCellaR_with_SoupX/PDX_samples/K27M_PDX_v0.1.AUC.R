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
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

pdf(file="AUCell-PDX-K27M-Cancer.pdf", width=18, height=2, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 5, nrow=1, as.table = FALSE)
dev.off()


p <- list()

p[[1]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Astrocyte-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.2,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,point.size = 0.1,selected.sampleID=c("1_Astro110fh","1_SJ_DIPG_X7_XP1","1_SU_DIPG_XIII","1_Tr041"),
                                        IsShowOnlySampleIDs = T, IsDownsample = T,downsample.size = 13262
                                        ,showLegend = F)

pdf(file="AUCell-PDX-K27M-Cancer.pdf", width=18, height=2, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 5, nrow=1, as.table = FALSE)
dev.off()

