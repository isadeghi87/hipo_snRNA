#############################################
library(SingCellaR)
library(AUCell)
library(gridExtra)
library(ggplot2)
library(grid)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")
Build_AUCell_Rankings(human,AUCell_buildRankings.file="human_K27M_MET_cells_rankings.AUCells.rdata")

set.seed(1)
AUC.scores<-Run_AUCell(human,"human_K27M_MET_cells_rankings.AUCells.rdata",
                                      "../../Human/../../GeneSets/human.K27M_MET.gmt")

p <- list()

p[[1]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Astrocyte-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl1"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 3203,point.size = 0.1
                                        ,showLegend = T)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl2"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 3203,point.size = 0.1
                                        ,showLegend = T)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl3","cl4","cl5"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 3203,point.size = 0.1
                                        ,showLegend = T)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl6"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 3203,point.size = 0.1
                                        ,showLegend = T)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl7"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 3203,point.size = 0.1
                                        ,showLegend = T)

pdf(file="AUCell-Human-Cancel-K27M_Cell_number.pdf", width=30, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 5, nrow=1, as.table = FALSE)
dev.off()

p <- list()

p[[1]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Astrocyte-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl1"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 5428,point.size = 0.1
                                        ,showLegend = T)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl2"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 5428,point.size = 0.1
                                        ,showLegend = T)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl3","cl4","cl5"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 5428,point.size = 0.1
                                        ,showLegend = T)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl6"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 5428,point.size = 0.1
                                        ,showLegend = T)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl7"),
                                        IsLimitedAUCscoreByClusters = T,selected.sampleID = c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                        IsShowOnlySampleIDs = T,IsDownsample = T,downsample.size = 5428,point.size = 0.1
                                        ,showLegend = T)

pdf(file="AUCell-Human-Cancel-MET_cell_number.pdf", width=30, height=6, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 5, nrow=1, as.table = FALSE)
dev.off()
