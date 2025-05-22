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
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[2]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "OPC/NPC-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl2"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[3]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "MES-like",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl3","cl4","cl5"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[4]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Intermediate_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl6"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[5]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl7"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[6]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Mitotic_cancer_cell_state",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl7"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[7]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Microglia",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl8"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[8]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Neuron",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl9"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[9]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Plasma",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl10"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[10]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Oligodendrocyte",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl11"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[11]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Endothelial",AUCell_cutoff = 0.1,
                                        AUCell_score = AUC.scores,selected.limited.clusters = c("cl12"),
                                        IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                        ,showLegend = F)

p[[12]]<-plot_umap_label_by_AUCell_score(human,AUCell_gene_set_name = "Pericyte",AUCell_cutoff = 0.1,
                                         AUCell_score = AUC.scores,selected.limited.clusters = c("cl13"),
                                         IsLimitedAUCscoreByClusters = T,point.size = 0.1
                                         ,showLegend = F)

pdf(file="AUCell-Human-Primary.pdf", width=20, height=10, onefile=T, bg="transparent",fonts = NULL,useDingbats=FALSE)
grid.arrange(grobs = p, ncol = 4, nrow=3, as.table = FALSE)
dev.off()
