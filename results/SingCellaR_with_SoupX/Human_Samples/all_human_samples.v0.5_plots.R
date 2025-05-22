#############################################
library(SingCellaR)
library(SingleCellExperiment)
library(sceasy)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

K27M<-subset(human@meta.data,sample_type=="K27M")
MET<-subset(human@meta.data,sample_type=="MET")

plot_umap_label_by_selected_sampleID(human,c("1_ICGC_GBM27","1_ICGC_GBM60","1_ICGC_GBM96"),
                                         point.front.color = "red",point.size = 1)

plot_umap_label_by_selected_sampleID(human,c("1_ICGC_GBM15","1_ICGC_GBM71"),
                                     point.front.color = "blue",point.size = 1)

plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",
                                         point.size = 1,mark.feature = F)



