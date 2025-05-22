#############################################
library(SingCellaR)
library(SingleCellExperiment)
library(sceasy)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

K27M<-subset(human@meta.data,sample_type=="K27M")
MET<-subset(human@meta.data,sample_type=="MET")

k27m<-as.data.frame(table(K27M$assigned_cell_anno_layer3))
colnames(k27m)<-c("cell_annotation","K27M")
met<-as.data.frame(table(MET$assigned_cell_anno_layer3))
colnames(met)<-c("cell_annotation","MET")

my.freq<-merge(k27m,met,all.x=T,all.y=T)

#write.csv(my.freq,file="cell_fraction.csv",quote = F,row.names = F)

K27_p1_GBM27<-subset(human@meta.data,sample_name=="K27M_p1_GBM27")
K27M_p2_GBM60<-subset(human@meta.data,sample_name=="K27M_p2_GBM60")
K27M_p3_GBM96<-subset(human@meta.data,sample_name=="K27M_p3_GBM96")
MET_p2_GBM15<-subset(human@meta.data,sample_name=="MET_p2_GBM15")
MET_p2_GBM71<-subset(human@meta.data,sample_name=="MET_p2_GBM71")

gbm27<-as.data.frame(table(K27_p1_GBM27$assigned_cell_anno_layer3))
colnames(gbm27)<-c("cell_annotation","K27_p1_GBM27")

gbm60<-as.data.frame(table(K27M_p2_GBM60$assigned_cell_anno_layer3))
colnames(gbm60)<-c("cell_annotation","K27M_p2_GBM60")

gbm96<-as.data.frame(table(K27M_p3_GBM96$assigned_cell_anno_layer3))
colnames(gbm96)<-c("cell_annotation","K27M_p3_GBM96")

gbm15<-as.data.frame(table(MET_p2_GBM15$assigned_cell_anno_layer3))
colnames(gbm15)<-c("cell_annotation","MET_p2_GBM15")

gbm71<-as.data.frame(table(MET_p2_GBM71$assigned_cell_anno_layer3))
colnames(gbm71)<-c("cell_annotation","MET_p2_GBM71")

freq2<-merge(gbm27,gbm60,all.x=T,all.y=T)
freq2.1<-merge(freq2,gbm96,all.x=T,all.y=T)
freq2.2<-merge(freq2.1,gbm15,all.x=T,all.y=T)
freq2.3<-merge(freq2.2,gbm71,all.x=T,all.y=T)

write.csv(freq2.3,file="cell_fraction_per_sample.csv",quote = F,row.names = F)


