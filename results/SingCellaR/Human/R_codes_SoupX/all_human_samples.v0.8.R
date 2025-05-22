#############################################
library(SingCellaR)

human<-readRDS(file="Human_all_samples_SingCellaR.Harmony.v2.with.final_annotation.rds")

######Qeustion1##########
table(human@meta.data$assigned_cell_anno_layer3)

#####Question2###########
microglia<-subset(human@meta.data,assigned_cell_anno_layer3=="Microglia")
non_microglia<-subset(human@meta.data,assigned_cell_anno_layer3 !="Microglia")

DE.genes<-identifyDifferentialGenes(objectA = human, objectB = human, 
                                    cellsA = microglia$Cell, cellsB = non_microglia$Cell)
head(DE.genes)

####Question3############
my_umi<-get_umi_count(human)
#########################

