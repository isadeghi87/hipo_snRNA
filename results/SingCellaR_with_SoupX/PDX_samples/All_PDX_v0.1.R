###################################
library(SingCellaR)
library(harmony)

human<-readRDS(file="PDX_SoupX.SingCellaR.No_batch_correction.rds")

plot_umap_label_by_qc(human)

normalize_UMIs(human,use.scaled.factor = T)

remove_unwanted_confounders(human,
                            residualModelFormulaStr="~UMI_count+percent_mito")

runPCA(human,use.components=50,use.regressout.data = T)

plot_PCA_Elbowplot(human)
#########################################################
SrunHarmony <- function(object,n.dims.use=30,covariates=c("data_set"),harmony.theta = NULL,harmony.lambda = NULL,harmony.sigma = 0.1,harmony.nclust = NULL,
                        harmony.max.iter = 10,n.seed=1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  ######################################################
  orig.pca<-get_pca.result(object)$x[,1:n.dims.use]
  cell.meta<-get_cells_annotation(object)
  cell.meta<-subset(cell.meta,IsPassed==TRUE)
  ######################################################
  ###Run Harmony########################################
  set.seed(seed=n.seed)
  ######################################################
  harmony_embeddings <- harmony::RunHarmony(orig.pca, 
                                            cell.meta,  
                                            covariates,
                                            theta = harmony.theta,
                                            lambda = harmony.lambda,
                                            sigma = harmony.sigma,
                                            nclust = harmony.nclust,
                                            max_iter = harmony.max.iter,
                                            verbose=TRUE)
  ######################################################
  object@Harmony.embeddings<-harmony_embeddings
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Harmony analysis is done!.")
}

################Harmony2##################
SrunHarmony(human,covariates = c("sampleID","sample_type","tenx_version","sequencer"),
            n.dims.use = 20,harmony.max.iter = 10,
            n.seed = 1,harmony.theta = c(3,2,2,2))

SingCellaR::runUMAP(human,useIntegrativeEmbeddings = T, 
                    integrative_method = "harmony",n.dims.use = 20,
                    n.neighbors = 30,uwot.metric = "euclidean")

plot_umap_label_by_qc(human)
plot_umap_label_by_a_feature_of_interest(human,feature = "sampleID",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "sample_name",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "sample_type",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "sequencer",point.size = 1)
plot_umap_label_by_a_feature_of_interest(human,feature = "tenx_version",point.size = 1)


SingCellaR::identifyClusters(human,useIntegrativeEmbeddings = T,
                             integrative_method = "harmony", 
                             n.dims.use = 10,n.neighbors = 40, 
                             knn.metric = "euclidean")

plot_umap_label_by_clusters(human,show_method = "louvain")

saveRDS(human,file="all_PDX_samples.SoupX.SingCellaR.Harmony_v1.rds")
################CCA######################