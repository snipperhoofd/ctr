#' Louvain clustering of Modules
#'
#' This function identifies clusters of bins with CTR based on the composite Z-scores 
#' normalized by the Jaccard distances of a cluster.
#'
#'
#' @param RNAseq_Annotated_Matrix The original matrix
#' @param Composite_Z_Score An array of all the composite Z-scores for each KO.
#' @param module_list A list of vectors.
#' @export
#' @return a list containing the following: a matrix (#bins x #modules in list) storing the cluster information,
#' a list of Z score matrices for each module (ave_Z_score_matrix),
#' a list of of matrices with distances after noralization by Jaccard Distances (JPE_distance),

#' @examples M_I_clustering_results_P_NRED <-cluster_func(RNAseq_Annotated_Matrix_BR_default_bin, all_pairwise_KO_distances$Composite_Z_Score, matrix_features_BR, M_I_module_list)
#' @examples PHA_clustering_results_P_NRED <-cluster_func(RNAseq_Annotated_Matrix_BR_default_bin, all_pairwise_KO_distances$Composite_Z_Score, matrix_features_BR, PHA_module)
# # ave_Z_score_matrix_PHA <- ave_Z_score_Func(all_pairwise_KO_distances$Composite_Z_Score[,,which(All_KOs%in%PHA_module)])
# J_PHA <-Jaccard_Distance_Function(RNAseq_Annotated_Matrix_BR_default_bin,matrix_features_BR, PHA_module)

cluster_func<-function(RNAseq_Annotated_Matrix, Composite_Z_Score, matrix_features, module_list) {

  if (is.list(module_list)) {
    newList <- list()
    All_KOs<-unique(unlist(module_list, use.names=FALSE))
    
    for (i in 1:length(module_list)) {
    Jaccard_Distance <- Jaccard_Distance_Function(RNAseq_Annotated_Matrix, 
                                                  matrix_features, 
                                                  module_list[[i]])
    if (sum(All_KOs%in%CCM_module_list[[i]])==1) {next} else { 
    ave_Z_score_matrix <- ave_Z_score_Func(Composite_Z_Score[, , which(All_KOs%in%module_list[[i]])])
    JPE_distance<-ave_Z_score_matrix*(1-Jaccard_Distance)
    rownames(JPE_distance)<-colnames(JPE_distance)
    JPE_distance_Table <- subset(melt(JPE_distance), value!=0)
    edge_frame<-as.data.frame(JPE_distance_Table)
    colnames(edge_frame)<-c("genome_1", "genome_2", "Score")
    edge_frame_reduced<-edge_frame[which(edge_frame[, 3] < 0), ] # remove negative numbers as these cannot be used in clustering algorithm
    edge_frame_reduced[,3]<-(-edge_frame_reduced[, 3])
    graph<-graph_from_data_frame(edge_frame_reduced, directed=FALSE)
    cl<-cluster_louvain(graph, weights=E(graph)$Score)
    
    minilist <- list("ave_Z_score_matrix" = ave_Z_score_matrix,
                       "JPE_distance"= JPE_distance,
                       "JPE_distance_Table" = JPE_distance_Table,
                       "cl" = cl)
    newList[[i]] <- minilist   
    }
    }
    
    } else {
      Jaccard_Distance <- Jaccard_Distance_Function(RNAseq_Annotated_Matrix, 
                                                    matrix_features, 
                                                    module_list)
      ave_Z_score_matrix <- ave_Z_score_Func(all_pairwise_KO_distances$Composite_Z_Score[,,which(matrix_features@All_KOs%in%module_list)])
      JPE_distance<-ave_Z_score_matrix*(1-Jaccard_Distance)
      rownames(JPE_distance)<-colnames(JPE_distance)
      JPE_distance_Table <- subset(melt(JPE_distance), value!=0)
      edge_frame<-as.data.frame(JPE_distance_Table)
      colnames(edge_frame)<-c("genome_1","genome_2","Score")
      edge_frame_reduced<-edge_frame[which(edge_frame[,3]<0),] # remove negative numbers as these cannot be used in clustering algorithm
      edge_frame_reduced[,3]<-(-edge_frame_reduced[,3])
      graph<-graph_from_data_frame(edge_frame_reduced,directed=FALSE)
      cl<-cluster_louvain(graph,weights=E(graph)$Score)

      newList <- list("ave_Z_score_matrix" = ave_Z_score_matrix,
                      "JPE_distance" = JPE_distance,
                      "JPE_distance_Table" = JPE_distance_Table,
                      "cl" = cl)
      }
  
  return(newList)  
  }
