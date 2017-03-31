#' Louvain clusters
#'
#' This function uses the average composite Z-score normalized by the Jaccard distances to cluster
#' genome bins to identify clusters with congruent transcriptional responses.
#'
#' @param Zscore_Matrix An array of all the pairwise Z-scores for a module
#' @param Jaccard_Distance An array of all the pairwise Jaccard Distances for a module
#' @export
#' @return a list containing the average Z score matrix (ave_Z_score_matrix),
#' the normalized distances after being noralized by Jaccard Distances (JPE_distance),
#' the JPE distances formatted for a network representation (e.g. cytoscpe),
#' and the cluster assignments (cl).

#' @examples PHA_clustering_results_P_NRED <-cluster_func(PHA_module_P_NRED$combined,Jaccard_Distance_PHA)

cluster_func<-function(Zscore_Matrix,Jaccard_Distance) {
  ave_Z_score_matrix <-ave_Z_score_Func(Zscore_Matrix)
  JPE_distance<-ave_Z_score_matrix*(1-Jaccard_Distance)
  rownames(JPE_distance)<-colnames(JPE_distance)
  JPE_distance_Table <- subset(melt(JPE_distance), value!=0)
  edge_frame<-as.data.frame(JPE_distance_Table)
  colnames(edge_frame)<-c("genome_1","genome_2","Score")
  edge_frame_reduced<-edge_frame[which(edge_frame[,3]<0),] # remove negative numbers as these cannot be used in clustering algorithm
  edge_frame_reduced[,3]<-(-edge_frame_reduced[,3])
  graph<-graph_from_data_frame(edge_frame_reduced,directed=FALSE)
  cl<-cluster_louvain(graph,weights=E(graph)$Score)

  newList <- list("ave_Z_score_matrix" = ave_Z_score_matrix,"JPE_distance"=JPE_distance,"JPE_distance_Table"=JPE_distance_Table,"cl"=cl)

  return(newList)

}
