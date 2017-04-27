#' Build Association Matrix
#'
#' Combine all cluster results into a single matrix
#'
#' @param clustering_results The output from cluster_func
#' @param matrix_features The minimum number of times a KO term must be present to be included in the Matrix
#' @param module_names, The names of the modules (ordered)
#' @export
#' @return A matrix of clustering results. bins are rows and modules are the columns. Values represent the cluster for ecah module
#' @examples M_I_association_matrix<-fill_association_matrix(M_I_clustering_results_P_NRED,matrix_features_BR,names(M_I_module_list))


fill_association_matrix<-function(clustering_results,matrix_features,module_names) {
  association_matrix<-matrix(data = NA,
                             nrow = length(matrix_features@high_quality_bins),
                             ncol = length(clustering_results),
                             dimnames = list(matrix_features@high_quality_bins,
                                             module_names))
  
  for (i in 1:length(clustering_results)) {
    
    if (length((clustering_results[[i]]$cl$membership))!=0) {
      
      for (j in 1:length(clustering_results[[i]]$cl)) {
        for (k in 1:length(clustering_results[[i]]$cl[[j]])) {
          bin<-as.numeric(clustering_results[[i]]$cl[[j]][k])
          row<- which(matrix_features@high_quality_bins==bin)
          association_matrix[row,i]<-j
        }
      }
    }
  }
  return(association_matrix)
}