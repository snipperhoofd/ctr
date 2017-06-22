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
fill_association_matrix <- function(clustering_results, matrix_features, module_names){

  association_matrix <- matrix(data = NA,
                               nrow = length(matrix_features@high_quality_bins),
                               ncol = length(clustering_results$clusters),
                               dimnames = list(matrix_features@high_quality_bins,
                                               names(module_names)))


  for(m_idx in 1: length(module_names)) {
    cluster      <- clustering_results$clusters[[m_idx]]$cl
    present_bins <- cluster$names
    memberships  <- cluster$membership
    bin_indices  <- which(matrix_features@high_quality_bins %in% present_bins)

    for ( i in 1: length(bin_indices) ) {
      bin_idx <- bin_indices[i]
      member  <- memberships[i]

      if (length(member) != 0) {
        if (! is.na(member) ) {
          association_matrix[bin_idx, m_idx] <- member
        }
      }
    }
  }
  return(association_matrix)
}
