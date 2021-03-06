# Psuedo code for how this function should work
# 1) identify the N in the module being tested and get the correct background distribution for that 'N'
# 2) identify the bins in a cluster (so each module will have multiple clusters)
# 3) get the Z-scores of all the pairwise comparisons in this cluster
# 4) t.test (one-sided) to determine if the Z-scores in the cluster is significanlty
#     less than expected based on the background distribution
# 5) if it is significant, keep it




#' ClusterPrune
#'
#' This function calculates p-values for each cluster and returns the pairwise distances in a cluster.
#' In the future this function will prune the clusters to remove.
#'
#' @param clustering_results The results from clustering.
#' @param Background_Module_Distances The background distribution for a module of the same size.
#' @param module_names A list containing all module names as index names
#' @param p_cutoff A pvalue cutoff for cluster significance (Default = 0.05)

#' @export
#' @examples PHA_pvalues_BR <- ClusterPrune(PHA_clustering_results_P_NRED_BR, Random_Background_Module_Distances_6_BR)
ClusterPrune <- function(clustering_results,
                         Background_Module_Distances,
                         module_names,
                         p_cutoff = 0.05){
  cluster_out <- clustering_results
  #Vector of the module names
  m_names <- names(module_names)
  p_values <- rep(NA, length(module_names))
  for(m_idx in 1: length(module_names)){
    tryCatch({
      #1
      n_terms    <- length(module_names[[m_idx]])
      #2
      m_clusters <- clustering_results[[m_idx]]
      bin_names  <- colnames(m_clusters$JPE_distance)

      clusters   <- levels(as.factor(m_clusters$cl$membership))

      result = tryCatch({
        for(cluster in clusters){
          cluster_Zscores <- cluster_info(cluster, m_clusters$cl$membership,
                                          m_clusters$ave_Z_score_matrix,
                                          m_clusters$cl$names)

          pval <- t.test(cluster_Zscores,
                         Background_Module_Distances[[as.character(m_idx)]],
                         alternative = "less")$p.value
          p_values[m_idx] <- pval

          if(pval > p_cutoff){
            del_positions <- which(cluster_out[[m_idx]]$cl$membership == as.numeric(cluster))

            cluster_out[[m_idx]]$cl$membership <- cluster_out[[m_idx]]$cl$membership[-del_positions]
            cluster_out[[m_idx]]$cl$names <- cluster_out[[m_idx]]$cl$names[-del_positions]
          } else{

          }

        }
      }, warning = function(w){
          return(w)
      }, error = function(e){
          # Expecting this to happen if there is no bg_distance_module calculated
          #  for the module size
      })
  }, warn = function(w)return(w)
   , err = function(e) return(w)
  )
  }
  return(list("clusters" = cluster_out, 'pvalues' = p_values))
}



cluster_info <- function(cluster, membership,ave_Z_score_matrix, names){
  scores <- ave_Z_score_matrix[names[which(membership==cluster)],
                               names[which(membership==cluster)]]
  return(scores)

}

