# Psuedo code for how this function should work
# 1) identify the N in the module being tested and get the correct background distribution for that 'N'
# 2) identify the bins in a cluster (so each module will have multiple clusters)
# 3) get the Z-scores of all the pairwise comparisons in this cluster
# 4) t.test (one-sided) to determine if the Z-scores in the cluster is significanlty
#     less than expected based on the background distribution
# 5) if it is significant, keep it
# 6) Bonus Round, if it is not significant, drop the lowest scoring genome from the comparison and re-test
# (repeat until significant or <3 genomes)



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

  for(m_idx in 6:6){# length(module_names)){
    #1
    n_terms    <- length(module_names[[m_idx]])
    #2
    m_clusters <- clustering_results[[m_idx]]
    bin_names  <- colnames(m_clusters$JPE_distance)

    clusters   <- levels(as.factor(m_clusters$cl$membership))

    #create pruning columns
    pruned_membername <- paste('pruned_membership_P',
                               as.character(p_cutoff), sep = "")
    pruned_names <- paste("pruned_names_P",
                          as.character(p_cutoff), sep = "")

    cluster_out[[m_idx]]$cl[[pruned_membername]] <- cluster_out[[m_idx]]$cl$membership
    cluster_out[[m_idx]]$cl[[pruned_names]] <- cluster_out[[m_idx]]$cl$names

    result = tryCatch({
      for(cluster in clusters){
        cluster_Zscores <- cluster_info(cluster, m_clusters$cl$membership,
                                        m_clusters$ave_Z_score_matrix,
                                        m_clusters$cl$names)

        pval <- t.test(cluster_Zscores,
                       transcriptional_responses$bg_distance_modules[[m_idx]],
                       alternative = "less")$p.value
        if(pval > p_cutoff){

          del_positions <- which(cluster_out[[m_idx]]$cl$pruned_membership == as.numeric(cluster))

          cluster_out[[m_idx]]$cl$pruned_membership <-
            cluster_out[[m_idx]]$cl$pruned_membership[-del_positions]
          cluster_out[[m_idx]]$cl$pruned_names <-
            cluster_out[[m_idx]]$cl$pruned_names[-del_positions]
        } else{

        }

      }
    }, warning = function(w){
      return(w)
    }, error = function(e){
      return(e)
    })
  }
  return(cluster_out)
}



cluster_info <- function(cluster, membership,
                        ave_Z_score_matrix,
                         names){
  scores <- ave_Z_score_matrix[names[which(membership==cluster)],
                            names[which(membership==cluster)]]
  return(scores)

}



# #testing
 pruned_cluster <- ClusterPrune(transcriptional_responses$clustering_results_P_NRED,
                             transcriptional_responses$bg_distance_modules,
                             list_of_all_modules
                             )

# t.test(cluster$ave_Z_score_matrix[cluster$cl$names[which(cluster$cl$membership==6)],
#      cluster$cl$names[which(cluster$cl$membership==6)]],
#      transcriptional_responses$bg_distance_modules[[6]],
#      alternative = "less")
