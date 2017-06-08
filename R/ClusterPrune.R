#' ClusterPrune
#'
#' This function calculates p-values for each cluster and returns the pairwise distances in a cluster.
#' In the future this function will prune the clusters to remove.
#'
#' @param clustering_results The results from clustering.
#' @param Background_Module_Distances The background distribution for a module of the same size.

#' @export
#' @examples PHA_pvalues_BR<-ClusterPrune(PHA_clustering_results_P_NRED_BR,Random_Background_Module_Distances_6_BR)


# Psuedo code for how this function should work
# 1) identify the N in the module being tested and get the correct background distribution for that 'N'
# 2) identify the bins in a cluster (so each module will have multiple clusters)
# 3) get the Z-scores of all the pairwise comparisons in this cluster
# 4) t.test (one-sided) to determine if the Z-scores in the cluster is significanlty 
#     less than expected based on the background distribution
# 5) if it is significant, keep it
# 6) Bonus Round, if it is not significant, drop the lowest scoring genome from the comparison and re-test 
# (repeat until significant or <3 genomes)


ClusterPrune <- function(clustering_results,Background_Module_Distances) {
  pvalues<-rep(NA,length(clustering_results$cl))
  distances<-list(NA)
  
  for (i in 1:length(clustering_results$cl)) {
    bin_names<- colnames(clustering_results$JPE_distance)
    bin_membership <- clustering_results$cl$names[which(clustering_results$cl$membership==i)]
    bin_in_cluster<-which(bin_names %in% bin_membership)
    clusterX<-(clustering_results$JPE_distance[bin_in_cluster,bin_in_cluster])
    ttest_temp<-t.test(clusterX,Background_Module_Distances)
    
    pvalues[i]<-ttest_temp$p.value
    distances[[i]]<-clusterX
  }
  newList<-list("pvalues"=pvalues,"distances"=distances)
  return(newList)
}
