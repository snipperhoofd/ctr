#' Make figures from arrays
#'
#' This function creates standardized figures to visualize differences between clusters
#' 
#' @param array An array of all the pairwise Z-scores for a module
#' @param module the name of the KOs
#' @param xlimits a vector with the limits of the X axis
#' @param ylimits a vector with the limits of the Y axis
#' @param background The background distribution for individual KOs
#' @param cluster_results the results from the clustering algorithm

#' @return a multi-panel figure. The first column is the distributions, each subsequent column is the cluster. Rows represent each individual KO
#' 
#' @examples for(i in 1:length(PHA_clustering_results_P_NRED$cl)) {
#' array_fig(PHA_module_P_NRED$combined,cluster_text_matrix[,i],c(-4,4),c(0,1),I_KOs_Background$random_pairwise_gene_correlation,PHA_clustering_results_P_NRED$cl)
#' }

array_fig<-function(array,module,xlimits,ylimits,background,cluster_results) {
  colors<-rainbow(dim(array)[3])
  plot(density(background,na.rm=TRUE),ylim=ylimits,xlim=xlimits,lwd=2,main=paste("Cluster",i))
  legend("topright",c("background",module),col=c("black",colors),lwd=2)
  for (j in 1:dim(array)[3])
    if(sum(!is.na(array[which(rownames(array)%in%cluster_results$names[which(cluster_results$membership==i)]),which(colnames(array)%in%cluster_results$names[which(cluster_results$membership==i)]),j]))<=2) {next} else {
      points(density(array[which(rownames(array)%in%cluster_results$names[which(cluster_results$membership==i)]),which(colnames(array)%in%cluster_results$names[which(cluster_results$membership==i)]),j],na.rm=TRUE),type="l",col=colors[j],lwd=2)
    }
}


#USAGE: array_fig(glutamate_aspartate_transport_module_P_NRED$combined,cluster_text_matrix[,i],c(-4,4),c(0,1.5),I_KOs_Background$random_pairwise_gene_correlation,glutamate_aspartate_transport_clustering_results_P_NRED$cl)