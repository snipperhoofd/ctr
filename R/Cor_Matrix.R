#' Calculates a correlation matrix between two vectors of rows
#'
#' This function calculates a pairwise correlation matrix.
#'
#' @param vector1,vevctor2 The two lists rows in which a particular KO is identified
#' @param RNAseq_Annotation_Matrix_no_sd_of_zero The annotated matrix of RNAseq data with rows that
#' have a standard deviation of zero removed.
#' @export
#' @return a pairwise correlation matrix,
#' @examples Cor_Matrix(KO_position_of_A,KO_position_of_B,RNAseq_Annotation_Matrix_no_sd_of_zero)

Cor_Matrix <- function(vector1, vector2, input_matrix){
  pairwise_gene_correlation<-matrix(NA,nrow=length(vector1),ncol=length(vector2))
  pairwise_gene_euclidean<-pairwise_gene_correlation
  for (m in 1:length(vector1)){
    for (n in 1:length(vector2)){
      # cannot calculate Corrleation with no standard error. If there is a sd, proceed with calculations
      if (sd(as.numeric(input_matrix[vector1[m],2:7]))!=0 &
          sd(as.numeric(input_matrix[vector2[n],2:7]))!=0) {
        pairwise_gene_correlation[m,n]<-(cor(as.numeric(input_matrix[vector1[m],2:7]),
                                             as.numeric(input_matrix[vector2[n],2:7])))
        subtracted_lists<- input_matrix[vector1[m],10:15]-
          input_matrix[vector2[n],10:15]
        pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
      } else {
        # If there is no standard deviation, the correlation is NA
        pairwise_gene_correlation[m,n]<-NA
        subtracted_lists<- input_matrix[position_of_kegg_enzyme_A[m],10:15]-
          input_matrix[position_of_kegg_enzyme_B[n],10:15]
        pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
      }
    }
  }
  return(array(c(pairwise_gene_correlation,pairwise_gene_euclidean),dim=c(length(vector1),length(vector2),2)))
}
