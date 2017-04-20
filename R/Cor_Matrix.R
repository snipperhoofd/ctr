#' Calculates a correlation matrix between two vectors of rows
#'
#' This function calculates a pairwise correlation matrix for both pearson and NRED.
#' Returns an array containing two matrices.
#'
#' @param vector1,vector2 The two lists rows in which a particular KO is identified
#' @param RNAseq_Annotation_Matrix_no_sd_of_zero The annotated matrix of RNAseq data with rows that
#' have a standard deviation of zero removed.
#' @export
#' @return a pairwise correlation matrix,
#' @examples Cor_Matrix(KO_position_of_A,KO_position_of_B,RNAseq_Annotation_Matrix_no_sd_of_zero)

Cor_Matrix <- function(vector1, vector2, input_matrix, matrix_features){
  pairwise_gene_correlation<-matrix(NA,nrow=length(vector1),ncol=length(vector2))
  pairwise_gene_euclidean<-pairwise_gene_correlation
  for (m in 1:length(vector1)){
    for (n in 1:length(vector2)){
      # Pearson correlation
      pairwise_gene_correlation[m,n]<-(cor(as.numeric(input_matrix[vector1[m],matrix_features@SS:matrix_features@SE]),
                                             as.numeric(input_matrix[vector2[n],matrix_features@SS:matrix_features@SE])))
      # Euclidean distance
      subtracted_lists<- input_matrix[vector1[m],matrix_features@RS:matrix_features@RE] - input_matrix[vector2[n],matrix_features@RS:matrix_features@RE]
      pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))

    }
  }
  return(array(c(pairwise_gene_correlation,pairwise_gene_euclidean),dim=c(length(vector1),length(vector2),2)))
}
