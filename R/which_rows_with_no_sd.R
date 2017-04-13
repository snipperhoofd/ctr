#' Removes rows with no standard deviation
#'
#' This function removes all rows with a standard deviation of 0.
#' This is necessary because Pearson Correlations cannot be calcualted in these cases.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @param matrix_features The General_features object containing matrix feature information
#' @export
#' @return The matrix with rows that have a standard deviation of 0 removed
#' @examples RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)
which_rows_with_no_sd<- function(RNAseq_Annotated_Matrix, matrix_features){

  sample_cols <- matrix_features@SS:matrix_features@SE
  data_matrix <- matrix(apply(RNAseq_Annotated_Matrix[, sample_cols], 2, as.numeric), ncol = length(sample_cols))
  print(data_matrix[1:5,])
  stdev_vector <- which_rows_with_no_sd_cpp(data_matrix)
  RNAseq_Annotated_Matrix <- cbind(RNAseq_Annotated_Matrix, stdev_vector)

  return(RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix[,ncol(RNAseq_Annotated_Matrix)] != 0),])


}
