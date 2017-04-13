#' Removes rows with no standard deviation
#'
#' This function removes all rows with a standard deviation of 0.
#' This is necessary because Pearson Correlations cannot be calcualted in these cases.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @export
#' @return The matrix with rows that have a standard deviation of 0 removed
#' @examples RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)

which_rows_with_no_sd <- function(RNAseq_Annotated_Matrix, matrix_features){
  sd_of_zero <- c()

  for (i in 1 : nrow(RNAseq_Annotated_Matrix)) {
    if(sd(RNAseq_Annotated_Matrix[i, matrix_features@SS : matrix_features@SE]) == 0) {
      sd_of_zero<-c(sd_of_zero, i)
    }
  }
  RNAseq_Annotation_Matrix_no_sd_of_zero <- RNAseq_Annotated_Matrix[-sd_of_zero,]
  return(RNAseq_Annotation_Matrix_no_sd_of_zero)
}

which_rows_with_no_sd2 <- function(RNAseq_Annotated_Matrix, matrix_features){
  #sourceCpp("src/which_rows_with_no_sd.cpp")
  #filtered_matrix still contains empty row
  sample_cols <- matrix_features@SS:matrix_features@SE
  num_matrix <- matrix(apply(RNAseq_Annotated_Matrix[, sample_cols], 2, as.numeric), ncol = length(sample_cols))
  filtered_matrix <- which_rows_with_no_sd_cpp(num_matrix)
 # print('ja')
 # return(filtered_matrix)
#  return(filtered_matrix[which(filtered_matrix[, ncol(filtered_matrix)] != "0"), 1:ncol(filtered_matrix)])

}
