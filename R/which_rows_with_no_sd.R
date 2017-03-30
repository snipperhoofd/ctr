#' Removes rows with no standard deviation
#'
#' This function removes all rows with a standard deviation of 0.
#' This is necessary because Pearon Correlations cannot be calcualted in these cases.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#'
#' @return The matrix with rows that have a standard deviation of 0 removed
#' @examples RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)

which_rows_with_no_sd <- function(RNAseq_Annotated_Matrix){
  sd_of_zero<-NULL
  for (i in 1:dim(RNAseq_Annotated_Matrix)[1]) {if(sd(RNAseq_Annotated_Matrix[i,SS:SE])==0) {sd_of_zero<-c(sd_of_zero,i)}}
  RNAseq_Annotation_Matrix_no_sd_of_zero<-RNAseq_Annotated_Matrix[-sd_of_zero,]
  return(RNAseq_Annotation_Matrix_no_sd_of_zero)
}
