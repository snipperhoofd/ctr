#' Build Presence Absence Matrix of Annotations
#'
#' Using the annotation column(s), this function will build a presence absence matrix
#' that will be used for calculating the Jaccard Distance. The matrix will be X*Y, where
#' X is the number of bins, and Y is the number of annotated features.
#'
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @param N, The minimum number of times a KO term must be present to be included in the Matrix
#' @export
#' @return An X*Y matrix where X = # bins, and Y = # of annotated features.
#' @examples Pairwise_Bin_Array_Presence	<- Presence_Absence_Matrix(RNAseq_Annotated_Matrix,5)



Presence_Absence_Matrix <- function(RNAseq_Annotated_Matrix, N) {
  if (missing(N)) {N = 1}
  dim_matrix <-length(table(RNAseq_Annotated_Matrix$Bin))
  no_annotation <- which(names((table(RNAseq_Annotated_Matrix$KO))) == "")
  All_KOs <- names(which(table(RNAseq_Annotated_Matrix$KO) >= N))[-no_annotation] #list of all KOs which apear greater than 5 times **
  Pairwise_Bin_Array_Presence <- matrix(0, nrow = dim_matrix,
                                           ncol =length(All_KOs))

  rownames(Pairwise_Bin_Array_Presence) <- names(table(RNAseq_Annotated_Matrix$Bin))[order(as.numeric(names(table(RNAseq_Annotated_Matrix$Bin))))]

  for (y in 1 : dim(Pairwise_Bin_Array_Presence)[2]) {
    # Fill in presence
    if (length(RNAseq_Annotated_Matrix$Bin[which(RNAseq_Annotated_Matrix$KO%in%All_KOs[y])]) >= N) {
      Pairwise_Bin_Array_Presence[which(rownames(Pairwise_Bin_Array_Presence)%in%RNAseq_Annotated_Matrix$Bin[which(RNAseq_Annotated_Matrix$KO%in%All_KOs[y])]),y] <-1
    }
  }

  return(Pairwise_Bin_Array_Presence)
}

