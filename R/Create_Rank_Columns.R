#' Calculate Normalized Ranks (e.g. rank / total genes)
#'
#' This function calculates the Normalized Ranks (0-1 scale) of each gene in each genome.
#' The ranks must be normalized so that they may be compared between genomes of different size.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @export
#' @return The normalized ranks of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)


Create_Rank_Columns <- function(RNAseq_Annotated_Matrix, matrix_features){

  Bin_Column <- which(colnames(RNAseq_Annotated_Matrix) == "Bin")
  Rank_column_name <- rep(NA, matrix_features@sample_size)
  for (i in 1 : matrix_features@sample_size) {
    Rank_column_name[i] <- paste("Rank", i, sep = "")
    RNAseq_Annotated_Matrix[, Bin_Column + i] <- c(NA)
  }

  colnames(RNAseq_Annotated_Matrix)<-c("Locus_ID", matrix_features@sample_names, "KO", "Bin", Rank_column_name)

  # add to the matrix the ranks
  for (i in 1:matrix_features@sample_size) {
    for (s in 1:length(matrix_features@high_quality_bins)) {
      matrix_HQ_bins <- RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),]
      matrix_HQ_bins[,Bin_Column+i]<-rank(-matrix_HQ_bins[,i+1], na.last=TRUE, ties.method="random")/max(rank(-matrix_HQ_bins[,i+1], na.last=TRUE, ties.method="random"))
      RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),] <- matrix_HQ_bins
    }
  }
  return(RNAseq_Annotated_Matrix)
}
