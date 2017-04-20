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


  for (i in 1:matrix_features@sample_size) {
    Rank_column_name<-paste("rank",i,sep="")
    RNAseq_Annotated_Matrix[,matrix_features@Bin_Column+i]<-c(NA)
  }

#  colnames(RNAseq_Annotated_Matrix)<-c("Locus_ID",matrix_features@sample_names,"KO","Bin",Rank_column_name)

  # add to the matrix the ranks
  for (i in 1:matrix_features@sample_size) {
    for (s in 1:length(matrix_features@high_quality_bins)) {
      matrix_HQ_bins <- RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[s]),]
      matrix_HQ_bins[,matrix_features@Bin_Column+i]<-rank(-matrix_HQ_bins[,i+1], na.last=TRUE, ties.method="random")/max(rank(-matrix_HQ_bins[,i+1], na.last=TRUE, ties.method="random"))
      RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),] <- matrix_HQ_bins
    }

  }
  return(RNAseq_Annotated_Matrix)
}


#########
#Testing
#########
Create_Rank_Columns2 <- function(RNAseq_Annotated_Matrix, matrix_features){


  for (i in 1:matrix_features@sample_size) {
    Rank_column_name<-paste("rank",i,sep="")
    RNAseq_Annotated_Matrix[,matrix_features@Bin_Column+i]<-c(NA)
  }

  #  colnames(RNAseq_Annotated_Matrix)<-c("Locus_ID",matrix_features@sample_names,"KO","Bin",Rank_column_name)

  # add to the matrix the ranks
  # for (i in 1:matrix_features@sample_size) {
  #   for (s in 1:length(matrix_features@high_quality_bins)) {
  #     matrix_HQ_bins <- RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[s]),]
  #     ranks <- rank(-matrix_HQ_bins[,i+1], na.last=TRUE, ties.method="random")
  #     matrix_HQ_bins[,matrix_features@Bin_Column+i]<-ranks/max(ranks)
  #     RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),] <- matrix_HQ_bins
  #   }
  #
  # }

  sample_cols <- matrix_features@SS:matrix_features@SE
  data_matrix <- matrix(apply(RNAseq_Annotated_Matrix[, sample_cols], 2, as.numeric), ncol = length(sample_cols))
  vals <- GetRanksPerBin(data_matrix,
                 sapply(RNAseq_Annotated_Matrix[, matrix_features@Bin_Column], as.numeric),
                 matrix_features@high_quality_bins)

  return(vals)
}
