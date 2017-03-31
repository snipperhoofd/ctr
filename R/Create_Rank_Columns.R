#' Calculate Normalized Ranks (e.g. rank / total genes)
#'
#' This function calculates the Normalized Ranks (0-1 scale) of each gene in each genome.
#' The ranks must be normalized so that they may be compared between genomes of different size.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @export
#' @return The normalized ranks of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)


Create_Rank_Columns <- function(RNAseq_Annotated_Matrix){

  Bin_Column<-which(colnames(RNAseq_Annotated_Matrix)=="Bin")
  #SS<-2
  #SE<-length(sample_names)+1
  #RS<-Bin_Column+1
  #RE<-Bin_Column+length(sample_names)
  for (i in 1:length(sample_names)) {
    Rank_column_name<-paste("rank",i,sep="")
    RNAseq_Annotated_Matrix[,Bin_Column+i]<-c(NA)
  }

  Rank_column_name<-rep(NA,length(sample_names))
  for (i in 1:length(sample_names)) {Rank_column_name[i]<-paste("Rank",i,sep="")}

  colnames(RNAseq_Annotated_Matrix)<-c("Locus_ID",sample_names,"KO","Bin",Rank_column_name)

  # add to the matrix the ranks
  for (i in 1:length(sample_names)) {
    for (s in 1:length(high_quality_bins)) {
      RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),Bin_Column+i]<-rank(-RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),i+1], na.last=TRUE, ties.method="random")/max(rank(-RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[s]),i+1], na.last=TRUE, ties.method="random"))
    }
  }
  return(RNAseq_Annotated_Matrix)
}
