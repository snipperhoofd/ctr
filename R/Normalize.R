#' Normalized RNAseq raw read counts
#'
#' RNAseq raw read counts may by normalized based on various parameters including reads per sample,
#' reads mapped per genome, gene length, log2 RPKM
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @param no_feature,ambiguous,not_aligned  A set of vectors equal to the number of samples,
#' containing the number of reads that had no feature,
#' where ambiguously mapped, or not aligned in their  (obtained from the mapping output).
#' @param gene_lengths A matrix with the length of each gene (genes must be in same order as input RNAseq_Annotated_Matrix)
#'
#' @return The normalized read counts  of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned)
#' @note \preformatted{To remove rows that have a 0 for its read counts:}
#' \code{RNAseq_Annotated_Matrix[apply(RNAseq_Annotated_Matrix[, SS:SE], 1, function(x) !any(x == 0)), ]}
#' \preformatted{Where SS and SE are the start and end columns of the samples (raw counts).}

RNAseq_Normalize <- function(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned){

  Bin_Column<-which(colnames(RNAseq_Annotated_Matrix) == "Bin")
  SS<-2 # start column for samples
  SE<-length(sample_names) + 1 # end column of samples
  RS<-Bin_Column + 1 # start column for ranks
  RE<-Bin_Column + length(sample_names) # end column for ranks

  # Calculate the number of reads mapped to each bin in each sample (This may be a separate function)
    sum_reads_per_genome_matrix<-matrix(NA,nrow=length(high_quality_bins),ncol=length(sample_names))
  for (i in 1:length(high_quality_bins)) {
    for (j in 2:(length(sample_names)+1)) {
      sum_reads_per_genome_matrix[i,j-1]<-sum(RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==names(table(RNAseq_Annotated_Matrix$Bin))[i]),j])
    }
  }

  # calculate max per bin
  normalized_sum_reads_per_genome_matrix<-sum_reads_per_genome_matrix/apply(sum_reads_per_genome_matrix,1,max)

  # normalize reads by max mapped to a genome
  for (i in 1:length(high_quality_bins)) {
    RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==names(table(RNAseq_Annotated_Matrix$Bin))[i]),SS:SE]<-RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==names(table(RNAseq_Annotated_Matrix$Bin))[i]),SS:SE]/normalized_sum_reads_per_genome_matrix[i,]
  }

  # normalized by total of non-rRNA reads per sample mapped
  sum_aligned<-apply(sum_reads_per_genome_matrix,2,sum)
  total_nonRNA_reads<-sum_aligned+no_feature+ ambiguous+ not_aligned
  normalized_by_total<- total_nonRNA_reads/max(total_nonRNA_reads)

  # Finale Normalization
  RNAseq_Annotated_Matrix[,SS:SE]<-RNAseq_Annotated_Matrix[,SS:SE]/normalized_by_total

  # convert to log base 2
  RNAseq_Annotated_Matrix[,SS:SE]<-log(RNAseq_Annotated_Matrix[,SS:SE],2)

  # replace -Inf with 0
  for (i in 2:(length(sample_names)+1)) {
    RNAseq_Annotated_Matrix[,i][is.infinite(RNAseq_Annotated_Matrix[,i])] <- 0
  }

  return(RNAseq_Annotated_Matrix)
}

