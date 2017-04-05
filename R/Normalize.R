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
#' @param method A string containing the method to use, either one of: ["default", "TMM", "RLE"].  In addition to the described default method, TMM and RLE from bioconductors edgeR
#' package are implemented as well
#' @export
#' @return The normalized read counts  of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Normalize(RNAseq_Annotated_Matrix, no_feature,ambiguous, not_aligned)
#' @note \preformatted{To remove rows that have a 0 for its read counts:}
#' \code{RNAseq_Annotated_Matrix[apply(RNAseq_Annotated_Matrix[, SS:SE], 1, function(x) !any(x == 0)), ]}
#' \preformatted{Where SS and SE are the start and end columns of the samples (raw counts).}
RNAseq_Normalize <- function(RNAseq_Annotated_Matrix, no_feature, ambiguous, 
			     not_aligned, method = "default"){

  SS<-2 # start column for samples
  SE<-length(sample_names) + 1 # end column of samples
  if(method == "default"){
    return(defaultRNA_Normalize(SS, SE, RNAseq_Annotated_Matrix, no_featureRNAseq_Annotated_Matrix, ambiguous,not_aligned))
    }
  else if(method == "TMM" || method == "RLE"){
    return(edgeRmethods(SS, SE, method, RNAseq_Annotated_Matrix))
    }
}

defaultRNA_Normalize <- function(SS, SE, RNAseq_Annotated_Matrix,no_featureRNAseq_Annotated_Matrix, ambiguous,not_aligned){
  # Calculate the number of reads mapped to each bin in each sample (This may be a separate function)
  sum_reads_per_genome_matrix<-matrix(NA,nrow=length(high_quality_bins),ncol=length(sample_names))
  for (i in 1:length(high_quality_bins)) {
    for (j in 2:(length(sample_names)+1)) {
      sum_reads_per_genome_matrix[i,j-1]<-sum(RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),j])
    }
  }

  # calculate max per bin.
  # Devide each column (sample) per row in normalized_sum_reads_per_genome_matrix (each bin) by the max count per bin
  normalized_sum_reads_per_genome_matrix<-sum_reads_per_genome_matrix/apply(sum_reads_per_genome_matrix,1,max)


  # normalize reads by max mapped to a genome
  for (i in 1:length(high_quality_bins)) {
    RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),SS:SE]<-RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),SS:SE]/normalized_sum_reads_per_genome_matrix[i,]
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

edgeRmethods <- function(SS, SE, method_name, RNAseq_Annotated_Matrix){
  library(edgeR)
  norm_factors <- calcNormFactors(RNAseq_Annotated_Matrix[, SS:SE], method=method_name)
  RNAseq_Annotated_Matrix[, SS:SE] <- as.matrix(RNAseq_Annotated_Matrix[,SS:SE]) * norm_factors
  return(RNAseq_Annotated_Matrix)
}




