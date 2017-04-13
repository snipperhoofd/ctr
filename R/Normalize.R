#' Normalized RNAseq raw read counts
#'
#' RNAseq raw read counts may by normalized based on various parameters including reads per sample,
#' reads mapped per genome, gene length etc. Here we implement the normalization of edgeR (citation)
#' which accounts for differences in both Sequencing Detph and RNA composition (see edgeR documentation 
#' page 2.7.2 & 2.7.3). However, in metatranscritpomic studies, it may also be beneficial to 
#' adjust for an additional source of compositional bias in which a single organisms may contribute
#' a high relative abundance of transcrtipts, resulting in an undersampling of other organisms. 
#' Therefore, we provide an additional normalization step to normalize on a per genome/bin basis.
#'
#' @param RNAseq_Annotated_Matrix The original count matrix (See X for format details).
#' @param matrix_features General_features object that contains various metrics about the RNAseq data in the matrix.
#'                          the values SS, SE and sample_names are used, together with  no_feature, ambiguous and not_aligned. A set of vectors equal to the number of samples,
#' containing the number of reads that had no feature,
#' where ambiguously mapped, or not aligned in their  (obtained from the mapping output).
#' @param gene_lengths <NOT IMPLEMENTED> A matrix with the length of each gene (genes must be in same order as input RNAseq_Annotated_Matrix)
#' @param method A string containing the method to use, either one of: ["default", "TMM", "RLE"].  In addition to the described default method, TMM and RLE from bioconductors edgeR
#' package are implemented as well
#' @export
#' @return The normalized read counts  of \code{Sample 1} ... \code{Sample N}.
#' @examples RNAseq_Normalize(RNAseq_Annotated_Matrix, no_feature,ambiguous, not_aligned)
#' @note \preformatted{To remove rows that have a 0 for its read counts:}
#' \code{RNAseq_Annotated_Matrix[apply(RNAseq_Annotated_Matrix[, SS:SE], 1, function(x) !any(x == 0)), ]}
#' \preformatted{Where SS and SE are the start and end columns of the samples (raw counts).}

RNAseq_Normalize <- function(RNAseq_Annotated_Matrix, matrix_features, method = "default2"){


  if(method == "default"){

    return(defaultRNA_Normalize(RNAseq_Annotated_Matrix, matrix_features))

    }
  else if(method == "TMM" || method == "RLE"){
    return(edgeRmethods(RNAseq_Annotated_Matrix, method, matrix_features))
  }
  else if(method == "default2"){
    return(test_method(RNAseq_Annotated_Matrix, matrix_features))
  }
}


defaultRNA_Normalize <- function(SS, SE, RNAseq_Annotated_Matrix,no_feature, ambiguous,not_aligned){
  
  # normalized by total of non-rRNA reads per sample mapped
  sum_aligned<-apply(sum_reads_per_genome_matrix,2,sum)
  total_nonRNA_reads<-sum_aligned+no_feature+ ambiguous+ not_aligned
  normalized_by_total<- total_nonRNA_reads/max(total_nonRNA_reads) 
  
  # An alternative would be to use the library size. 
  # However, since the pool of rRNA is often both physically and bioinformatically removed, 
  # the count of mRNA reads per sample is more intuitively relavent.
  # normalized_by_total <- library_size/max(library_size)
  
  # Finale Normalization
  RNAseq_Annotated_Matrix[,SS:SE]<-t(t(RNAseq_Annotated_Matrix[,SS:SE])/normalized_by_total)

  # convert to log base 2
  # RNAseq_Annotated_Matrix[,SS:SE]<-log(RNAseq_Annotated_Matrix[,SS:SE],2)


  # replace -Inf with 0
  for (i in 2 : matrix_features@SE) {
    RNAseq_Annotated_Matrix[,i][is.infinite(RNAseq_Annotated_Matrix[,i])] <- 0
  }
  return(RNAseq_Annotated_Matrix)
}


edgeRmethods <- function(RNAseq_Annotated_Matrix, method_name, matrix_features){
  library(edgeR)


  RNAseq_Annotated_Matrix1<-DGEList(counts=RNAseq_Annotated_Matrix[, SS:SE],lib.size=library_size)
  norm_factors <- calcNormFactors(RNAseq_Annotated_Matrix1, method=method_name)
  RNAseq_Annotated_Matrix[, SS:SE] <- t(t(as.matrix(norm_factors)) / norm_factors$samples[,3])
  return(RNAseq_Annotated_Matrix)
}

Normalize_by_bin <- function(RNAseq_Annotated_Matrix){
  
  # Step 1: Calculate the number of reads mapped to each bin in each sample (This may be a separate function)
  sum_reads_per_genome_matrix<-matrix(NA,nrow=length(high_quality_bins),ncol=length(sample_names))
  for (i in 1:length(high_quality_bins)) {
    for (j in 2:(length(sample_names)+1)) {
      sum_reads_per_genome_matrix[i,j-1]<-sum(RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),j])
    }
  }

  # Step 2: Calculate max per bin.
  # Divide each column (sample) per row in normalized_sum_reads_per_genome_matrix (each bin) by 
  # the max count per bin
  
  normalized_sum_reads_per_genome_matrix<-sum_reads_per_genome_matrix/apply(sum_reads_per_genome_matrix,1,max)

  # Step 3: normalize reads by max mapped to a genome
  for (i in 1:length(high_quality_bins)) {
  RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),SS:SE]<-t(t(RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[i]),SS:SE])/normalized_sum_reads_per_genome_matrix[i,])
  }

return(RNAseq_Annotated_Matrix)
}

  for (i in 1:length(matrix_features@high_quality_bins)) {
    for (j in 2:matrix_features@SE) {
      sum_reads_per_genome_matrix[i,j-1]<-sum(RNAseq_Annotated_Matrix[which(
        RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[i]),j])
    }
  }
  return(sum_reads_per_genome_matrix)
}


