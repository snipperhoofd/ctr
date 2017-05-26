#' Calculating the background distribution for an indiviudal gene
#'
#' This function generates background distributions of the pairwise
#' similarity in the expression profile of 1) randomly selected genes and 2) randomly
#' selected genes that share the same function (annontation) in two different genomes.
#' It is common to have numerous genes with the same annotation in a genome.
#' This is handled with two approaches, in the first, a single randomly selected pair
#' of genes with the same function are selected. In the second case, the highest scoring
#' pair is used. To correct for multiple comparisons when the highest scoring pair is used
#' the random background distibution is recalculed to include an equal number of multiple
#' comparisons.
#'
#'
#' @param RNAseq_Annotation_Matrix_no_sd_of_zero, the original matrix with rows that have a standard deviation of zero removed.
#' @param N the number of iterations (at least 10000 is suggested)
#' @param range is a vector containing the minimum and maximum number of copies a KO must be found in the dataset to be used in calculating the background.
#' @export
#' @return A list containting eight vectors of N distances: pearson and euclidean distances for
#' random, multiple comparison corrected random, random KO, highest scoring pair KO.
#' @examples I_KOs_Background_B <- Individual_KOs_Background(transcriptional_responses$pairwise_KO_distances$combined,matrix_features_B,10000, language = 'R')
#' @examples I_KOs_Background_B <- Individual_KOs_Background(RNAseq_Annotation_Matrix_no_sd_of_zero_B,matrix_features_B,10000, language = 'C')

# first remove rows with standard deviations of 0


Pairwise_KOs_Calculation <- function(combined_z_scores_array, matrix_features){

  dim_matrix<- dim(combined_z_scores_array)[1]
  Pairwise_comparison_ave<- matrix(NA, dim_matrix, dim_matrix)
  rownames(Pairwise_comparison_ave) <- dimnames(combined_z_scores_array[1,,])[[1]]
  colnames(Pairwise_comparison_ave) <- dimnames(combined_z_scores_array[1,,])[[1]]
  Pairwise_comparison_length <- Pairwise_comparison_ave
  Pairwise_comparison_jaccard <- Pairwise_comparison_ave

# iterate over all pairwise combinations of genome bins
  for (i in 1:(dim_matrix-1)) {
    for (j in (i+1):dim_matrix) {
      values<- which(!is.na(combined_z_scores_array[i,j,])==TRUE)
      Pairwise_comparison_ave[i,j]<-ave(combined_z_scores_array[i,j,values])[1]
      Pairwise_comparison_length[i,j]<-length(values)
      Pairwise_comparison_jaccard[i,j]<-Calc_Jaccard(matrix_features@Pairwise_Bin_Array_Presence[i,],matrix_features@Pairwise_Bin_Array_Presence[j,])
    }
  }
  newList <- list("averages" = Pairwise_comparison_ave,
                  "lengths" =  Pairwise_comparison_length,
                  "jaccard" = Pairwise_comparison_jaccard)
  return (newList)
}
