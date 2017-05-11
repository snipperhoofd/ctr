#' Calculates pairwise composite PC & NRED scores for a module
#'
#' This function calculates all pairwise NRED & PC for all genes, or a subset of of genes in a given module.
#' When two genomes contain multiple elements for comparison, this is handled either by takng the
#' score with the maximum similarity.
#'
#' @param RNAseq_Annotated_Matrix The normalized count data with annotations
#' @param Z_scores The Zscores obtained from the background distributions
#' @param matrix_features The matrix features associated with the annoated normalized count data
#' @param Subset_KOs A vector of the KOs that form a module
#'
#' @export
#' @return a list of arrays including pairwise PC ($pearsons), NRED ($nred), the row for genome A and B
#' used in the calculation ($positionA and $positionB), the Z-scores ($Zscore_pearson and $Zscore_nred)
#' and a composite Z-score ($combined). Each array has the dimensions #bins x #bins x #KOs.
#'
#' @examples
#' PHA_module_P_NRED <- P_NRED_Distance_Function(RNAseq_Annotated_Matrix_BR_default_bin, Z_scores_B, matrix_features_B, PHA_module)
#' All_KO_module_P_NRED_test <- P_NRED_Distance_Function(RNAseq_Annotated_Matrix_BR_default_bin, Z_scores_B, matrix_features_B)

P_NRED_Distance_Function <- function(RNAseq_Annotated_Matrix, Z_scores, matrix_features, Subset_KOs) {

  # if no subset of KOs is introduced, than by default, all calcuations are conducted
  if(missing(Subset_KOs)) {
    no_annotation <- which(names((table(RNAseq_Annotated_Matrix$KO))) == "")
    Subset_KOs<- names((table(RNAseq_Annotated_Matrix$KO)))[-no_annotation]
  }

  # build empty arrays, label columns and rows, for each:
  # Pearson Correlation, NRED and the position of genome A and B

  dim_matrix<- length(matrix_features@high_quality_bins)
  # Pearson
  Pairwise_Bin_Array_Pearson<-array(data=NA,
                                    dim = c(dim_matrix, dim_matrix, length(Subset_KOs)),
                                    dimnames = list(sort(matrix_features@high_quality_bins),
                                                    sort(matrix_features@high_quality_bins)))

 # NRED
  Pairwise_Bin_Array_Euclidean<- Pairwise_Bin_Array_Pearson
  # Positions
  Pairwise_PositionsA<- Pairwise_Bin_Array_Euclidean
  Pairwise_PositionsB<- Pairwise_Bin_Array_Euclidean

  # save names in variable to reduce the number of functions called
  bin_names<- rownames(Pairwise_Bin_Array_Pearson)

  ######### This is the maximum pairwise Pearson correlation between genome bins for all KOs,
  ######### converted to Z score, and keeping the pair with the highest sum z-score

  for (x in 1:(dim(Pairwise_Bin_Array_Pearson)[1]-1)) { # Iterate over the Pairwise bin array collumns
    for (y in (x+1):dim(Pairwise_Bin_Array_Pearson)[2]) {
      for (z in 1:dim(Pairwise_Bin_Array_Pearson)[3]) { #iterate over array

        # Identify the rows in the original matrix corresponding to each genome
        position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==bin_names[x])
        position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==bin_names[y])
        # Second identify the rows in the original matrix corresponding to a KO
        position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]), position_of_genome_A)
        position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]), position_of_genome_B)

        # Make sure the KO is present in both genomes
        if (!length(position_of_kegg_enzyme_A) == 0 && !length(position_of_kegg_enzyme_B) == 0) {
          # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores

          # First define two empty matrices and then fill them with the PCC and NRED distances (This could be moved to improve speed)
          max_pairwise_gene_correlation<- matrix(NA,
                                                nrow=length(position_of_kegg_enzyme_A),
                                                ncol=length(position_of_kegg_enzyme_B))
          max_pairwise_gene_euclidean<- max_pairwise_gene_correlation

            # loop through the possible KO pairs
            for (m in 1:length(position_of_kegg_enzyme_A)){
              for (n in 1:length(position_of_kegg_enzyme_B)){
                # make sure there is always a standard deviation, or else cor gives an error. If there is a sd, proceed with calculations
                if (sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                          matrix_features@SS:matrix_features@SE]))!=0 &&
                    sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                          matrix_features@SS:matrix_features@SE]))!=0) {

                  max_pairwise_gene_correlation[m,n]<- (cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                                                               matrix_features@SS:matrix_features@SE]),
                                                            as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                                                               matrix_features@SS:matrix_features@SE])))
                  subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                             matrix_features@RS:matrix_features@RE] -
                                     RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                             matrix_features@RS:matrix_features@RE]

                  max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
                } else {
                  # If there is no standard deviation, the correlation is NA
                  max_pairwise_gene_correlation[m,n]<-NA
                  subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                             matrix_features@RS:matrix_features@RE] -
                                     RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                             matrix_features@RS:matrix_features@RE]

                  max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
              }
            }
          }
          # Convert to Z scores
          Zscore_pairwise_gene_correlation<-((max_pairwise_gene_correlation-Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
          Zscore_pairwise_gene_euclidean<-((max_pairwise_gene_euclidean-Z_scores$mu[6]) / Z_scores$mu[6])

          best_scoring_pair<-which.min((1-max_pairwise_gene_correlation) + (Zscore_pairwise_gene_euclidean))
          rownames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_A
          colnames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_B

          if (length(best_scoring_pair)>0) {
            Pairwise_Bin_Array_Pearson[x,y,z]<- max_pairwise_gene_correlation[best_scoring_pair]
            Pairwise_Bin_Array_Euclidean[x,y,z]<- max_pairwise_gene_euclidean[best_scoring_pair]
            Pairwise_PositionsA[x,y,z]<- rownames(max_pairwise_gene_correlation)[best_scoring_pair]
            Pairwise_PositionsB[x,y,z]<- colnames(max_pairwise_gene_correlation)[best_scoring_pair]
          } else {
            scoring_pair<- which.min(max_pairwise_gene_euclidean)
            Pairwise_Bin_Array_Pearson[x,y,z]<- max_pairwise_gene_correlation[scoring_pair]
            Pairwise_Bin_Array_Euclidean[x,y,z]<- max_pairwise_gene_euclidean[scoring_pair]
            Pairwise_PositionsA[x,y,z]<- rownames(max_pairwise_gene_correlation)[scoring_pair]
            Pairwise_PositionsB[x,y,z]<- colnames(max_pairwise_gene_correlation)[scoring_pair]}
        }
        else {next}
        #  print(c(x,y,z)) # Unhash to monitor progress
      }
    }
  }

  Zscore_pairwise_gene_correlation<- ((Pairwise_Bin_Array_Pearson-Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
  Zscore_pairwise_gene_euclidean<- ((Pairwise_Bin_Array_Euclidean-Z_scores$sd[6]) / Z_scores$sd[6])
  Combined_Pairwise_Z_Score_Array<- ((-Zscore_pairwise_gene_correlation) + Zscore_pairwise_gene_euclidean)

  newList <- list("pearsons" = Pairwise_Bin_Array_Pearson,
                  "nred" =  Pairwise_Bin_Array_Euclidean,
                  "combined" = Combined_Pairwise_Z_Score_Array,
                  "Zscore_pearson" = Zscore_pairwise_gene_correlation,
                  "Zscore_nred" = Zscore_pairwise_gene_euclidean,
                  "positionsA" = Pairwise_PositionsA,
                  "positionsB" = Pairwise_PositionsB)

  return(newList)
}

P_NRED_Distance_Function <- function(RNAseq_Annotated_Matrix, Z_scores, matrix_features, Subset_KOs) {

  # if no subset of KOs is introduced, than by default, all calcuations are conducted
  if(missing(Subset_KOs)) {
    no_annotation <- which(names((table(RNAseq_Annotated_Matrix$KO))) == "")
    Subset_KOs<- names((table(RNAseq_Annotated_Matrix$KO)))[-no_annotation]
  }

  # build empty arrays, label columns and rows, for each:
  # Pearson Correlation, NRED and the position of genome A and B

  dim_matrix<- length(matrix_features@high_quality_bins)
  # Pearson
  Pairwise_Bin_Array_Pearson<-array(data=NA,
                                    dim = c(dim_matrix, dim_matrix, length(Subset_KOs)),
                                    dimnames = list(sort(matrix_features@high_quality_bins),
                                                    sort(matrix_features@high_quality_bins)))

  # NRED
  Pairwise_Bin_Array_Euclidean<- Pairwise_Bin_Array_Pearson
  # Positions
  Pairwise_PositionsA<- Pairwise_Bin_Array_Euclidean
  Pairwise_PositionsB<- Pairwise_Bin_Array_Euclidean

  # save names in variable to reduce the number of functions called
  bin_names<- rownames(Pairwise_Bin_Array_Pearson)

  ######### This is the maximum pairwise Pearson correlation between genome bins for all KOs,
  ######### converted to Z score, and keeping the pair with the highest sum z-score

  for (x in 1:(dim(Pairwise_Bin_Array_Pearson)[1]-1)) { # Iterate over the Pairwise bin array collumns
    for (y in (x+1):dim(Pairwise_Bin_Array_Pearson)[2]) {
      for (z in 1:dim(Pairwise_Bin_Array_Pearson)[3]) { #iterate over array

        # Identify the rows in the original matrix corresponding to each genome
        position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==bin_names[x])
        position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==bin_names[y])
        # Second identify the rows in the original matrix corresponding to a KO
        position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]), position_of_genome_A)
        position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]), position_of_genome_B)

        # Make sure the KO is present in both genomes
        if (!length(position_of_kegg_enzyme_A) == 0 && !length(position_of_kegg_enzyme_B) == 0) {
          # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores

          # First define two empty matrices and then fill them with the PCC and NRED distances (This could be moved to improve speed)
          max_pairwise_gene_correlation<- matrix(NA,
                                                 nrow=length(position_of_kegg_enzyme_A),
                                                 ncol=length(position_of_kegg_enzyme_B))
          max_pairwise_gene_euclidean<- max_pairwise_gene_correlation

          # loop through the possible KO pairs
          for (m in 1:length(position_of_kegg_enzyme_A)){
            for (n in 1:length(position_of_kegg_enzyme_B)){
              # make sure there is always a standard deviation, or else cor gives an error. If there is a sd, proceed with calculations
              if (sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                        matrix_features@SS:matrix_features@SE]))!=0 &&
                  sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                        matrix_features@SS:matrix_features@SE]))!=0) {

                max_pairwise_gene_correlation[m,n]<- (cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                                                             matrix_features@SS:matrix_features@SE]),
                                                          as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                                                             matrix_features@SS:matrix_features@SE])))
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                           matrix_features@RS:matrix_features@RE] -
                  RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                          matrix_features@RS:matrix_features@RE]

                max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
              } else {
                # If there is no standard deviation, the correlation is NA
                max_pairwise_gene_correlation[m,n]<-NA
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                           matrix_features@RS:matrix_features@RE] -
                  RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                          matrix_features@RS:matrix_features@RE]

                max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
              }
            }
          }
          # Convert to Z scores
          Zscore_pairwise_gene_correlation<-((max_pairwise_gene_correlation-Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
          Zscore_pairwise_gene_euclidean<-((max_pairwise_gene_euclidean-Z_scores$mu[6]) / Z_scores$mu[6])

          best_scoring_pair<-which.min((1-max_pairwise_gene_correlation) + (Zscore_pairwise_gene_euclidean))
          rownames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_A
          colnames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_B

          if (length(best_scoring_pair)>0) {
            Pairwise_Bin_Array_Pearson[x,y,z]<- max_pairwise_gene_correlation[best_scoring_pair]
            Pairwise_Bin_Array_Euclidean[x,y,z]<- max_pairwise_gene_euclidean[best_scoring_pair]
            Pairwise_PositionsA[x,y,z]<- rownames(max_pairwise_gene_correlation)[best_scoring_pair]
            Pairwise_PositionsB[x,y,z]<- colnames(max_pairwise_gene_correlation)[best_scoring_pair]
          } else {
            scoring_pair<- which.min(max_pairwise_gene_euclidean)
            Pairwise_Bin_Array_Pearson[x,y,z]<- max_pairwise_gene_correlation[scoring_pair]
            Pairwise_Bin_Array_Euclidean[x,y,z]<- max_pairwise_gene_euclidean[scoring_pair]
            Pairwise_PositionsA[x,y,z]<- rownames(max_pairwise_gene_correlation)[scoring_pair]
            Pairwise_PositionsB[x,y,z]<- colnames(max_pairwise_gene_correlation)[scoring_pair]}
        }
        else {next}
        #  print(c(x,y,z)) # Unhash to monitor progress
      }
    }
  }

  Zscore_pairwise_gene_correlation<- ((Pairwise_Bin_Array_Pearson-Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
  Zscore_pairwise_gene_euclidean<- ((Pairwise_Bin_Array_Euclidean-Z_scores$sd[6]) / Z_scores$sd[6])
  Combined_Pairwise_Z_Score_Array<- ((-Zscore_pairwise_gene_correlation) + Zscore_pairwise_gene_euclidean)

  newList <- list("pearsons" = Pairwise_Bin_Array_Pearson,
                  "nred" =  Pairwise_Bin_Array_Euclidean,
                  "combined" = Combined_Pairwise_Z_Score_Array,
                  "Zscore_pearson" = Zscore_pairwise_gene_correlation,
                  "Zscore_nred" = Zscore_pairwise_gene_euclidean,
                  "positionsA" = Pairwise_PositionsA,
                  "positionsB" = Pairwise_PositionsB)

  return(newList)
}
