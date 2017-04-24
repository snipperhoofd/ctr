#' Calculates pairwise composite PC & NRED scores for a module
#'
#' This function calculates all pairwise NRED & PC for a given module. When two genomes contain multiple
#' elements for comparison, this is handled either by takng a random pairwise comparison,
#' or by taking the score with the maximum similarity. A background distribution using each method
#' is provided for comparison. (Currently it only provides the score that minimizes the distances)
#'
#' @param Subset_KOs A list of KOs that form a module
#' @export
#' @return a list of vectors containing pairwise PC (pearsons), NRED (nred) and their composite Z-scores (Zscore),
#' and the position in the matrix for the highest scoring pair between genomes A and B
#' (positionA & positionB respectively)
#' @examples PHA_module_P_NRED <- P_NRED_Distance_Function(RNAseq_Annotated_Matrix_B, Z_scores_B, matrix_features_B, PHA_module)

P_NRED_Distance_Function <- function(RNAseq_Annotated_Matrix, Z_scores, matrix_features, Subset_KOs) {

  # Define two congruent arrays to be filled during the second step. Name the columns and rows based on the genome bins
  dim_matrix<-length(matrix_features@high_quality_bins)
  Pairwise_Bin_Array_Pearson<-array(NA,c(dim_matrix,dim_matrix, length(Subset_KOs)))
  colnames(Pairwise_Bin_Array_Pearson)<-names(table(RNAseq_Annotated_Matrix$Bin))[order(as.numeric(names(table(RNAseq_Annotated_Matrix$Bin))))]
  rownames(Pairwise_Bin_Array_Pearson)<-colnames(Pairwise_Bin_Array_Pearson)
  Pairwise_Bin_Array_Euclidean<-Pairwise_Bin_Array_Pearson
  Pairwise_PositionsA<-Pairwise_Bin_Array_Euclidean
  Pairwise_PositionsB<-Pairwise_Bin_Array_Euclidean
  Pairwise_Bin_Array_Presence	<- Presence_Absence_Matrix(RNAseq_Annotated_Matrix,5) # This should be removed from this function and stored globally, perhaps in matrix features


  ######### This is the maximum pairwise Pearson correlation between genome bins for all KOs, converted to Z score, and keeping the pair with the highest sum z-score

  for (x in 1:(dim(Pairwise_Bin_Array_Pearson)[1]-1)) {
    for (y in (x+1):dim(Pairwise_Bin_Array_Pearson)[2]) {
      for (z in 1:dim(Pairwise_Bin_Array_Pearson)[3]) { #iterate over array

        # Identify the rows in the original matrix corresponding to each genome
        position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_Bin_Array_Presence)[x])
        position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_Bin_Array_Presence)[y])
        # Second identify the rows in the original matrix corresponding to a KO
        position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]),position_of_genome_A)
        position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]),position_of_genome_B)
        # Make sure the KO is present in both genomes
        if (!length(position_of_kegg_enzyme_A)==0 && !length(position_of_kegg_enzyme_B)==0) {
          # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores
          # First define two empty matrices and then fill them with the PCC and Euc distances
          max_pairwise_gene_correlation<-matrix(NA,nrow=length(position_of_kegg_enzyme_A),ncol=length(position_of_kegg_enzyme_B))
          max_pairwise_gene_euclidean<-max_pairwise_gene_correlation
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
        } else {next}
      }
    }
  }

  Zscore_pairwise_gene_correlation<- ((Pairwise_Bin_Array_Pearson-Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
  Zscore_pairwise_gene_euclidean<- ((Pairwise_Bin_Array_Euclidean-Z_scores$sd[6]) / Z_scores$sd[6])
  Combined_Pairwise_Z_Score_Array<- ((-Zscore_pairwise_gene_correlation) + Zscore_pairwise_gene_euclidean)

  newList <- list("pearsons" = Zscore_pairwise_gene_correlation, 
                  "nred" = Zscore_pairwise_gene_euclidean,
                  "combined"=Combined_Pairwise_Z_Score_Array,
                  "positionsA"=Pairwise_PositionsA,
                  "positionsB"=Pairwise_PositionsB)

  return(newList)
}
