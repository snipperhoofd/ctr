#' Generate a background distibution of composite Z scores
#'
#' This function generates a background distibution based on randomly generated modules of size 'N'.
#' For each element of the module, pairwise Pearson Correlation (PC), and the Normalized Rank
#' Euclidean Distance (NRED) are calculated between all genomes. When two genomes contain multiple
#' elements for comparison, this is handled either by takng a random pairwise comparison,
#' or by taking the score with the maximum similarity. A background distribution using each method
#' is provided for comparison. (Currently it only provides the score that minimizes the distances)
#'
#' @param RNAseq_Annotated_Matrix The annotated matrix
#' @param N The number of elements to be included in the randomly generated module
#' @param Z The number of iterations used to calculatea background distribution


#' @export
#' @return a list of vectors containing XXX & YYY
#' @examples Random_Background_Module_Distances_6<-Background_Distribution_Modules(RNAseq_Annotated_Matrix,6,1000)

Background_Distribution_Modules <- function(RNAseq_Annotated_Matrix,matrix_features,Z_scores,N,Z) {

  Pairwise_Bin_Array_Presence	<- Presence_Absence_Matrix(RNAseq_Annotated_Matrix,5) # This should probably only be calculated once per dataset. Currently it is calculated in numerous functions. add another variable to replace 5 (e.g. The minimum number of times a KO term must be present to be included in the matrix)
  Random_Jaccard_Distances<-rep(NA,Z)
  Random_Composite_Distances<-rep(NA,Z)
  Random_Pearson_Distances<-rep(NA,N)
  Random_Euclidean_Distances<-rep(NA,N)
  Random_Zscore_Pearson_Distances<-rep(NA,N)
  Random_Zscore_Euclidean_Distances<-rep(NA,N)
  All_KOs<-names(which(table(RNAseq_Annotated_Matrix$KO)>=5))[-1] # This was originally a global variable but was moved so that it can change depending on the annotation matrix used

  # iterate Z times
  for (i in 1:Z) {


    random_genomes<- sample(length(matrix_features@high_quality_bins), 2)
    Random_Module<- Generate_Random_Module(All_KOs,N)
    All_position_KOs<- which(All_KOs%in%Random_Module)
    # Calculate Jaccard Distance
    PA_position_of_genome_A<-which(rownames(Pairwise_Bin_Array_Presence)==matrix_features@high_quality_bins[random_genomes[1]])
    PA_position_of_genome_B<-which(rownames(Pairwise_Bin_Array_Presence)==matrix_features@high_quality_bins[random_genomes[2]])

    Random_Jaccard_Distances[i]<-Calc_Jaccard(Pairwise_Bin_Array_Presence[PA_position_of_genome_A,
                                                                          All_position_KOs],
                                              Pairwise_Bin_Array_Presence[PA_position_of_genome_B,
                                                                          All_position_KOs])

    # Next calculate Pearson and NRED
    for (j in 1:N) {
      # Identify the rows in the original matrix corresponding to each genome
      position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[random_genomes[1]])
      position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[random_genomes[2]])
      # Identify the rows in the original matrix corresponding to a KO
      KO_positions<- which(RNAseq_Annotated_Matrix$KO==Random_Module[j])
      # Find intersection between the genome lists and the KO list
      position_of_kegg_enzyme_A = intersect(KO_positions,position_of_genome_A)
      position_of_kegg_enzyme_B = intersect(KO_positions,position_of_genome_B)
      # Make sure the KO is present in both genomes
      l_position_of_kegg_enzyme_A<-length(position_of_kegg_enzyme_A)
      l_position_of_kegg_enzyme_B<-length(position_of_kegg_enzyme_B)
      # check if both are zero by multiplication
      if (!(l_position_of_kegg_enzyme_A * l_position_of_kegg_enzyme_B)==0) { # may be
        # Define two empty matrices and then fill them with the PCC and Euc distances
        max_pairwise_gene_correlation<-matrix(NA,nrow=l_position_of_kegg_enzyme_A,ncol=l_position_of_kegg_enzyme_B)
        max_pairwise_gene_euclidean<-max_pairwise_gene_correlation
        # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances
        for (m in 1:l_position_of_kegg_enzyme_A){
          for (n in 1:l_position_of_kegg_enzyme_B){
            # use the no_sd dataset so that calculating a Pearson correlation never gives an error
              max_pairwise_gene_correlation[m,n]<- cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                                                          matrix_features@SS:matrix_features@SE]),
                                                       as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                                                          matrix_features@SS:matrix_features@SE]))

              subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@RS:matrix_features@RE] -
                                 RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@RS:matrix_features@RE]
              max_pairwise_gene_euclidean[m, n]<- sqrt(sum(subtracted_lists * subtracted_lists))
          }
        } # End Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances

        # Convert to Z scores
        Zscore_pairwise_gene_correlation<- (max_pairwise_gene_correlation-Z_scores$mu[2])/Z_scores$sd[2] # need to inverse PCC
        Zscore_pairwise_gene_euclidean<- (max_pairwise_gene_euclidean-Z_scores$mu[6])/Z_scores$sd[6]
        best_scoring_pair<- which.min((1-Zscore_pairwise_gene_correlation)+(Zscore_pairwise_gene_euclidean))
        if (length(best_scoring_pair)>0) {
          Random_Pearson_Distances<-max_pairwise_gene_correlation[best_scoring_pair]
          Random_Euclidean_Distances<-max_pairwise_gene_euclidean[best_scoring_pair]
        } else {
          scoring_pair<-which.min(Zscore_pairwise_gene_euclidean)
          Random_Pearson_Distances<-max_pairwise_gene_correlation[scoring_pair]
          Random_Euclidean_Distances<-max_pairwise_gene_euclidean[scoring_pair]
        }
        # if one genomes does not contain the KO (as checked by multiplication above)
      } else {
        Random_Pearson_Distances<- NA
        Random_Euclidean_Distances<- NA
      }
      Random_Zscore_Pearson_Distances[j]<-((Random_Pearson_Distances-Z_scores$mu[2])/Z_scores$sd[2]) # Need to input the Z_scores matrix
      Random_Zscore_Euclidean_Distances[j]<-((Random_Euclidean_Distances-Z_scores$mu[6])/Z_scores$sd[6])

    }
    Random_Composite_Distances[i]<-mean((-Random_Zscore_Pearson_Distances)+Random_Zscore_Euclidean_Distances,na.rm=TRUE)[1]
    Random_Background_Module_Distances<-Random_Composite_Distances*(1-Random_Jaccard_Distances)
  }
  return(Random_Background_Module_Distances)
}

dothestuff <- function(RNAseq_Annotated_Matrix,matrix_features,Z_scores,N,i){


  random_genomes<- sample(length(matrix_features@high_quality_bins), 2)
  Random_Module<- Generate_Random_Module(matrix_features@All_KOs, N)
  All_position_KOs<- which(matrix_features@All_KOs%in%Random_Module)
  # Calculate Jaccard Distance
  PA_position_of_genome_A<-which(rownames(matrix_features@Pairwise_Bin_Array_Presence)==matrix_features@high_quality_bins[random_genomes[1]])
  PA_position_of_genome_B<-which(rownames(matrix_features@Pairwise_Bin_Array_Presence)==matrix_features@high_quality_bins[random_genomes[2]])

  Random_Jaccard_Distances[i]<-Calc_Jaccard(matrix_features@Pairwise_Bin_Array_Presence[PA_position_of_genome_A,
                                                                                        All_position_KOs],
                                            matrix_features@Pairwise_Bin_Array_Presence[PA_position_of_genome_B,
                                                                                        All_position_KOs])

  # Next calculate Pearson and NRED
  for (j in 1:N) {
    # Identify the rows in the original matrix corresponding to each genome
    position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[random_genomes[1]])
    position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==matrix_features@high_quality_bins[random_genomes[2]])
    # Identify the rows in the original matrix corresponding to a KO
    KO_positions<- which(RNAseq_Annotated_Matrix$KO==Random_Module[j])
    # Find intersection between the genome lists and the KO list
    position_of_kegg_enzyme_A = intersect(KO_positions,position_of_genome_A)
    position_of_kegg_enzyme_B = intersect(KO_positions,position_of_genome_B)
    # Make sure the KO is present in both genomes
    l_position_of_kegg_enzyme_A<-length(position_of_kegg_enzyme_A)
    l_position_of_kegg_enzyme_B<-length(position_of_kegg_enzyme_B)
    # check if both are zero by multiplication
    if (!(l_position_of_kegg_enzyme_A * l_position_of_kegg_enzyme_B)==0) { # may be
      # Define two empty matrices and then fill them with the PCC and Euc distances
      max_pairwise_gene_correlation<-matrix(NA,nrow=l_position_of_kegg_enzyme_A,ncol=l_position_of_kegg_enzyme_B)
      max_pairwise_gene_euclidean<-max_pairwise_gene_correlation
      # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances
      for (m in 1:l_position_of_kegg_enzyme_A){
        for (n in 1:l_position_of_kegg_enzyme_B){
          # use the no_sd dataset so that calculating a Pearson correlation never gives an error
          max_pairwise_gene_correlation[m,n]<- cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                                                      matrix_features@SS:matrix_features@SE]),
                                                   as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                                                      matrix_features@SS:matrix_features@SE]))

          subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@RS:matrix_features@RE] -
            RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@RS:matrix_features@RE]
          max_pairwise_gene_euclidean[m, n]<- sqrt(sum(subtracted_lists * subtracted_lists))
        }
      } # End Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances

      # Convert to Z scores
      Zscore_pairwise_gene_correlation<- (max_pairwise_gene_correlation-Z_scores$mu[2])/Z_scores$sd[2] # need to inverse PCC
      Zscore_pairwise_gene_euclidean<- (max_pairwise_gene_euclidean-Z_scores$mu[6])/Z_scores$sd[6]
      best_scoring_pair<- which.min((1-Zscore_pairwise_gene_correlation)+(Zscore_pairwise_gene_euclidean))
      if (length(best_scoring_pair)>0) {
        Random_Pearson_Distances<-max_pairwise_gene_correlation[best_scoring_pair]
        Random_Euclidean_Distances<-max_pairwise_gene_euclidean[best_scoring_pair]
      } else {
        scoring_pair<-which.min(Zscore_pairwise_gene_euclidean)
        Random_Pearson_Distances<-max_pairwise_gene_correlation[scoring_pair]
        Random_Euclidean_Distances<-max_pairwise_gene_euclidean[scoring_pair]
      }
      # if one genomes does not contain the KO (as checked by multiplication above)
    } else {
      Random_Pearson_Distances<- NA
      Random_Euclidean_Distances<- NA
    }
    Random_Zscore_Pearson_Distances[j]<-((Random_Pearson_Distances-Z_scores$mu[2])/Z_scores$sd[2]) # Need to input the Z_scores matrix
    Random_Zscore_Euclidean_Distances[j]<-((Random_Euclidean_Distances-Z_scores$mu[6])/Z_scores$sd[6])
  }
  Random_Composite_Distances[i]<-mean((-Random_Zscore_Pearson_Distances)+Random_Zscore_Euclidean_Distances,na.rm=TRUE)[1]
  Random_Background_Module_Distances<-Random_Composite_Distances*(1-Random_Jaccard_Distances)
}

bgModule_par <- function(RNAseq_Annotated_Matrix,matrix_features,Z_scores,N,Z) {


  Random_Jaccard_Distances<-rep(NA,Z)
  Random_Composite_Distances<-rep(NA,Z)
  Random_Pearson_Distances<-rep(NA,N)
  Random_Euclidean_Distances<-rep(NA,N)
  Random_Zscore_Pearson_Distances<-rep(NA,N)
  Random_Zscore_Euclidean_Distances<-rep(NA,N)


  # iterate Z times
  for (i in 1:Z) {
    dothestuff(RNAseq_Annotated_Matrix,matrix_features,Z_scores,N,i)
  }
}
