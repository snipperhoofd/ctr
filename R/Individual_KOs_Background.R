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
#' @export
#' @return A list containting eight vectors of N distances: pearson and euclidean distances for
#' random, multiple comparison corrected random, random KO, highest scoring pair KO.
#' @examples I_KOs_Background_B <- Individual_KOs_Background(RNAseq_Annotation_Matrix_no_sd_of_zero_B,matrix_features_B,10000, language = 'R')
#' @examples I_KOs_Background_B <- Individual_KOs_Background(RNAseq_Annotation_Matrix_no_sd_of_zero_B,matrix_features_B,10000, language = 'C')

# first remove rows with standard deviations of 0


Individual_KOs_Background <- function(RNAseq_Annotation_Matrix_no_sd_of_zero, matrix_features, N, language = 'C'){

  # build empty vectors for each Pearson Correlation and NRED
  random_pairwise_gene_pearson<- rep(NA, N)
  H_random_pairwise_gene_pearson<- rep(NA, N)
  random_pairwise_gene_euclidean<- rep(NA, N)
  H_random_pairwise_gene_euclidean<- rep(NA, N)

  KO_pairwise_gene_pearson<- rep(NA, N)
  H_KO_pairwise_gene_pearson<- rep(NA, N)

  KO_pairwise_gene_euclidean<- rep(NA, N)
  H_KO_pairwise_gene_euclidean<- rep(NA, N)

  dim_matrix<- length(table(RNAseq_Annotation_Matrix_no_sd_of_zero$Bin))
  All_KOs<- names(which(table(RNAseq_Annotation_Matrix_no_sd_of_zero$KO) > 5)) [-1] #list of all KOs which appear greater than 5 times
  Pairwise_Bin_Array_Presence<- matrix(0, dim_matrix, length(All_KOs))
  rownames(Pairwise_Bin_Array_Presence)<- names(table(RNAseq_Annotation_Matrix_no_sd_of_zero$Bin))[
                                          order(as.numeric(names(table(RNAseq_Annotation_Matrix_no_sd_of_zero$Bin))))]


  RNAseqExpresssionCounts <- as.matrix(RNAseq_Annotation_Matrix_no_sd_of_zero[, matrix_features@SS:matrix_features@SE])
  RNAseqExpressionRanks <- as.matrix(RNAseq_Annotation_Matrix_no_sd_of_zero[, matrix_features@RS:matrix_features@RE])

  for (x in 1:N) {
    random_genomes<- sample(length(matrix_features@high_quality_bins), 2) # grab 2 genomes, no replacing
    position_of_genome_A<- which(RNAseq_Annotation_Matrix_no_sd_of_zero$Bin == matrix_features@high_quality_bins[random_genomes[1]])
    position_of_genome_B<- which(RNAseq_Annotation_Matrix_no_sd_of_zero$Bin == matrix_features@high_quality_bins[random_genomes[2]])
    position_of_A<- sample(position_of_genome_A, 1)
    position_of_B<- sample(position_of_genome_B, 1)

    # m0
    # Calculate pearson correlation and euclidean distance for random genes in random genomes.
    # The results per iteration are saved in a vector
    random_pairwise_gene_pearson[x]<- cor(as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_A,
                                                                                            matrix_features@SS : matrix_features@SE]),
                                          as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_B,
                                                                                            matrix_features@SS : matrix_features@SE]))

    random_pairwise_gene_euclidean[x]<- Calc_Norm_Euc(as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_A,
                                                                                            matrix_features@RS : matrix_features@RE]),
                                                      as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[position_of_B,
                                                                                            matrix_features@RS : matrix_features@RE]))

    # m1
    # Same as m0, except now the comparison is done for random genes with the same KO term.
    # Sample_KO_Position_of_A is a single value from all possible positions shared between genome and KO match.
    sample_KO_positions <- sample_KO(RNAseq_Annotation_Matrix_no_sd_of_zero, position_of_genome_A, position_of_genome_B)
    KO_pairwise_gene_pearson[x]<- cor(as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[sample_KO_positions$sample_KO_position_of_A,
                                                                                            matrix_features@SS : matrix_features@SE]),
                                      as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[sample_KO_positions$sample_KO_position_of_B,
                                                                                            matrix_features@SS : matrix_features@SE]))

    KO_pairwise_gene_euclidean[x]<- Calc_Norm_Euc(as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[sample_KO_positions$sample_KO_position_of_A,
                                                                                            matrix_features@RS : matrix_features@RE]),
                                                  as.numeric(RNAseq_Annotation_Matrix_no_sd_of_zero[sample_KO_positions$sample_KO_position_of_B,
                                                                                            matrix_features@RS : matrix_features@RE]))
    # m1.1
    Array_Pearson_Euclidean<-NULL
    if(language == 'R'){
      Array_Pearson_Euclidean <-   Cor_Matrix(sample_KO_positions$KO_positions_of_A,
                                              sample_KO_positions$KO_positions_of_B,
                                              RNAseq_Annotation_Matrix_no_sd_of_zero,
                                              matrix_features)
    } else if(language == 'C'){
      Array_Pearson_Euclidean <- Cor_Matrix_C(sample_KO_positions$KO_positions_of_A,
                                              sample_KO_positions$KO_positions_of_B,
                                              RNAseqExpresssionCounts,
                                              RNAseqExpressionRanks
                                              )
    }
    H_KO_pairwise_gene_pearson[x]<- Array_Pearson_Euclidean$correlation[which.min((1-Array_Pearson_Euclidean$correlation) +
                                                                             Array_Pearson_Euclidean$euclidean)] # max(Array_Pearson_Euclidean[random_row,,1],na.rm=TRUE)
    H_KO_pairwise_gene_euclidean[x]<- Array_Pearson_Euclidean$euclidean[which.min((1-Array_Pearson_Euclidean$correlation)+
                                                                               Array_Pearson_Euclidean$euclidean)] # min(Array_Pearson_Euclidean[random_row,,2],na.rm=TRUE)

    # m0.1
    position_of_As<- sample(position_of_genome_A,
                            length(sample_KO_positions$KO_positions_of_A))
    position_of_Bs<- sample(position_of_genome_B,
                            length(sample_KO_positions$KO_positions_of_B))

    Array_Pearson_Euclidean_random<-NULL
    if(language == 'R'){
      Array_Pearson_Euclidean_random<- Cor_Matrix(position_of_As,
                                               position_of_Bs,
                                               RNAseq_Annotation_Matrix_no_sd_of_zero,
                                               matrix_features)
    } else if(language == 'C'){
      Array_Pearson_Euclidean_random <- Cor_Matrix_C(position_of_As,
                                                     position_of_Bs,
                                                     RNAseqExpresssionCounts,
                                                     RNAseqExpressionRanks
      )
    }

    random_position_matrix <- sample(Array_Pearson_Euclidean_random$correlation, 1)
    random_pairwise_gene_pearson[x]<- Array_Pearson_Euclidean_random$correlation[1]
    random_pairwise_gene_euclidean[x]<- Array_Pearson_Euclidean_random$euclidean[1]

    H_random_pairwise_gene_pearson[x]<- Array_Pearson_Euclidean_random$correlation[which.min((1-Array_Pearson_Euclidean_random$correlation)+
                                                                                        Array_Pearson_Euclidean_random$euclidean)]
    H_random_pairwise_gene_euclidean[x]<- Array_Pearson_Euclidean_random$euclidean[which.min((1-Array_Pearson_Euclidean_random$correlation)+
                                                                                          Array_Pearson_Euclidean_random$euclidean)]
  }

  newList<- list("random_pairwise_gene_pearson" = random_pairwise_gene_pearson,
                  "H_random_pairwise_gene_pearson" = H_random_pairwise_gene_pearson,
                  "KO_pairwise_gene_pearson" = KO_pairwise_gene_pearson,
                  "H_KO_pairwise_gene_pearson"= H_KO_pairwise_gene_pearson,
                  "random_pairwise_gene_euclidean" = random_pairwise_gene_euclidean,
                  "H_random_pairwise_gene_euclidean" = H_random_pairwise_gene_euclidean,
                  "KO_pairwise_gene_euclidean" = KO_pairwise_gene_euclidean,
                  "H_KO_pairwise_gene_euclidean"=H_KO_pairwise_gene_euclidean)

  return(newList)
}

sample_KO <- function(RNAseq_Annotation_Matrix_no_sd_of_zero, position_of_genome_A, position_of_genome_B){

  shared_KO<- intersect(RNAseq_Annotation_Matrix_no_sd_of_zero$KO[position_of_genome_A],
                        RNAseq_Annotation_Matrix_no_sd_of_zero$KO[position_of_genome_B])

  shared_KO<- shared_KO[shared_KO != ""]
  random_KO<- shared_KO[sample(length(shared_KO), 1)]

  KO_positions_of_A<- intersect(which(RNAseq_Annotation_Matrix_no_sd_of_zero$KO == random_KO), position_of_genome_A)
  KO_positions_of_B<- intersect(which(RNAseq_Annotation_Matrix_no_sd_of_zero$KO == random_KO), position_of_genome_B)


  sample_KO_position_of_A<- KO_positions_of_A[sample(length(KO_positions_of_A), 1)]
  sample_KO_position_of_B<- KO_positions_of_B[sample(length(KO_positions_of_B), 1)]
  return(list("sample_KO_position_of_A" = sample_KO_position_of_A,
              "sample_KO_position_of_B" = sample_KO_position_of_B,
              "KO_positions_of_A" = KO_positions_of_A,
              "KO_positions_of_B" = KO_positions_of_B))
}
