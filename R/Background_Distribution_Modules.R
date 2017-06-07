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
#' @param Z The number of iterations used to calculate  background distribution


#' @export
#' @return a list of vectors containing XXX & YYY
#' @examples Background_Distribution_Modules(RNAseq_Annotated_Matrix,6,1000)

Background_Distribution_Modules <- function(RNAseq_Annotated_Matrix, matrix_features,
                                            Z_scores, N, Z, P = 2) {
  library(doParallel)



  Random_Jaccard_Distances <- rep(NA, Z)
  Random_Composite_Distances <- rep(NA, Z)
  All_KOs<-names(which(table(RNAseq_Annotated_Matrix$KO) >= 2))[-1] # This was originally a global variable but was moved so that it can change depending on the annotation matrix used

  #For paralellization
  cl <-makeCluster(P)
  registerDoParallel(cl)
  #Exporting existing functions to be available in the "worker" nodes
  functionNames <- c("GetFeatures","RandomDistances", "comparePairwise",
                     "matrix_features")
  clusterExport(cl, varlist = functionNames, envir = environment())


  # iterate Z times
  RandomDistList <- foreach(i = 1:Z, .combine = 'comb', .multicombine = TRUE,
                           .init = list(list(), list())) %dopar% {
    library(ctr)

    #Initializing empty vectors
    Random_Zscore_Pearson_Distances <- rep(NA,N)
    Random_Zscore_Euclidean_Distances <- rep(NA,N)

    #Pick two random genomes
    random_genomes <- sample(length(matrix_features@high_quality_bins), 2)

    KO_A<-unique(RNAseq_Annotated_Matrix$KO[which(RNAseq_Annotated_Matrix$Bin == random_genomes[1])])
    KO_B<-unique(RNAseq_Annotated_Matrix$KO[which(RNAseq_Annotated_Matrix$Bin == random_genomes[2])])
    intersection_AB <- intersect(KO_A,KO_B)
    random_module <- Generate_Random_Module(intersection_AB, N)
    All_position_KOs <- which(colnames(matrix_features@Pairwise_Bin_Array_Presence) %in% random_module)

    # Calculate Jaccard Distance
    PA_position_of_genome_A <- which(rownames(matrix_features@Pairwise_Bin_Array_Presence)==
                                     random_genomes[1])
    PA_position_of_genome_B <- which(rownames(matrix_features@Pairwise_Bin_Array_Presence)==
                                     random_genomes[2])

    # Next calculate Pearson and NRED
    for (j in 1:N) {
      Random_Pearson_Distances <- NA
      Random_Euclidean_Distances <- NA

      features <- GetFeatures(RNAseq_Annotated_Matrix, matrix_features,
                              random_genomes, random_module, j)

        # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances
      pairwiseDistances <- comparePairwise(features$position_of_kegg_enzyme_A,
                                             features$position_of_kegg_enzyme_B,
                                             RNAseq_Annotated_Matrix,
                                             matrix_features)

      dist <- RandomDistances(pairwiseDistances, Z_scores)

      Random_Zscore_Pearson_Distances[j] <- dist$Zscore_Pearson
      Random_Zscore_Euclidean_Distances[j] <- dist$Zscore_Euclidean
    }

    composite_exp <- mean((-Random_Zscore_Pearson_Distances) +
                          Random_Zscore_Euclidean_Distances,
                          na.rm = TRUE)[1]

    jaccard_exp <- Calc_Jaccard(matrix_features@Pairwise_Bin_Array_Presence[PA_position_of_genome_A,
                                                            All_position_KOs],
                                matrix_features@Pairwise_Bin_Array_Presence[PA_position_of_genome_B,
                                                            All_position_KOs])
    list(jaccard_exp, composite_exp)
  }

  Random_Jaccard_Distances <- unlist(RandomDistList[[1]])
  Random_Composite_Distances <- unlist(RandomDistList[[2]])

  Random_Background_Module_Distances <- Random_Composite_Distances *
                                        (1 - Random_Jaccard_Distances)

  stopImplicitCluster()
  stopCluster(cl)
  return(Random_Background_Module_Distances)
}



GetFeatures <- function(RNAseq_Annotated_Matrix, matrix_features,
                        random_genomes, random_module, j){
  # Identify the rows in the original matrix corresponding to each genome
  position_of_genome_A <- which(RNAseq_Annotated_Matrix$Bin == random_genomes[1])
  position_of_genome_B <- which(RNAseq_Annotated_Matrix$Bin == random_genomes[2])

  # Identify the rows in the original matrix corresponding to a KO
  KO_positions <- which(RNAseq_Annotated_Matrix$KO == random_module[j])

  # Find intersection between the genome lists and the KO list
  position_of_kegg_enzyme_A <- intersect(KO_positions, position_of_genome_A)
  position_of_kegg_enzyme_B <- intersect(KO_positions, position_of_genome_B)

  # Make sure the KO is present in both genomes
  l_position_of_kegg_enzyme_A <- length(position_of_kegg_enzyme_A)
  l_position_of_kegg_enzyme_B <- length(position_of_kegg_enzyme_B)
  return(list("position_of_genome_A" = position_of_genome_A,
              "position_of_genome_B" = position_of_genome_B,
              "position_of_kegg_enzyme_A" = position_of_kegg_enzyme_A,
              "position_of_kegg_enzyme_B" = position_of_kegg_enzyme_B,
              "l_position_of_kegg_enzyme_A" = l_position_of_kegg_enzyme_A,
              "l_position_of_kegg_enzyme_B" = l_position_of_kegg_enzyme_B))
}

RandomDistances <- function(pairwiseDistances, Z_scores){
  # Convert to Z scores
  Zscore_pairwise_gene_correlation <- (pairwiseDistances$pairwise_correlation-Z_scores$mu[2])/Z_scores$sd[2] # need to inverse PCC
  Zscore_pairwise_gene_euclidean <- (pairwiseDistances$pairwise_euclidean-Z_scores$mu[6])/Z_scores$sd[6]

  #If there is a tie, the first minimum is taken
  best_scoring_pair<- which.min((1-Zscore_pairwise_gene_correlation)+(Zscore_pairwise_gene_euclidean))


  Random_Pearson_Distances <- pairwiseDistances$pairwise_correlation[best_scoring_pair]
  Random_Euclidean_Distances <- pairwiseDistances$pairwise_euclidean[best_scoring_pair]


  return(list("Random_Pearson_Distances" = Random_Pearson_Distances,
              "Random_Euclidean_Distances" = Random_Euclidean_Distances,
              "Zscore_Pearson" = Zscore_pairwise_gene_correlation[best_scoring_pair],
              "Zscore_Euclidean"= Zscore_pairwise_gene_euclidean[best_scoring_pair]))
}


comb <- function(x, ...){
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}



comparePairwise <- function(position_of_kegg_enzyme_A, position_of_kegg_enzyme_B,
                            RNAseq_Annotated_Matrix, matrix_features){

    l_position_of_kegg_enzyme_A <- length(position_of_kegg_enzyme_A)
    l_position_of_kegg_enzyme_B <- length(position_of_kegg_enzyme_B)
    max_pairwise_gene_correlation <- matrix(NA, nrow = l_position_of_kegg_enzyme_A,
                                            ncol = l_position_of_kegg_enzyme_B)
    max_pairwise_gene_euclidean <-max_pairwise_gene_correlation
    for (m in 1:l_position_of_kegg_enzyme_A){
      for (n in 1:l_position_of_kegg_enzyme_B){

        # use the no_sd dataset so that calculating a Pearson correlation never gives an error
        correlation <- cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],
                                                              matrix_features@SS:matrix_features@SE]),
                           as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],
                                                              matrix_features@SS:matrix_features@SE]))

        max_pairwise_gene_correlation[m, n]<- correlation

        subtracted_lists <- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@RS:matrix_features@RE] -
          RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@RS:matrix_features@RE]
        max_pairwise_gene_euclidean[m, n] <- sqrt(sum(subtracted_lists * subtracted_lists))
      }
    }
    return(list("pairwise_correlation" = max_pairwise_gene_correlation,
                "pairwise_euclidean" = max_pairwise_gene_euclidean))
}






