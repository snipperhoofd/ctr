#' Generate a background distibution of composite Z scores
#'
#' This function generates a background distibution based on randomly generated modules of size 'N'.
#' For each element of the module, pairwise Pearson Correlation (PC), and the Normalized Rank
#' Euclidean Distance (NRED) are calculated between all genomes. When two genomes contain multiple
#' elements for comparison, this is handled either by takng a random pairwise comparison, 
#' or by taking the score with the maximum similarity. A background distribution using each method
#' is provided for comparison. (Currently it only provides the score that minimizes the distances)
#' 
#' @param N The number of elements to be included in the randomly generated module
#' @param Z The number of iterations used to calculatea background distribution
#'
#' @return a list of vectors containing XXX & YYY
#' @examples Random_Background_Module_Distances_6<-Background_Distribution_Modules(6,1000)

Background_Distribution_Modules <- function(N,Z) {
  
  Random_Jaccard_Distances<-rep(NA,Z)
  Random_Composite_Distances<-rep(NA,Z)
  Random_Pearson_Distances<-rep(NA,N)
  Random_Euclidean_Distances<-rep(NA,N)
  Random_Zscore_Pearson_Distances<-rep(NA,N)
  Random_Zscore_Euclidean_Distances<-rep(NA,N)
  dim_matrix<-length(table(RNAseq_Annotated_Matrix$Bin))
  
  for (i in 1:Z) {
    two_random_genomes<-sample(length(high_quality_bins),2)
    Random_Module<-Generate_Random_Module(All_KOs,N)
    
    # Calculate Jaccard Distance
    genome1<-which(rownames(Pairwise_Bin_Array_Presence)==high_quality_bins[two_random_genomes[1]])
    genome2<-which(rownames(Pairwise_Bin_Array_Presence)==high_quality_bins[two_random_genomes[2]])
    Random_Jaccard_Distances[i]<-Calc_Jaccard(Pairwise_Bin_Array_Presence[genome1, which(All_KOs%in%Random_Module)],Pairwise_Bin_Array_Presence[genome2, which(All_KOs%in%Random_Module)])
    
    # Next calculate Pearson and NRED
    for (j in 1:length(Random_Module)) {          
      # Identify the rows in the original matrix corresponding to each genome
      position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[two_random_genomes[1]])
      position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==high_quality_bins[two_random_genomes[2]])
      # Second identify the rows in the original matrix corresponding to a KO
      position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==Random_Module[j]),position_of_genome_A)
      position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==Random_Module[j]),position_of_genome_B)
      # Make sure the KO is present in both genomes
      if (!length(position_of_kegg_enzyme_A)==0 & !length(position_of_kegg_enzyme_B)==0) {
        # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores
        # First define two empty matrices and then fill them with the PCC and Euc distances
        max_pairwise_gene_correlation<-matrix(NA,nrow=length(position_of_kegg_enzyme_A),ncol=length(position_of_kegg_enzyme_B))
        max_pairwise_gene_euclidean<-max_pairwise_gene_correlation
        for (m in 1:length(position_of_kegg_enzyme_A)){
          for (n in 1:length(position_of_kegg_enzyme_B)){
            # make sure there is always a standard deviation, or else cor gives an error. If there is a sd, proceed with calculations
            if (sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],2:7]))!=0 & sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],2:7]))!=0) {
              max_pairwise_gene_correlation[m,n]<-(cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],2:7]),as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],2:7])))
              subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],10:15]-RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],10:15]
              max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
            } else {
              # If there is no standard deviation, the correlation is NA
              max_pairwise_gene_correlation[m,n]<-NA
              subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],10:15]-RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],10:15]
              max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
            }
          }
        }
        # Convert to Z scores							
        Zscore_pairwise_gene_correlation<-((max_pairwise_gene_correlation-mu_pearson)/sd_pearson) # need to inverse PCC
        Zscore_pairwise_gene_euclidean<-((max_pairwise_gene_euclidean-mu_euclidean)/sd_euclidean)
        best_scoring_pair<-which.min((1-Zscore_pairwise_gene_correlation)+(Zscore_pairwise_gene_euclidean))							
        if (length(best_scoring_pair)>0) {
          Random_Pearson_Distances<-max_pairwise_gene_correlation[best_scoring_pair]
          Random_Euclidean_Distances<-max_pairwise_gene_euclidean[best_scoring_pair]
        } else {
          scoring_pair<-which.min(Zscore_pairwise_gene_euclidean)
          Random_Pearson_Distances<-max_pairwise_gene_correlation[scoring_pair]
          Random_Euclidean_Distances<-max_pairwise_gene_euclidean[scoring_pair]
        }
      } else {
        Random_Pearson_Distances<-NA
        Random_Euclidean_Distances<-NA
      }
      Random_Zscore_Pearson_Distances[j]<-((Random_Pearson_Distances-mu_pearson)/sd_pearson) # need to inverse PCC
      Random_Zscore_Euclidean_Distances[j]<-((Random_Euclidean_Distances-mu_euclidean)/sd_euclidean)
      
    }
    Random_Composite_Distances[i]<-mean((-Random_Zscore_Pearson_Distances)+Random_Zscore_Euclidean_Distances,na.rm=TRUE)[1]
    Random_Background_Module_Distances<-Random_Composite_Distances*(1-Random_Jaccard_Distances)
  }
  return(Random_Background_Module_Distances)
}