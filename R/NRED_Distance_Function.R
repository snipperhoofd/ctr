#' Calculates pairwise NRED for a module
#'
#' This function calculates all pairwise NRED for a given module. When two genomes contain multiple
#' elements for comparison, this is handled either by takng a random pairwise comparison,
#' or by taking the score with the maximum similarity. A background distribution using each method
#' is provided for comparison. (Currently it only provides the score that minimizes the distances)
#'
#' @param Subset_KOs A list of KOs that form a module
#' @export
#' @return a list of vectors containing pairwise NRED scores (nred), their Z-scores (Zscore),
#' and the position in the matrix for the highest scoring pair between genomes A and B
#' (positionA & positionB respectively)
#' @examples PHA_module_NRED <- NRED_Distance_Function(PHA_module)

NRED_Distance_Function <- function(Subset_KOs) {

  # Define two congruent arrays to be filled during the second step. Name the columns and rows based on the genome bins
  dim_matrix<-length(table(RNAseq_Annotated_Matrix$Bin))
  Pairwise_NRED<-array(NA,c(dim_matrix,dim_matrix,length(Subset_KOs)))
  colnames(Pairwise_NRED)<-names(table(RNAseq_Annotated_Matrix$Bin))[order(as.numeric(names(table(RNAseq_Annotated_Matrix$Bin))))]
  rownames(Pairwise_NRED)<-colnames(Pairwise_NRED)
  Pairwise_PositionsA<-Pairwise_NRED
  Pairwise_PositionsB<-Pairwise_NRED

  ######### This is the maximum pairwise NRED between genome bins for all KOs, converted to Z score, and keeping the pair with the highest sum z-score

  for (x in 1:(dim(Pairwise_NRED)[1]-1)) {
    for (y in (x+1):dim(Pairwise_NRED)[2]) {
      for (z in 1:dim(Pairwise_NRED)[3]) { #iterate over array

        # Identify the rows in the original matrix corresponding to each genome
        position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_NRED)[x])
        position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_NRED)[y])
        # Second identify the rows in the original matrix corresponding to a KO
        position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]),position_of_genome_A)
        position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==Subset_KOs[z]),position_of_genome_B)
        # Make sure the KO is present in both genomes
        if (!length(position_of_kegg_enzyme_A)==0 & !length(position_of_kegg_enzyme_B)==0) {
          # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores
          # First define two empty matrices and then fill them with the PCC and Euc distances
          max_pairwise<-matrix(NA,nrow=length(position_of_kegg_enzyme_A),ncol=length(position_of_kegg_enzyme_B))
          for (m in 1:length(position_of_kegg_enzyme_A)){
            for (n in 1:length(position_of_kegg_enzyme_B)){
              # make sure there is always a standard deviation, or else cor gives an error. If there is a sd, proceed with calculations
              if (sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],2:7]))!=0 & sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],2:7]))!=0) {
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],10:15]-RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],10:15]
                max_pairwise[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
              } else {
                # If there is no standard deviation, the correlation is NA
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],10:15]-RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],10:15]
                max_pairwise[m,n]<-sqrt(sum(subtracted_lists* subtracted_lists))
              }
            }
          }
          # Convert to Z scores

          if (length(max_pairwise)>0) {
            Pairwise_NRED[x,y,z]<-min(max_pairwise)
            rownames(max_pairwise)<-position_of_kegg_enzyme_A
            colnames(max_pairwise)<-position_of_kegg_enzyme_B
            Pairwise_PositionsA[x,y,z]<-rownames(max_pairwise)[which(max_pairwise == min(max_pairwise), arr.ind = TRUE)[1]]
            Pairwise_PositionsB[x,y,z]<-colnames(max_pairwise)[which(max_pairwise == min(max_pairwise), arr.ind = TRUE)[2]]
          }else {next}
        }
      }
    }
  }

  Zscore_pairwise_gene_euclidean<-((Pairwise_NRED-mu_euclidean)/sd_euclidean)

  newList <- list("nred" = Zscore_pairwise_gene_euclidean,"Zscore"=Zscore_pairwise_gene_euclidean,"positionsA"=Pairwise_PositionsA,"positionsB"=Pairwise_PositionsB)

  return(newList)
}
