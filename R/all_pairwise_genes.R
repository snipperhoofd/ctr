#' Calculating all pairwise gene distances
#'
#' This function calculates all pairwise distances across each bin and gene (KO) combination
#' It is common to have numerous genes with the same annotation in a genome. When this occurs
#' the highest scoring pair between two bins is used.
#'
#'
#' @param RNAseq_Annotation_Matrix, the original matrix
#' @param Z_scores_BR, Z-scores obtained from the individual background distribution 
#' @param matrix_features, the matrix features
#' @export
#' @return An array containing the pairwise distances between bins for all KOs and their positions in the original RNAseq_Annotation_Matrix
#' @examples 	all_pairwise_genes(RNAseq_Annotation_Matrix_BR,Z_scores_BR,matrix_features_BR)
# first remove rows with standard deviations of 0

all_pairwise_genes <- function(RNAseq_Annotation_Matrix, Z_scores, matrix_features){

  # build empty arrays for each Pearson Correlation and NRED
  no_annotation <- which(names((table(RNAseq_Annotation_Matrix$KO)))=="")
  All_KOs<- names((table(RNAseq_Annotation_Matrix$KO)))[-no_annotation]
  
  dim_matrix<-length(matrix_features@high_quality_bins)
  H_KO_pairwise_gene_pearson<- array(NA, dim=c(dim_matrix,
                                             dim_matrix,
                                             length(All_KOs)))
  colnames(H_KO_pairwise_gene_pearson)<- names(table(RNAseq_Annotated_Matrix$Bin))[order(as.numeric(names(table(RNAseq_Annotated_Matrix$Bin))))]
  rownames(H_KO_pairwise_gene_pearson)<- colnames(Pairwise_Bin_Array_Pearson)
  
  H_KO_pairwise_gene_euclidean<- H_KO_pairwise_gene_pearson
  Pairwise_PositionsA<- H_KO_pairwise_gene_pearson
  Pairwise_PositionsB<- H_KO_pairwise_gene_pearson
  
  ######### Calculate the minimum pairwise distance between genome bins for all KOs, converted to Z score, and keeping the pair with the highest sum z-score
  
  for (x in 1:(dim(H_KO_pairwise_gene_pearson)[1]-1)) {
    for (y in (x+1):dim(H_KO_pairwise_gene_pearson)[2]) {
      for (z in 1:(dim(H_KO_pairwise_gene_pearson)[3])) { #iterate over array
        
        # Identify the rows in the original matrix corresponding to each genome
        position_of_genome_A = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_Bin_Array_Presence)[x])
        position_of_genome_B = which(RNAseq_Annotated_Matrix$Bin==rownames(Pairwise_Bin_Array_Presence)[y])
        # Second identify the rows in the original matrix corresponding to a KO
        position_of_kegg_enzyme_A = intersect(which(RNAseq_Annotated_Matrix$KO==All_KOs[z]), position_of_genome_A)
        position_of_kegg_enzyme_B = intersect(which(RNAseq_Annotated_Matrix$KO==All_KOs[z]), position_of_genome_B)
        # Make sure the KO is present in both genomes
        if (!length(position_of_kegg_enzyme_A)==0 && !length(position_of_kegg_enzyme_B)==0) {
          # Conduct all pairwise comparisons between Pearson Correlations and Normalized Euclidean Distances, converting to Z scores
          # First define two empty matrices and then fill them with the PCC and Euc distances
          max_pairwise_gene_correlation<- matrix(NA, 
                                                nrow=length(position_of_kegg_enzyme_A), 
                                                ncol=length(position_of_kegg_enzyme_B))
          max_pairwise_gene_euclidean<- max_pairwise_gene_correlation
          
          for (m in 1:length(position_of_kegg_enzyme_A)){
            for (n in 1:length(position_of_kegg_enzyme_B)){
              # make sure there is always a standard deviation, or else cor gives an error. If there is a sd, proceed with calculations
              if (sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m],matrix_features@SS:matrix_features@SE]))!=0 && 
                  sd(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n],matrix_features@SS:matrix_features@SE]))!=0) {
                
                max_pairwise_gene_correlation[m,n]<- (cor(as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@SS:matrix_features@SE]),
                                                          as.numeric(RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@SS:matrix_features@SE])))
                
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@RS:matrix_features@RE] - 
                                   RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@RS:matrix_features@RE]
                
                max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
              } else {
                # If there is no standard deviation, the correlation is NA
                max_pairwise_gene_correlation[m,n]<- NA
                subtracted_lists<- RNAseq_Annotated_Matrix[position_of_kegg_enzyme_A[m], matrix_features@RS:matrix_features@RE] - 
                                   RNAseq_Annotated_Matrix[position_of_kegg_enzyme_B[n], matrix_features@RS:matrix_features@RE]
                
                max_pairwise_gene_euclidean[m,n]<-sqrt(sum(subtracted_lists * subtracted_lists))
              }
            }
          }
          # Convert to Z scores
          Zscore_pairwise_gene_correlation<-((max_pairwise_gene_correlation-Z_scores$mu[2])/Z_scores$sd[2]) # need to inverse PCC
          Zscore_pairwise_gene_euclidean<-((max_pairwise_gene_euclidean-Z_scores$mu[6])/Z_scores$mu[6])
          
          best_scoring_pair<-which.min((1-max_pairwise_gene_correlation)+(Zscore_pairwise_gene_euclidean))
          rownames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_A
          colnames(max_pairwise_gene_correlation)<-position_of_kegg_enzyme_B
          
          if (length(best_scoring_pair)>0) {
            H_KO_pairwise_gene_pearson[x,y,z]<-max_pairwise_gene_correlation[best_scoring_pair]
            H_KO_pairwise_gene_euclidean[x,y,z]<-max_pairwise_gene_euclidean[best_scoring_pair]
            Pairwise_PositionsA[x,y,z]<-rownames(max_pairwise_gene_correlation)[best_scoring_pair]
            Pairwise_PositionsB[x,y,z]<-colnames(max_pairwise_gene_correlation)[best_scoring_pair]
          } else {
            scoring_pair<-which.min(max_pairwise_gene_euclidean)
            H_KO_pairwise_gene_pearson[x,y,z]<-max_pairwise_gene_correlation[scoring_pair]
            H_KO_pairwise_gene_euclidean[x,y,z]<-max_pairwise_gene_euclidean[scoring_pair]
            Pairwise_PositionsA[x,y,z]<-rownames(max_pairwise_gene_correlation)[scoring_pair]
            Pairwise_PositionsB[x,y,z]<-colnames(max_pairwise_gene_correlation)[scoring_pair]}
        } else {next}
      #  print(c(x,y,z))
      }
    }
  }
  
  Zscore_pearson<- ((H_KO_pairwise_gene_pearson - Z_scores$mu[2]) / Z_scores$sd[2]) # need to inverse PCC
  Zscore_nred<- ((H_KO_pairwise_gene_euclidean - Z_scores$sd[6]) / Z_scores$sd[6])
  Composite_Z_Score<- ((-Zscore_pearson) + Zscore_nred)
  
  newList <- list("pearsons" = H_KO_pairwise_gene_pearson, 
                  "nred" = H_KO_pairwise_gene_euclidean, 
                  "positionsA" = Pairwise_PositionsA, 
                  "positionsB" = Pairwise_PositionsB,
                  "Zscore_pearson" = Zscore_pearson,
                  "Zscore_nred" = Zscore_nred,
                  "Composite_Z_Score" = Composite_Z_Score)
  
  return(newList)
}
