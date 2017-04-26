#' Calculte the Jaccard Distance
#'
#' This function calculates the Jaccard Distances between all bins for a given module
#'
#' @param Module
#' @export
#' @return A matrix of pairwise Jaccard distnces
#' @examples 	J_PHA <-Jaccard_Distance_Function(RNAseq_Annotated_Matrix_BR_default_bin, matrix_features_BR, PHA_module)

Jaccard_Distance_Function <- function(RNAseq_Annotated_Matrix, matrix_features, Module) {
#  Pairwise_Bin_Array_Presence	<- Presence_Absence_Matrix(RNAseq_Annotated_Matrix)
  Jaccard_Module_Distance <- matrix(data = NA, 
                                    nrow = length(high_quality_bins),
                                    ncol = length(high_quality_bins),
                                    dimnames = list(sort(matrix_features_BR@high_quality_bins),
                                                    sort(matrix_features_BR@high_quality_bins)))
  
  for (x in 1:37) { #iterate over lower part of matrix
    for (y in (x+1):38) {
      Jaccard_Module_Distance[x,y]<-Calc_Jaccard(matrix_features@Pairwise_Bin_Array_Presence[x, 
                                                                             which(matrix_features@All_KOs%in%Module)],
                                                 matrix_features@Pairwise_Bin_Array_Presence[y, 
                                                                            which(matrix_features@All_KOs%in%Module)])
    }
  }
  return(Jaccard_Module_Distance)
}
