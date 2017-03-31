#' Calculte the Jaccard Distance
#'
#' This function calculates the Jaccard Distances between all bins for a given module
#'
#' @param Module
#' @export
#' @return A matrix of pairwise Jaccard distnces
#' @examples 	Jaccard_Distance_Function(PHA_Module)

Jaccard_Distance_Function <- function(Module) {
  Jaccard_Module_Distance <- matrix(NA,nrow=length(high_quality_bins),ncol=length(high_quality_bins))
  colnames(Jaccard_Module_Distance)<-names(table(RNAseq_Annotated_Matrix$Bin))[order(as.numeric(names(table(RNAseq_Annotated_Matrix$Bin))))]
  rownames(Jaccard_Module_Distance)<-colnames(Pairwise_Bin_Array_Presence)
  for (x in 1:37) { #iterate over lower part of matrix
    for (y in (x+1):38) {
      Jaccard_Module_Distance[x,y]<-Calc_Jaccard(Pairwise_Bin_Array_Presence[x, which(All_KOs%in%Module)],Pairwise_Bin_Array_Presence[y, which(All_KOs%in%Module)])
      # Module<-return(Module)
    }
  }
  return(Jaccard_Module_Distance)
}
