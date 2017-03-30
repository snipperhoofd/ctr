#' Calculate Jaccard Distance
#'
#' Calculates the Jaccard Distance between two vectors
#' \deqn{J(A,B) = {{|A \cap B|}\over{|A \cup B|}} = {{|A \cap B|}\over{|A| + |B| - |A \cap B|}}}
#'
#' @param vector1,vector2  Two vectors from which the Jaccard disance is to be calculated
#'
#' @return The Jaccard Distance between two two vectors.
#' @examples Jaccard<-Calc_Jaccard(vector1,vector2)


Calc_Jaccard <- function(vector1, vector2){
  Jaccard_Distance<-1-(sum(which(vector1==1)%in%which(vector2==1)) / (sum(which(vector1==0)%in%which(vector2==1)) + sum(which(vector1==1)%in%which(vector2==0)) + sum(which(vector1==1)%in%which(vector2==1))))
  return(Jaccard_Distance)
}
