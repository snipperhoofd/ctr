#' Generate a random module of size 'N'
#'
#' This function generates a vector of length 'N' composed of randomly selected annotated features
#'
#' @param All_KOs  A vector of all annotated features
#' @param N The number of elements to be included in the randomly generated module
#'
#' @return a vector containing N random annoated features (without replacement).
#' @examples Generate_Random_Module(All_KOs,6)

Generate_Random_Module <- function(All_KOs,N) {
  Random_Module<-sample(All_KOs,N,replace=FALSE)
  return(Random_Module)
}