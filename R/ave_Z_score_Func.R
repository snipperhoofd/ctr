#' Calculates the average Z Score
#'
#' This function converts an array of Z-scores into a matrix of average scores.
#' (may be used on outputs for F11-F13)
#'
#' @export
#' @param combined_scores_array An array of all the pairwise Z-scores in a module
#' @return a pairwise correlation matrix of average Z-scores
#' @examples ave_Z_score_matrix <-ave_Z_score_Func(PHA_module_P_NRED$combined)

ave_Z_score_Func <- function(combined_scores_array){
  ave_Z_score<-matrix(NA,
                      ncol=length(high_quality_bins),
                      nrow=length(high_quality_bins),
                      dimnames = list(sort(matrix_features_BR@high_quality_bins),
                                      sort(matrix_features_BR@high_quality_bins)))
  for (x in 1:(dim(combined_scores_array)[1]-1)) {
    for (y in (x+1):dim(combined_scores_array)[2]) {
      ave_Z_score[x,y]<-mean(combined_scores_array[x,y,],na.rm=TRUE)}}
  return(ave_Z_score)
}
