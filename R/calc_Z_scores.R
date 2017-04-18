#' Calcultes mu, sd and Z scores
#'
#' This function calculates the Jaccard Distances between all bins for a given module
#'
#' @param I_KOs_Background a background distribution for KOs
#' @export
#' @return a list of mu, sd and Z scores of the background distribution
#' @examples 	Z_scores1<-calc_Z_scores(I_KOs_Background)

calc_Z_scores <- function (I_KOs_Background) {
  mu_values<-sapply(I_KOs_Background,mean)
  sd_values<-sapply(I_KOs_Background,sd)
  
  Z_H_random_pairwise_gene_correlation<-((I_KOs_Background[[2]])-mu_values[2])/sd_values[2]
  Z_H_random_pairwise_gene_euclidean<-((I_KOs_Background[[6]])-mu_values[6])/sd_values[6]
  
  Z_H_KO_pairwise_gene_correlation<-((I_KOs_Background[[4]])-mu_values[2])/sd_values[2]
  Z_H_KO_pairwise_gene_euclidean<-((I_KOs_Background[[8]])-mu_values[6])/sd_values[6]
  
  # df_random<-as.data.frame(cbind((Z_random_pairwise_gene_correlation),(Z_random_pairwise_gene_euclidean)))
  # df_KO<-as.data.frame(cbind((KO_pairwise_gene_correlation),(KO_pairwise_gene_euclidean)))
  # df_H_KO<-as.data.frame(cbind((Z_H_KO_pairwise_gene_euclidean),(H_KO_pairwise_gene_euclidean)))
  
  Z_df_random<-as.data.frame(cbind(-(Z_H_random_pairwise_gene_correlation),(Z_H_random_pairwise_gene_euclidean)))
  Z_df_H_KO<-as.data.frame(cbind(-(Z_H_KO_pairwise_gene_correlation),(Z_H_KO_pairwise_gene_euclidean)))
  
  newList<-list("mu"=mu_values,"sd"=sd_values,"Z_H_random"=Z_df_random,"Z_H_KO"=Z_df_H_KO)
  return(newList)
}