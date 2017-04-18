#' Plot background distribution
#'
#' A 4-panel plot
#'
#' @param I_KOs_Background, the background distributions (output from Individual_KOs_Background function)
#' @export
#' @return TA 4-panel plot
#' @examples plot_individual_bd(I_KOs_Background_d_BR)

plot_individual_bd <- function(I_KOs_Background) {
  t.test_KO_random_pearson<-t.test(I_KOs_Background$KO_pairwise_gene_correlation,I_KOs_Background$random_pairwise_gene_correlation, alternative="greater") # x > y (NULL)
  t.test_H_KO_H_random_pearson<-t.test(I_KOs_Background$H_KO_pairwise_gene_correlation,I_KOs_Background$H_random_pairwise_gene_correlation, alternative="greater")
  t.test_H_KO_KO_pearson<-t.test(I_KOs_Background$H_KO_pairwise_gene_correlation,I_KOs_Background$KO_pairwise_gene_correlation, alternative="greater")
  
  t.test_KO_random_euclidean<-t.test(I_KOs_Background$KO_pairwise_gene_euclidean,I_KOs_Background$random_pairwise_gene_euclidean, alternative="less") # x > y (NULL)
  t.test_H_KO_H_random_euclidean<-t.test(I_KOs_Background$H_KO_pairwise_gene_euclidean,I_KOs_Background$H_random_pairwise_gene_euclidean, alternative="less")
  t.test_H_KO_KO_euclidean<-t.test(I_KOs_Background$H_KO_pairwise_gene_euclidean,I_KOs_Background$KO_pairwise_gene_euclidean, alternative="less")
  
  # t.test_KO_random_pearson$p.value
  # t.test_H_KO_H_random_pearson$p.value
  # t.test_H_KO_KO_pearson # This is for comparing the two KO distributions
  # t.test_KO_random_euclidean$p.value
  # t.test_H_KO_H_random_euclidean$p.value
  # t.test_H_KO_KO_euclidean # This is for comparing the two KO distributions
  # if necessary, add a line to return the resuts from the t.tests
  
  par(mfrow=c(2,2),mar=c(3,3,3,1))
  # plot 1
  plot(density(I_KOs_Background_d$random_pairwise_gene_correlation,adjust = 2,na.rm=TRUE),ylim=c(0,1),xlab="",ylab="",main="")
  points(density(I_KOs_Background_d$KO_pairwise_gene_correlation,adjust = 2),typ="l",col="blue")
  mtext(paste("p-value = ",signif(t.test_KO_random_pearson$p.value,2)),side=3,col="blue",padj=2,cex=.75)
  title(ylab="Density", line=2, cex.lab=1)
  title(xlab="PC", line=2, cex.lab=1)
  
  # plot 2
  plot(density(I_KOs_Background$H_random_pairwise_gene_correlation,adjust = 2),ylim=c(0,1),xlab="",ylab="",main=" ")
  points(density(I_KOs_Background$H_KO_pairwise_gene_correlation,adjust = 2),typ="l",col="red")
  mtext(paste("p-value = ",signif(t.test_H_KO_H_random_pearson$p.value,2)),side=3,col="red",padj=2,cex=.75)
  title(ylab="Density", line=2, cex.lab=1)
  title(xlab="PC", line=2, cex.lab=1)
  
  # plot 3
  plot(density(I_KOs_Background$random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
  points(density(I_KOs_Background$KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="blue")
  title(ylab="Density", line=2, cex.lab=1)
  title(xlab="NRED", line=2, cex.lab=1)
  mtext(paste("p-value = ",signif(t.test_KO_random_euclidean$p.value,2)),side=3,col="blue",padj=2,cex=.75)
  
  # plot 4
  plot(density(I_KOs_Background$H_random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
  points(density(I_KOs_Background$H_KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="red")
  title(ylab="Density", line=2, cex.lab=1)
  title(xlab="NRED", line=2, cex.lab=1)
  title(" \n\nComparison of random & functional \n pairwise comparisons", outer=TRUE) 
  mtext(paste("p-value = ",signif(t.test_H_KO_H_random_euclidean$p.value,2)),side=3,col="red",padj=2,cex=.75)
  
}