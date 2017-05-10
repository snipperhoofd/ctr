CTR <- setRefClass("CongruentTranscriptionalResponse",
                   fields = list(
                     RNAseq_Annotated_Matrix = "data.frame",
                     matrix_features = "General_features",
                     I_KOs_Background = "vector",
                     bg_distance_modules = "vector",
                     Z_scores = "list"
                    ),
                    methods = list(
                      Run = function(show_plots = TRUE){
                        print("Steps indicated with an asterisk(*) are relatively slow")

                        print("__Normalizing RNAseq data__")
                        RNAseq_Annotated_Matrix <<- RNAseq_Normalize(RNAseq_Annotated_Matrix,
                                                                     matrix_features,
                                                                     method = "default")
                        RNAseq_Annotated_Matrix <<- Normalize_by_bin(RNAseq_Annotated_Matrix,
                                                                     matrix_features)


                        print("__Creating rank columns__")
                        set.seed(1396558)
                        RNAseq_Annotated_Matrix <<- Create_Rank_Columns(RNAseq_Annotated_Matrix,
                                                                        matrix_features)


                        print("__Removing rows with a stdev of 0__")
                        RNAseq_Annotated_Matrix <<- which_rows_with_no_sd(RNAseq_Annotated_Matrix,
                                                                          matrix_features)


                        print("__*Calculating background distribution for random individual KOs*__")
                        I_KOs_Background <<- Individual_KOs_Background(RNAseq_Annotated_Matrix,
                                                                       matrix_features,10000)

                        print("Calculating Z scores of individual KO background distributions")
                        Z_scores <<- calc_Z_scores(I_KOs_Background)

                        print("__*Calculating Backround Distributions for random Modules of sizes 7*__")
                        #bg_distance_modules_6 <- Background_Distribution_Modules(RNAseq_Annotated_Matrix,
                         #                                                        matrix_features,
                          #                                                       Z_scores,
                           #                                                      6, 10000, 4)
                        bg_distance_modules_7 <- Background_Distribution_Modules(RNAseq_Annotated_Matrix,
                                                                                 matrix_features,
                                                                                 Z_scores,
                                                                                 7, 10000, 4)
                        #bg_distance_modules_8 <- Background_Distribution_Modules(RNAseq_Annotated_Matrix,
                         #                                                        matrix_features,
                          #                                                       Z_scores,
                           #                                                      8, 10000, 4)



                        bg_distance_modules <<- bg_distance_modules_7
                      },
                      AssociationMatrix_perModule = function(module, module_list){
                      pairwise_KO_distances <- P_NRED_Distance_Function(RNAseq_Annotated_Matrix,
                                                                        Z_scores,
                                                                        matrix_features,
                                                                        module)

                      clustering_results_P_NRED <- cluster_func(RNAseq_Annotated_Matrix_default_bin,
                                                                AA_pairwise_KO_distances$combined,
                                                                matrix_features,
                                                                module_list)

                      association_matrix <- fill_association_matrix(clustering_results_P_NRED,
                                                                    matrix_features,
                                                                    names(module_list))

                                        },
                      plotIndividualBackgroundDist = function(){
                        rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))



                        plot(hexbin(Z_scores$Z_H_random), colramp=rf, mincnt=1, maxcnt=max(hexbin(Z_scores$Z_H_KO)@count),ylab="NRED",xlab="PC")
                        plot(hexbin(Z_scores$Z_H_KO), colramp=rf, mincnt=1, maxcnt=max(hexbin(Z_scores$Z_H_KO)@count),ylab="NRED",xlab="PC")


                        t_test_KO_random_pearson<-t.test(I_KOs_Background$KO_pairwise_gene_pearson,I_KOs_Background$random_pairwise_gene_pearson, alternative="greater") # x > y (NULL)
                        t_test_KO_random_euclidean<-t.test(I_KOs_Background$random_pairwise_gene_euclidean,I_KOs_Background$H_KO_pairwise_gene_euclidean, alternative="greater") # x > y (NULL)

                        par(mfrow=c(2,2),mar=c(3,3,3,1))
                        # plot 1
                        plot(density(I_KOs_Background$random_pairwise_gene_pearson,adjust = 2,na.rm=TRUE),ylim=c(0,1),xlab="",ylab="",main="")
                        points(density(I_KOs_Background$KO_pairwise_gene_pearson,adjust = 2),typ="l",col="blue")
                        mtext(paste("p-value = ",signif(t_test_KO_random_pearson$p.value,2)),side=3,col="blue",padj=2,cex=.75)
                        title(ylab="Density", line=2, cex.lab=1)
                        title(xlab="PC", line=2, cex.lab=1)

                        # plot 2
                        plot(density(I_KOs_Background$H_random_pairwise_gene_pearson,adjust = 2),ylim=c(0,1),xlab="",ylab="",main=" ")
                        points(density(I_KOs_Background$H_KO_pairwise_gene_pearson,adjust = 2),typ="l",col="red")
                        mtext(paste("p-value = ",signif(t_test_KO_random_pearson$p.value,2)),side=3,col="red",padj=2,cex=.75)
                        title(ylab="Density", line=2, cex.lab=1)
                        title(xlab="PC", line=2, cex.lab=1)

                        # plot 3
                        plot(density(I_KOs_Background$random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
                        points(density(I_KOs_Background$KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="blue")
                        title(ylab="Density", line=2, cex.lab=1)
                        title(xlab="NRED", line=2, cex.lab=1)
                        mtext(paste("p-value = ",signif(t_test_KO_random_euclidean$p.value,2)),side=3,col="blue",padj=2,cex=.75)

                        # plot 4
                        plot(density(I_KOs_Background$H_random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
                        points(density(I_KOs_Background$H_KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="red")
                        title(ylab="Density", line=2, cex.lab=1)
                        title(xlab="NRED", line=2, cex.lab=1)
                        title(" \n\nComparison of random & functional \n pairwise comparisons", outer=TRUE)
                        mtext(paste("p-value = ",signif(t_test_KO_random_euclidean$p.value,2)),side=3,col="red",padj=2,cex=.75)
                      },
                      plotModuleBackgroundDist = function(){
                        plot(density(bg_distance_modules,na.rm=TRUE),col="red",main="Background distributions of modules \nof size \"N\"\nBR",ylim=c(0,1),xlim=c(-4,4),cex.main=.75)
                        points(density(bg_distance_modules,na.rm=TRUE),main="N=7",col="blue",type="l")
                        points(density(bg_distance_modules,na.rm=TRUE),main="N=8",col="green",type="l")
                        legend("topright",legend=c("N=6","N=6","N=6"),col=c("red","blue","green"),lty=c(1,1,1))
                      }
                    )
                   )


