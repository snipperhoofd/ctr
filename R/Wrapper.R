#' @export CTR
#' @exportClass CTR
CTR <- setRefClass("CTR",
                   fields = list(
                     RNAseq_Annotated_Matrix = "data.frame",
                     matrix_features = "General_features",
                     I_KOs_Background = "vector",
                     bg_distance_modules = "list",
                     Z_scores = "list",
                     All_association_matrix = 'matrix'
                    ),
                    methods = list(
                      Run = function(iterations = 10000,
                                     random_module_sizes = c(6,7),
                                     parallel_cores = 2){
                        print("Steps indicated with an asterisk(*) are relatively slow")

                        print("Normalizing RNAseq data")
                        RNAseq_Annotated_Matrix <<- RNAseq_Normalize(RNAseq_Annotated_Matrix,
                                                                     matrix_features,
                                                                     method = "default")
                        RNAseq_Annotated_Matrix <<- Normalize_by_bin(RNAseq_Annotated_Matrix,
                                                                     matrix_features)


                        print("Creating rank columns")
                        set.seed(1396558)
                        RNAseq_Annotated_Matrix <<- Create_Rank_Columns(RNAseq_Annotated_Matrix,
                                                                        matrix_features)


                        print("Removing rows with a stdev of 0")
                        RNAseq_Annotated_Matrix <<- which_rows_with_no_sd(RNAseq_Annotated_Matrix,
                                                                          matrix_features)


                        print("*Calculating background distribution for random individual KOs*")
                        I_KOs_Background <<- Individual_KOs_Background(RNAseq_Annotated_Matrix,
                                                                       matrix_features,iterations)

                        print("Calculating Z scores of individual KO background distributions")
                        Z_scores <<- calc_Z_scores(I_KOs_Background)


                        bg_distance_modules <<- list()
                        print("Calculating Background distributions for random modules")
                        for(i in 1:length(random_module_sizes)){
                          m_size = random_module_sizes[i]
                          distance <-  Background_Distribution_Modules(RNAseq_Annotated_Matrix,
                                                                       matrix_features,
                                                                       Z_scores,
                                                                       m_size,
                                                                       iterations,
                                                                       parallel_cores
                                                                       )
                          m_char = as.character(m_size)
                          bg_distance_modules[[m_char]] <<- distance
                        }


                        ##Create association matrix
                        # All_modules_pairwise_KO_distances <- P_NRED_Distance_Function(RNAseq_Annotated_Matrix, Z_scores,
                        #                                                               matrix_features, list_of_all_modules_KOs)
                        #
                        # All_clustering_results_P_NRED <-cluster_func(RNAseq_Annotated_Matrix,
                        #                                              All_modules_pairwise_KO_distances$combined,
                        #                                              matrix_features, list_of_all_modules)
                        #
                        # All_association_matrix <<- fill_association_matrix(All_clustering_results_P_NRED,
                        #                                                 matrix_features,
                        #                                                 names(list_of_all_modules))


                      },
                      AssociationMatrix = function(KO_terms_in_module_list, module_list){
                        pairwise_KO_distances <<- P_NRED_Distance_Function(RNAseq_Annotated_Matrix,
                                                                          Z_scores,
                                                                          matrix_features,
                                                                          KO_terms_in_module_list)

                        clustering_results_P_NRED <<- cluster_func(RNAseq_Annotated_Matrix_default_bin,
                                                                  pairwise_KO_distances,
                                                                  matrix_features,
                                                                  module_list
                                                                  )

                        All_association_matrix <<- fill_association_matrix(clustering_results_P_NRED,
                                                                      matrix_features,
                                                                      names(module_list))

                      },
                      plotIndividualBackgroundDist = function(){
                        library("hexbin")
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
                        bg_distance_names <- names(bg_distance_modules)
                        legend <- sapply(bg_distance_names, function(x) paste('N=', x, sep=""))
                        colours <- c("red", "blue", "green", "pink", "black")
                        plot(density(bg_distance_modules[[bg_distance_names[1]]],na.rm=TRUE),col = colours[1],
                             main="Background distributions of modules \nof size \"N\"\nBR",ylim=c(0,1),xlim=c(-4,4),cex.main=.75)
                        for(i in 2:length(bg_distance_names)){
                          points(density(bg_distance_modules[[bg_distance_names[i]]],na.rm=TRUE),
                                 col=colours[i], type="l")
                        }

                        legend("topright",legend=legend,col=colours[1:length(bg_distance_names)],lty=c(1,1,1))
                      }

                    )
                   )


