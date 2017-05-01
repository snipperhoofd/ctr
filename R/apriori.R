#apriori
library(arules)
library(arulesViz)

Apriori <- setRefClass("AssociationRules",
                        fields = list(dataset = "data.frame",
                                     flatdata = "matrix",
                                     ruleset = "list"),
                        methods = list(
                          run_apriori = function(){
                            ruleset <<- list("rules" = apriori(flatdata, parameter=list(support=0.01, confidence=0.5)))
                          },
                          plot_graph = function(n = 50){
                            plot(head(sort(ruleset$rules, by = "lift"), n), method = "graph", control=list(cex=.8), interactive=T)
                          },
                          get_top10 = function(){
                            return(inspect(head(sort(ruleset$rules), n=10)))
                        },
                          Data_clean = function(split="sbs"){
                            if(split == "sbs"){
                              extendCollumns <- function(column, starting_matrix, column_name){

                                cluster_to_true_false <-  function(val, clusternum){
                                if(is.na(val) || val != clusternum) return(0)
                                  else return(1)
                                }

                                for(i in levels(as.factor(column))){
                                  changed_column <- sapply(column, cluster_to_true_false, clusternum=i)
                                  starting_matrix <- cbind(starting_matrix, changed_column)
                                  colnames(starting_matrix)[ncol(starting_matrix)] <- paste(column_name, ".", i, sep="")
                                }
                                return(starting_matrix)
                              }
                              starting_matrix <- matrix(nrow = nrow(dataset))
                              COLUMNS <- colnames(dataset)
                              for(i in 1:ncol(dataset)){
                                starting_matrix <- extendCollumns(dataset[,i], starting_matrix, COLUMNS[i])
                              }
                              flatdata <<- starting_matrix[,2: ncol(starting_matrix)]
                            }else #this else is for testing purposes with pre-formatted data
                              flatdata <<- as.matrix(dataset)
                          }
                       )
)




test <- Apriori$new(dataset = dat)

test$Data_clean("")
test$run_apriori()
test$get_top10()
test$plot_graph(20)



#load a table
#dat <- read.csv("/home/joris/Downloads/module_cluster_matrix2.csv")
dat <- read.csv("/home/joris/Downloads/testmatrix.csv")
rownames(dat) <- dat[,1]
dat[,1] <- NULL



