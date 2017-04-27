#apriori
library(arules)
library(arulesViz)
#rock
library(cba)

Apriori <- setRefClass("AssociationRules",
                       fields = list(dataset = "data.frame",
                                     flatdata = "matrix"),
                         methods = list(
                         run = function(){
                          return(apriori(dataset, parameter=list(support=0.3, confidence=0.5)))},
                         plot_graph = function(rules, n = 50 ){
                           plot(head(sort(rules, by = "lift"), n), method = "graph", control=list(cex=.8))},
                         get_top10 = function(rules){
                           return(inspect(head(sort(rules), n=10)))},
                         Data_clean = function(){

                            extendCollumns = function(column, starting_matrix){
                                cluster_to_true_false = function(val, clusternum){
                                  if(is.na(val) || val != clusternum){
                                    return(0)
                                  }else{
                                    return(1)
                                  }
                                }
                                
                              
                              for(i in levels(as.factor(column))){
                                vec <- column
                                starting_matrix <- cbind(starting_matrix, sapply(vec, cluster_to_true_false, clusternum=i))
                              }
                              return(starting_matrix)
                            }
                           starting_matrix <- matrix(nrow = nrow(dataset))
                           for(i in ncol(dataset)){
                             starting_matrix <- extendCollumns(dataset[,i], starting_matrix)
                           }
                           #starting_matrix <- starting_matrix[,2: ncol(starting_matrix)]
                           print(starting_matrix)

                         }
                       )
)




test <- Apriori$new(dataset = dat)
test$Data_clean()
levels(as.factor(test$cleandata[,2]))


#load a table
dat <- read.csv("/home/joris/Downloads/module_cluster_matrix2.csv")
rownames(dat) <- dat[,1]
dat$Bins <- NULL


data <- apply(dat,2, as.factor)
#testcluster1
c1_data <- Data_clean(data, 1)
ret <- AssociationRules(c1_data)

plot(itemsets$top10_rules)

