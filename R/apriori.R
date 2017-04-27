library(arules)
library(arulesViz)

Apriori <- setRefClass("AssociationRules",
                       fields = list(dataset = "matrix"),
                       methods = list(
                         run = function(){
                          return(apriori(dataset, parameter=list(support=0.3, confidence=0.5)))},
                         plot_graph = function(rules, n = 50 ){
                           plot(head(sort(rules, by = "lift"), n), method = "graph", control=list(cex=.8))},
                         get_top10 = function(rules){
                           return(inspect(head(sort(rules), n=10)))},
                         Data_clean = function(){
                           dataset <<- apply(dataset, 2, as.factor)
                           extendCollumns(dataset[,2])


                         }
                       )
)
extendCollumns = function(column){
    print(levels(column))
}


test <- Apriori$new(dataset = c1_data)
rules <- test$run()
test$plot_graph(rules)
test$get_top10(rules)


test <- Apriori$new(dataset = dat)
test$Data_clean()






#load a table
dat <- as.matrix(read.csv("/home/joris/Desktop/module_cluster_matrix.csv"))
rownames(dat) <- dat$Bins
dat$Bins <- NULL


data <- apply(dat,2, as.factor)
#testcluster1
c1_data <- Data_clean(data, 1)
ret <- AssociationRules(c1_data)

plot(itemsets$top10_rules)

