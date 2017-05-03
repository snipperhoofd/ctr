"
Example:

Antecedent 	Consequent
A 	          0
A           	0
A 	          1
A  	          0
B 	          1
B 	          0
B 	          1
                  Support      Confidence     Lift
rule 1: A -> 0      3/7           3/4       (3/4)/(4/7)
rule 2: B -> 1      2/7           2/3       (2/3)/(3/7)

Support:
 How frequent an itemset appears in the dataset.
 The proportion of transactions wich contains itemset X

Confidence:
 How often the rule has been found to be true
 The confidence value of a rule, X -> Y is the proportion of transactions that contains X
 which also contains Y.

Lift:
  The ratio of the observed support to the expcted support if X and Y are independent.
  If the lift is > 1, like it is here for Rules 1 and 2,
    that lets us know the degree to which those two occurrences are dependent on one another,
    and makes those rules potentially useful for predicting the consequent in future data sets

Conviction
  The frequency the rule makes an incorrect prediction.
  conv(X -> Y) = (1-supp(Y)) / (1-conf(X -> Y))
"
Apriori <- setRefClass("AssociationRules",
                        fields = list(dataset = "matrix",
                                     ruleset = "list",
                                     rules_dataframe = "data.frame"),
                        methods = list(
                          run_apriori = function(supp, conf,
                                                 targets = NA){
                            #testing
                            getLineNumbers = function(vec){
                              out = c()
                              for(i in 1:length(vec)){
                                if(vec[i]){
                                  out = c(out,i)
                                }
                              }
                              return(out)
                            }
                            filter_general_redundancy = function(rule){
                              redundant <- is.redundant(rule)
                              redundant_line_nums <- getLineNumbers(redundant)
                              return(rule[-redundant_line_nums])
                            }
                            filter_mirrored_redundancy = function(rule){
                              redundant <- duplicated(generatingItemsets(rule))
                              redundant_line_nums <- getLineNumbers(redundant)
                              return(rule[-redundant_line_nums])
                            }
                            rule = NULL
                            if(is.na(targets)){
                              rule = apriori(dataset,
                                             parameter=list(support=supp,
                                                            confidence=conf))
                            }else{
                              rule = apriori(dataset,
                                             parameter=list(support=supp,
                                                            confidence=conf),
                                             appearance = list(lhs=targets, default= 'rhs'))
                            }
                            original_size = length(rule)
                            #rule = filter_general_redundancy(rule)
                            #rule = filter_mirrored_redundancy(rule)

                            if(original_size > length(rule)){
                              print(paste("Removed", (original_size - length(rule)), "redundant rows"))
                            }
                            ruleset <<- list("rules" = rule)
                          },
                          plot_graph = function(n = 20){
                            plot(head(sort(ruleset$rules, by = "support"), n), method = "graph", control=list(cex=.8), interactive=F)
                          },
                          get_topN = function(n=20){
                            return(inspect(head(sort(ruleset$rules, by="confidence"), n)))
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
                              dataset <<- starting_matrix[,2: ncol(starting_matrix)]
                            }else #this else is for testing purposes with pre-formatted data
                              dataset <<- as.matrix(dataset)
                          }
                       )
)
library(arules)
library(arulesViz)


#load a table
dat <- read.csv("/home/joris/Downloads/expanded_matrix_all_modules.csv")
rownames(dat) <- dat[,1]
dat[,1] <- NULL

associationSearch <- Apriori$new(dataset = as.matrix(dat))
associationSearch$run_apriori(supp = 0, conf = 0.5, targets = c('M00580_1'))

associationSearch$get_topN(20)
associationSearch$plot_graph(20)





