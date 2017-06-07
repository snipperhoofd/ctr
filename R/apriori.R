#' @export Apriori
#' @exportClass AssociationRules
Apriori <- setRefClass(
  "AssociationRules",
  fields = list(dataset = "matrix",
                tmp = 'list',
                rules_dataframe = "data.frame"),
  methods = list(
    run_apriori = function(supp,
                           conf,
                           antecedents = NA,
                           consequents = NA) {
      library(arules)
      Data_clean()
      rule = NULL
      #Check if there are antecedents or consequents queried
      get_rules <- function(lhs = NA, rhs = NA) {
        if (is.na(lhs) &&
            is.na(rhs)) {
          rule = apriori(dataset,
                         parameter = list(support = supp,
                                          confidence = conf))
        } else if (!is.na(lhs)) {
          rule = apriori(
            dataset,
            parameter = list(support = supp,
                             confidence = conf),
            appearance = list(lhs = lhs, default = 'rhs')
          )
        } else if (!is.na(rhs)) {
          rule = apriori(
            dataset,
            parameter = list(support = supp,
                             confidence = conf),
            appearance = list(rhs = rhs, default = 'lhs')
          )
        }
        return(rule)
      }

      rules_dataframe <<- data.frame()
      if(!is.na(consequents) && !is.na(antecedents)){
        rules_dataframe <<- rbind(rules_dataframe,
                                  as(get_rules(lhs = consequents),
                                     'data.frame'))
        print(dim(rules_dataframe))
        print("FIRST")
        rules_dataframe <<- rbind(rules_dataframe,
                                 as(get_rules(rhs = antecedents),
                                    'data.frame'))
        print(dim(rules_dataframe))
        print("SECOND")
      }
      else{
        tmp[['rule']] <<- get_rules(lhs = consequents,
                          rhs = antecedents)
        rules_dataframe <<-rbind(rules_dataframe,
                                 as(tmp[['rule']], 'data.frame'))
      }
    },
    get_topN = function(n = 20, by = 'support') {
      return("Not implemented at this moment")
      #return(head(sort(
      #  rules_dataframe, by = by
      #), n))
    },
    Data_clean = function(split = "sbs") {
      if (split == "sbs") {
        extendCollumns <- function(column, starting_matrix, column_name) {
          cluster_to_true_false <-  function(val, clusternum) {
            if (is.na(val) || val != clusternum)
              return(0)
            else
              return(1)
          }
          for (i in levels(as.factor(column))) {
            changed_column <-
              sapply(column, cluster_to_true_false, clusternum = i)
            starting_matrix <-
              cbind(starting_matrix, changed_column)
            colnames(starting_matrix)[ncol(starting_matrix)] <-
              paste(column_name, ".", i, sep = "")
          }
          return(starting_matrix)
        }
        starting_matrix <- matrix(nrow = nrow(dataset))
        COLUMNS <- colnames(dataset)
        for (i in 1:ncol(dataset)) {
          starting_matrix <-
            extendCollumns(dataset[, i], starting_matrix, COLUMNS[i])
        }
        dataset <<-
          starting_matrix[, 2:ncol(starting_matrix)]
      } else
        #this else is for testing purposes with pre-formatted data
        dataset <<- as.matrix(dataset)
    }
  )
)
