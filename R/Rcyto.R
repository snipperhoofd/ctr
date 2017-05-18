Clean_Data <- function(rules_dataframe){

  replace_chars <- function(val){
    val <- gsub(" ", "", val)
    val <- gsub("\\{", "", val)
    val <- gsub("\\}", "", val)
    return(val)
  }

  lhs_rhs <- str_split_fixed(rules_dataframe[,1], "=> ", 2)
  rules_dataframe <- cbind(lhs_rhs, rules_dataframe[,2:4])
  rules_dataframe[,1] <- sapply(rules_dataframe[,1], replace_chars)
  rules_dataframe[,2] <- sapply(rules_dataframe[,2], replace_chars)
  return(rules_dataframe)
}

getPathwayName <- function(KO_modules_table2){
  terms <- list()
  for(i in 1:nrow(KO_modules_table2)){
    pathway <- KO_modules_table2[i,1]
    term <- KO_modules_table2[i,2]
    terms[[term]] = gsub(" ", "_", pathway)
  }
  return(terms)
}

#' @export
Draw_Network <- function(rules_dataframe, KO_modules_table2, N = 100){
  expanded_rules <- Clean_Data(rules_dataframe)
  expanded_rules <- expanded_rules[order(expanded_rules[,5]), ][1:N,]
  #create network
  g <- new ('graphNEL', edgemode='directed')
  g <-  initNodeAttribute(graph=g, attribute.name = 'pathway_module',
                          attribute.type = 'char',
                          default.value = 'undefined')
  g <-  initEdgeAttribute(graph=g, attribute.name = 'edgeType',
                          attribute.type = 'char',
                          default.value = 'Arrow')
  name_per_module <- getPathwayName(KO_modules_table2)

  for(i in 1:nrow(expanded_rules)){
    lhs <- strsplit(expanded_rules[i,1], ",")[[1]]
    rhs <- strsplit(expanded_rules[i,2], ",")[[1]]

    #Create all lhs nodes if more than one all are added seperate
    for(i in 1:length(lhs)){
      err = tryCatch({
        g <- graph::addNode(lhs, g)
        nodeData (g, lhs[i], 'pathway_module') <- name_per_module[[strsplit(lhs, "\\.")[[1]][1]]]
        },
        error = function(err) {
          return(NA)
        })

      #Create all rhs nodes if more than one all are added seperate
      for(j in 1:length(rhs)){
        err = tryCatch({
          g <- graph::addNode(rhs, g)
          nodeData (g, rhs[j], 'pathway_module') <- name_per_module[[strsplit(rhs, "\\.")[[1]][1]]]
        },
        error = function(err) {
          return(NA)
        })

        #Create the connecting edge
        err = tryCatch({
          g <- graph::addEdge(lhs, rhs, g)
          edgeData(g, lhs[i], rhs[j], 'edgeType') <- 'Associates_with'
        },
        error = function(err){
          return(NA)
        }
        )
      }
    }








  }

  cw <- CytoscapeWindow('vignette', graph = g, overwriteWindow = T)
  cw <-  setGraph (cw, g)
  displayGraph(cw)
  layoutNetwork (cw, layout.name = 'circular')



  top_10_modules <- c("Biosynthesis_of_secondary metabolites", "Cell_signaling",
                    "Two-component_regulatory_system", "Cofactor_and_vitamin_biosynthesis",
                    "Central_carbohydrate_metabolism", "ATP_synthesis", "Ribosome",
                    "Carbon_fixation", "Other_carbohydrate_metabolism", "Drug_resistance",
                    "Drug efflux transporter/pump", "Aromatic amino acid metabolism",
                    "RNA processing", "Bacterial secretion system", "Methane metabolism",
                    "Spliceosome", "Aromatics degradation", "Saccharide, polyol, and lipid transport system",
                    "Other amino acid metabolism", "Serine and threonine metabolism"

  )
  colors <- c("#070D00", "#00FF00", "#341D83", "#483D84", "#5E568D", "#746E9A",
              "#8A86A7", "#A09EB5", "#B6B5C2", "#CBCBCF", "#CBCBCF", "#B4B6C2",
              "#9D9FB5", "#8488A8", "#6C719A", "#525A8E", "#374284", "#062983",
              "#0000FF", "#0F0A00")

  setNodeColorRule (cw, node.attribute.name = 'pathway_module',
                    control.points = top_10_modules,  colors = colors,
                    mode = 'lookup', default.color='#AA0000')

  setEdgeTargetArrowRule(cw, 'edgeType', 'Associates_with', 'Arrow', default= 'Arrow')

  redraw(cw)
}


