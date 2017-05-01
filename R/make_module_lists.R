#' Build module lists
#'
#' A list of all modules in the KEGG database
#'
#' @param KO_modules_table the KEGG database in the structure : WWW, XXX, YYY, ZZZ 
#' @export
#' @return a list of all modules in the KEGG database.
#' @examples list_of_all_modules <- make_module_lists(KO_modules_table)


make_module_lists<- function(KO_modules_table) {
  
  for (i in 1:length(table(KO_modules_table$V1))) {
    module_list<-as.list(names(table(KO_modules_table$V2[which(KO_modules_table$V1%in%names(table(KO_modules_table$V1)))])))
    names(module_list)<-names(table(KO_modules_table$V2[which(KO_modules_table$V1%in%names(table(KO_modules_table$V1)))]))
  }
  for (j in 1:length(module_list)) {
    module_list[[j]]<-KO_modules_table$V4[which(KO_modules_table$V2==module_list[j])]
  }
  return(module_list)
}