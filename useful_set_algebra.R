
union_and_unique_parts<- function(rel_varNames_list){
  if (!all(unlist(lapply(rel_varNames_list,class))=="character"))
    stop("elements of rel_varNames_list should be character vectors (colnames of variables)")
  all_variables<-unique(unlist(rel_varNames_list))
  lapply(seq_along(rel_varNames_list), function(i) 
                setdiff(rel_varNames_list[[i]], 
                        intersect(rel_varNames_list[[i]],
                                  unlist(rel_varNames_list[-i]))
                        )
        )-> unique_parts
  names(unique_parts)<- names(rel_varNames_list)
  list( union=all_variables,
        unique_parts=unique_parts)
}

jaccard_index<- function(A,B) length(intersect(A,B))/length(union(A,B))

which_most_common_set<- function(rel_varNames_list){
  if (!all(unlist(lapply(rel_varNames_list,class))=="character"))
    stop("elements of rel_varNames_list should be character vectors (colnames of variables)")
  common<- Reduce(intersect,rel_varNames_list)
  if (!length(common))
    stop("intersection is empty!")
  which.max(unlist(lapply(rel_varNames_list,
                function(SET) jaccard_index(SET,common)
  )))
}
