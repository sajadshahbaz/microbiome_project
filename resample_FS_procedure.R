
U_test_FS<- FS_function(data, decision, lvl=0.05, p.adjust.method="holm"){
pv<- lapply(1:ncol(data),
		function(j)
		wilcox.test(  data[decision==1,j], data[decision==0,j] )$p.value
	   ) |> unlist()
pva<- p.adjust(pv, p.adjust.method)
return(list(relevant.variables=which(pva<lvl)
	    p.value=lv,
	    adjusted.p.value=pva,
	    rel_set= colnames(data)[which(pva<lvl)]
	   ))
}



resample_FS_procedure<- function(List_base_resamples,
				 List_bootrap_indexes,
				 FS_function,
				 FS_on_binary=FALSE,
				 lvl=0.05,
				  p.adjust.method="holm"){
		
	stopifnot(length(List_base_resamples)==length(List_bootrap_indexes))
	rel_sets_base<-list()
	rel_sets_bs<-list()
	taxa_field<- ifelse(FS_on_binary,"taxa_binary","taxa_abundance")
	for (i in seq_along(List_base_resamples))
		{
		   bIDX<- List_bootrap_indexes[[i]]
		   rel_sets_base[[i]]<- FS_function(data=List_base_resamples[[i]][[taxa_field]],
						    decision=List_base_resamples[[i]]$disease_status,
						    lvl=lvl, p.adjust.method=p.adjust.method)
		   rel_sets_bs[[i]]<- FS_function(data=List_base_resamples[[i]][[taxa_field]][bIDX,]
						    decision=List_base_resamples[[i]]$disease_status[bIDX],
						    lvl=lvl, p.adjust.method=p.adjust.method)
		}
	any_times_rel_set<-Reduce(union, lapply(rel_sets_base, function(x) x$rel_set))
	times_rel_base<-Reduce("+", lapply( rel_sets_base, function(x)  any_times_rel_set %in% x) )
	times_rel_bs<-Reduce("+", lapply(rel_sets_bs, function(x) any_times_rel_set %in% x  ))
	features_status<- data.frame( feature_name= any_times_rel_set,
				      n_rel_base= times_rel_base,
				      n_rel_bs=time_rel_bs) 
	return(features_status)
}
