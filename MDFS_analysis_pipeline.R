
base_MDFS_end2end<-function(baseResamples,
			    variant_name){
	MDFS_filename=sprintf("%s_rel_MDFS_group.RData" ,variant_name)
	if (!file.exists(MDFS_filename)) {
	MDFS_1D<-MDFS_2D<-MDFS_1D_relnames<-MDFS_2D_relnames<-MDFS_1D2D_relnames<-list()
	for (i in 1:30){current_seed=i
	gc()
	  set.seed(current_seed)
	  # Perform MDFS with 1D
	  MDFS_1D_result <- MDFS.discrete(
	    data = baseResamples[[i]]$taxa_binary,
	    decision = baseResamples[[i]]$disease_status,
	    dimensions = 1,
	    seed= current_seed,  # Set the seed for reproducibility
	    level = 0.05
	  )

	  MDFS_1D[[i]] <- MDFS_1D_result
	  MDFS_1D_relnames[[i]] <- colnames(baseResamples[[i]]$taxa_binary)[MDFS_1D_result$relevant.variables]
	gc()
	set.seed(current_seed)
	  # Perform MDFS with 2D
	  MDFS_2D_result <- MDFS.discrete(
	    data = baseResamples[[i]]$taxa_binary,
	    decision = baseResamples[[i]]$disease_status,
	    dimensions = 2,
	    seed = current_seed,
	    level = 0.05
	  )
	  MDFS_2D[[i]] <- MDFS_2D_result
	  MDFS_2D_relnames[[i]] <- colnames(baseResamples[[i]]$taxa_binary)[MDFS_2D_result$relevant.variables]

	  # Combine 1D and 2D relevant variables
	  MDFS_1D2D_relnames[[i]] <- union(MDFS_1D_relnames[[i]], MDFS_2D_relnames[[i]])
	message(sprintf("%d/30 done",i))
	}
	MDFS_object_lists <- list(MDFS_1D = MDFS_1D, MDFS_2D = MDFS_2D)
	relnames_lists <- list(MDFS_1D_rel = MDFS_1D_relnames, 
			       MDFS_2D_rel = MDFS_2D_relnames, 
			       MDFS_1D2D_rel = MDFS_1D2D_relnames)
	# Save results
	save(relnames_lists, MDFS_object_lists, 
	     file =MDFS_filename)
	} else { message("loading data");load(MDFS_filename)}
	message(sprintf("MDFS stage done"))	
	any_times_relevant_taxa<- Reduce(union, relnames_lists$MDFS_1D2D_rel)
	if (!length(any_times_relevant_taxa)) 
		stop(sprintf("no relevant variables found for %s in %d repeats",
			     variant_name,length(baseResamples)))
	rel_freq<-rel_freq_1D<-rel_freq_2D<-rep(0, length(any_times_relevant_taxa))
	names(rel_freq)<-names(rel_freq_1D)<-names(rel_freq_2D)<-any_times_relevant_taxa
	for (i in seq_along(relnames_lists[[1]])) {
	rel_freq<- rel_freq + ( any_times_relevant_taxa %in% relnames_lists$MDFS_1D2D_rel[[i]] )*1. 
	rel_freq_1D<- rel_freq_1D + ( any_times_relevant_taxa %in% relnames_lists$MDFS_1D_rel[[i]] )*1. 
	rel_freq_2D<- rel_freq_2D + ( any_times_relevant_taxa %in% relnames_lists$MDFS_2D_rel[[i]] )*1. 
	}
	which_stay<-(rel_freq==30)
	base_relevant_set<- any_times_relevant_taxa[which_stay]
	rel_freq<- rel_freq[which_stay]
	rel_freq_1D<-rel_freq_1D[which_stay]
	rel_freq_2D<-rel_freq_2D[which_stay]
	partners_per_run<- list()
	i_set<-1
	message("finding synergistic partners...")
	for (i_set in seq_along(baseResamples))
	{
	partners_per_run[[i_set]] <-
	get_interaction_partner(dataset=baseResamples[[i_set]]$taxa_binary,
				 decision=baseResamples[[i_set]]$disease_status,
			     interesting.vars=MDFS_object_lists[[2]][[i_set]]$relevant.variables ,
				 dimensions=2)

	}
	do.call(rbind, partners_per_run)-> partners_all_runs #all runs together
	aggregate(partners_all_runs$rel_IG, #value to aggregate
		  list(rel=partners_all_runs$rel_name,       #factors over which to group
		      partner=partners_all_runs$partner_name),
		  mean #how to aggregate
		  )-> partner_IG_tabulation


	lapply(base_relevant_set[rel_freq_2D==30], function(taxon)
	{
		taxon_subset<- partner_IG_tabulation$rel == taxon
		if (any(taxon_subset))
		{
			taxon_partners<-partner_IG_tabulation[ taxon_subset, ,drop=FALSE]
			taxon_partners$partner[[ which.max(taxon_partners$x) ]]
		} else NA

	}       ) %>% unlist() -> best_partners_baseRelSet

	names(best_partners_baseRelSet)<- base_relevant_set[rel_freq_2D==30]
	union(base_relevant_set, unname(best_partners_baseRelSet))-> overall_rel_list
	  
	#ugly, but necessary check:
	##ensure that all resamples contain all taxa in overall_rel_list
	stopifnot(
	 lapply(baseResamples, function(x) 
	   all(overall_rel_list %in% colnames(x$taxa_binary)) 
	   ) |> unlist()  |> all() 
	)
	message("computing mean IG")

	baseResamples_present_taxa<- lapply(baseResamples, function(x)
					      colnames(x$taxa_binary))

	overall_rel_IG1D<- mean_IG_overRuns(MDFS_object_lists[[1]],
					    overall_rel_list,
					    baseResamples_present_taxa
					    )
	overall_rel_IG2D<- mean_IG_overRuns(MDFS_object_lists[[2]],
					    overall_rel_list,
					    baseResamples_present_taxa
					    )
	print(length(overall_rel_IG1D))
	print(length(overall_rel_IG2D))
	print(length(overall_rel_list))
	names(overall_rel_IG1D)<-names(overall_rel_IG2D)<-overall_rel_list
	message("assembling data.frame")
	data.frame( name=overall_rel_list,
		    IG1D=overall_rel_IG1D, IG2D= overall_rel_IG2D)-> ig_df
	ig_df$rel_1D<- rel_freq_1D[ ig_df$name ]
	ig_df$rel_2D<- rel_freq_2D[ ig_df$name ]
	ig_df$rel_1D[ is.na(ig_df$rel_1D) ] = 0
	ig_df$rel_2D[ is.na(ig_df$rel_2D) ] = 0
	ig_df$partner_name<- best_partners_baseRelSet[ ig_df$name ]
	ig_df$partner_only<- ig_df$name %in% setdiff(best_partners_baseRelSet, base_relevant_set)
	synergy_thr<- min(ig_df$IG1D[ig_df$rel_1D==30])
	print("minimal increase of IG in 2D analysis for estabilishing synergy:")
	print(synergy_thr)
	IG_increase= ig_df$IG2D - ig_df$IG1D
	# find the true synergies using the thredhold and IG_increase of primary 2D rel. variable:
	true_synergistic_partners<-vector(mode="character")
	for (partner_only in ig_df$name[ig_df$partner_only])
	{
	  if (max(IG_increase[which(ig_df$partner_name== partner_only) ]) > synergy_thr)
	    true_synergistic_partners[[ length(true_synergistic_partners)+1 ]] = partner_only
	}

	#for the rest of the analysis:
	# keep taxa that were relevant 30 times in 2D or 1D and true synergistic partners.
	taxa_to_keep_mask<- (ig_df$name %in% true_synergistic_partners) | 
			      (ig_df$rel_1D==30 ) | 
			      (ig_df$rel_2D==30)
	ig_df<- ig_df[ taxa_to_keep_mask,]
	ig_df<- ig_df[ order(-ig_df$IG2D),]

	message("mean abundance and frequency calculation")
	  
	summary_df<- ig_df
	summary_df$name-> interesting_taxa
	#ensure all baseResamples have same number of rows 
	stopifnot( all(
		(lapply(baseResamples, function(x) nrow(x$taxa_abundance)) %>%
		    unlist() )== nrow(baseResamples[[1]]$taxa_abundance) 
		      )
		 )

	#common basis for all interesting taxa
	# (containers for sums)
	A_h<-A_d<-P_h<-P_d<-P<-A<- rep(0, length(interesting_taxa))
	names(P)<-names(A)<-interesting_taxa
	names(P_d)<-names(A_d)<-interesting_taxa
	names(P_h)<-names(A_h)<-interesting_taxa

	#sum means over resamples
	for (i in seq_along(baseResamples)){
	    
	    intr_present_resample<- intersect(interesting_taxa,
					colnames(baseResamples[[i]]$taxa_abundance)
					    )

	    P_rsmp<-baseResamples[[i]]$taxa_binary[, intr_present_resample]
	    A_rsmp<-baseResamples[[i]]$taxa_abundance[, intr_present_resample]
	    y_rsmp<-baseResamples[[i]]$disease_status

	    P[ intr_present_resample ] =  P[ intr_present_resample ] +
		    colMeans(P_rsmp)
	    A[ intr_present_resample ] =  A[ intr_present_resample ] +
		    colMeans(A_rsmp)

	    P_d[ intr_present_resample ] =  P_d[ intr_present_resample ] +
		    colMeans(P_rsmp[y_rsmp==1,])
	    A_d[ intr_present_resample ] =  A_d[ intr_present_resample ] +
		    colMeans(A_rsmp[y_rsmp==1,])

	    P_h[ intr_present_resample ] =  P_h[ intr_present_resample ] +
		    colMeans(P_rsmp[y_rsmp==0,])
	    A_h[ intr_present_resample ] =  A_h[ intr_present_resample ] +
		    colMeans(A_rsmp[y_rsmp==0,])
	}
	#.. and divide by 30
	P<-P/30
	A<-A/30
	P_h<-P_h/30
	A_h<-A_h/30
	P_d<-P_d/30
	A_d<-A_d/30
	summary_df$freq= P
	summary_df$ab= A
	summary_df$freq_d= P_d
	summary_df$ab_d= A_d
	summary_df$freq_h= P_h
	summary_df$ab_h= A_h
	message("all done!")
	return(summary_df)
}


bootstrap_MDFS<- function(b_idx, baseResamples,seeds, variant_name,base_allRel){
	n_trials= length(b_idx)
	f_pattern<-"MDFS_results/%s_MDFS_replicate_%d.rds"
	f_pattern2<-"MDFS_results/%s_FS_replicate_%d.rds"
	for (B in 1:n_trials){
		IDX_B<-  b_idx[[B]]
		set_B<- baseResamples[[B]]
		set.seed(seeds[[B]])
		f_B<- sprintf(f_pattern,variant_name,B)
		if (!file.exists(f_B)){
		#print(set_B$taxa_binary[IDX_B,] |> dim())
		gc()
		MDFS_1D<- MDFS.discrete(
		data = set_B$taxa_binary[IDX_B,],
		decision = set_B$disease_status[IDX_B],
		dimensions = 1,
		seed= seeds[[B]], 
		level = 0.05
		)
		#print("1D done")
		rel_1D<- colnames(set_B$taxa_binary)[MDFS_1D$relevant.variables]
		# Perform MDFS with 2D
		gc()
		MDFS_2D<- MDFS.discrete(
		data = set_B$taxa_binary[IDX_B,],
		decision = set_B$disease_status[IDX_B],
		dimensions = 2,
		seed = seeds[[B]],
		level = 0.05
		)
		rel_2D<- colnames(set_B$taxa_binary)[MDFS_2D$relevant.variables]
		rel_1or2D<-  union( rel_1D, rel_2D)  
		FS<- list(rel_1D= rel_1D,
		    rel_2D= rel_2D,
		    rel_1or2D=rel_1or2D)
		saveRDS( FS, sprintf(f_pattern2,variant_name, B ) )
		saveRDS( list(MDFS_1D=MDFS_1D,
			MDFS_2D=MDFS_2D),
		    sprintf(f_pattern,variant_name, B ) )
    		message(sprintf("%d/%d done", B, n_trials))
		}
	}; message("MDFS stage done")
	message("aggregating...")
	rel_2D_trials<-rel_1D_trials<- list()  #for indicators of relevance per each run
	for (B in 1:n_trials){
		B_fs_fname<- sprintf( f_pattern2,variant_name, B)
		B_FS<- readRDS( B_fs_fname)
		print(length(B_FS$rel_1D))
		print(length(B_FS$rel_2D))
		rel_1D_trials[[B]]= 1.*(base_allRel %in% B_FS$rel_1D)
		rel_2D_trials[[B]]= 1.*(base_allRel %in% B_FS$rel_2D)
	}
	Reduce("+", rel_1D_trials)->times_rel1D_BS
	Reduce("+", rel_2D_trials)->times_rel2D_BS
	message("done!")
	data.frame(name=base_allRel,
		   BS_rel1D=times_rel1D_BS,
		   BS_rel2D=times_rel2D_BS)
}





