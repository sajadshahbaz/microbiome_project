library(magrittr)
source("draw_once_per_host2.R")
library(matrixStats)

# create resamples according to decision variable
# and compute 1D+2D MDFS on deterministicly binarized data
resample_MDFS<- function(metadata,
			 taxonomy,
			 decision_colname,
			 seeds,
			 n_resamples=30,
			 host_id_colname="host.subject.id",
			 #passed to draw_once_per_host
                   	 min_minority_class_taxa =20, 
			 min_minority_class_dis = 30,
			 binarize_on_median=TRUE,
			 print_summary = FALSE,
                   	 report_progess = TRUE) {
stopifnot(length(seeds)==n_resamples)
y=metadata[,decision_colname]
host_id= metadata[, host_id_colname]
names(y)<-names(host_id)<-rownames(metadata)
f_name= sprintf("%s_base_resamples.rds",decision_colname)

if (!file.exists(f_name)){
resamples<-list()
for (i in 1:n_resamples){
  set.seed(seeds[[i]])
draw_once_per_host(taxonomy=taxonomy, disease = y, host_id = host_id,
                   min_minority_class_taxa =min_minority_class_taxa,
		   min_minority_class_dis = min_minority_class_dis,
		   binarize_on_median=binarize_on_median,
		   print_summary = FALSE,
                   report_progess =report_progess)->resamples[[i]]

if (report_progess) message(sprintf("%d/%d done",i,n_resamples))
}
saveRDS(base_resamples,f_name) } else {
  resamples<-readRDS(f_name)
}
RESULT<- list()
RESULT$resamples=resamples

mdfs_f_name=sprintf("%s_rel_MDFS_group.RData",decision_colname)

if (!file.exists(mdfs_f_name)) {
MDFS_1D<-MDFS_2D<-MDFS_1D_relnames<-MDFS_2D_relnames<-MDFS_1D2D_relnames<-list()
for (i in 1:n_resamples){current_seed=seeds[[i]]
  set.seed(current_seed)
  # Perform MDFS with 1D
  MDFS_1D_result <- MDFS(
    data = resamples[[i]]$taxa_binary,
    decision = resamples[[i]]$disease_status,
    dimensions = 1,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed= current_seed,  # Set the seed for reproducibility
    level = 0.05
  )

  MDFS_1D[[i]] <- MDFS_1D_result
  MDFS_1D_relnames[[i]] <- colnames(resamples[[i]]$taxa_binary)[MDFS_1D_result$relevant.variables]

  # Perform MDFS with 2D
  MDFS_2D_result <- MDFS(
    data = resamples[[i]]$taxa_binary,
    decision = resamples[[i]]$disease_status,
    dimensions = 2,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed = current_seed,
    level = 0.05
  )
  MDFS_2D[[i]] <- MDFS_2D_result
  MDFS_2D_relnames[[i]] <- colnames(resamples[[i]]$taxa_binary)[MDFS_2D_result$relevant.variables]

  # Combine 1D and 2D relevant variables
  MDFS_1D2D_relnames[[i]] <- union(MDFS_1D_relnames[[i]], MDFS_2D_relnames[[i]])
if (report_progess) message(sprintf("%d/%d done",i,n_resamples))
}

# Combine results into lists of lists
MDFS_object_lists <- list(MDFS_1D = MDFS_1D, MDFS_2D = MDFS_2D)
relnames_lists <- list(MDFS_1D_rel = MDFS_1D_relnames, 
                       MDFS_2D_rel = MDFS_2D_relnames, 
                       MDFS_1D2D_rel = MDFS_1D2D_relnames)
# Save results
save(relnames_lists, MDFS_object_lists, 
     file =mdfs_f_name)
} else {load(mdfs_f_name)}

RESULT$relnames_lists=relnames_lists
RESULT$MDFS_object_lists=MDFS_object_lists
return(RESULT)
}

#get_interaction_partner<- function(dataset, decision, interesting.vars, ...){
#
#if (!(all(interesting.vars %in% 1:ncol(dataset))))
#        stop("interesting.vars must be in 1:ncol(dataset)")
#
#if (is.null(colnames(dataset)))
#        stop("datatset should have column names (taxa names)")
#
#best.tuples <- ComputeInterestingTuples(dataset,
#                                        decision,
#                                        interesting.vars =interesting.vars,                                         ...
#                                          )
#do.call(rbind,
#       lapply(interesting.vars, function(V){
#
#               V_subset<- best.tuples$Var == V
#               V_which.max<- which.max( best.tuples[ V_subset,]$IG )
#              m_t1<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.1"]]
#              m_t2<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.2"]]
#              p_IDX<- if( V== m_t1) m_t2 else m_t1
#              data.frame( rel_IDX= V, partner_IDX= p_IDX, rel_name = colnames(dataset)[[V]],                                                                             partner_name= colnames(dataset)[[p_IDX]],
#                 rel_IG= best.tuples[ V_subset, ][V_which.max,][["IG"]] )
#              }
#             )
#       )
#}
