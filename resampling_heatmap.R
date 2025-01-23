N_TRIALS=30
load("randomdataset.RData"); agp_sets<-randomdataset; rm(randomdataset)
readRDS("Bootstrap_IDX.rds")-> bs_IDX_list

read.delim("filled_IG_table.csv", sep=",")-> ref_table
ref_table$name[ ref_table$rel_1D==30 ] -> consensus_1D
ref_table$name -> consensus_AllRelevant

library(matrixStats)
source("pairwise_t_test.R")

sign_matrices_list<-list()

for (B in 1:N_TRIALS) {

B_tnames<-colnames(agp_sets[[B]]$abundance_sampled)

stopifnot(all(consensus_AllRelevant %in% B_tnames))

B_P<- agp_sets[[B]]$binary_taxa_sampled #presence
B_A<- agp_sets[[B]]$abundance_sampled #presence
B_y<- agp_sets[[B]]$y_sampled
B_idx<- bs_IDX_list[[B]]

B_Psampled<- B_P[ B_idx, consensus_AllRelevant ]
B_Asampled<- B_A[ B_idx,consensus_AllRelevant ]
B_ysampled<- B_y[ B_idx ]

pairwise_means_sderrs(abundance = B_Asampled, presence = B_Psampled, 
                      disease = B_ysampled, report_progress = TRUE) -> interaction_components

difference_in_differences(interaction_components, lvl=0.01) -> interaction_matrices

get_sign_matrices(interaction_matrices) -> sign_matrices

for (gr in c("total","d","h"))
	colnames(sign_matrices[[ gr ]])<-rownames(sign_matrices[[ gr ]])<- sign_matrices$taxa_names

sign_matrices_list[[B]]= sign_matrices
}

consensus_heatmap<- function(matrices_list, type="total",min_freq=30){
sign_v<-do.call(rbind, lapply(matrices_list, function(x) 
                  as.vector(x[[type]] )
                ))
 n_p=colSums(sign_v==1)  
 n_n=colSums(sign_v==-1)   
 n_0=-n_p -n_n + nrow(sign_v)
 print(max(n_0))
 rbind(n_p,n_n,n_0)-> sign_counts
 colMaxs(sign_counts)->n_mostfreq
 print(max(n_mostfreq))
 apply(sign_counts,2,function(COL)  c(1,-1,0)[[which.max(COL)]])-> consensus_signs
 consensus_signs[ n_mostfreq < min_freq ] = 0
 map_<- matrix(nrow=nrow( matrices_list[[1]][[type]]),
               ncol=ncol(matrices_list[[1]][[type]])
              )
 map_[,]= consensus_signs
 rownames(map_)<-colnames(map_)<- rownames( matrices_list[[1]][[type]])
 map_
}


consensus_heatmap(sign_matrices_list, "total",27)-> map_t
consensus_heatmap(sign_matrices_list, "d",27)-> map_d
consensus_heatmap(sign_matrices_list, "h",27)-> map_h

saveRDS(list(t=map_t,
	     h=map_h,
	     d=map_d),"sign_matrices_consensus_over_resampling27.rds")
