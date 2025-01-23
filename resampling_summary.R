N_TRIALS=30
load("randomdataset.RData"); agp_sets<-randomdataset; rm(randomdataset)
readRDS("Bootstrap_IDX.rds")-> bs_IDX_list

read.delim("filled_IG_table.csv", sep=",")-> ref_table
ref_table$name[ ref_table$rel_1D==30 ] -> consensus_1D
ref_table$name -> consensus_AllRelevant

#frequencies of all taxa in consensus set

resampling_df<- ref_table

resampling_df[, colnames(resampling_df)!="name"]<- NA

fs_fpattern<-"MDFS_results/FS_replicate_%d.rds"
mdfs_fpattern<-"MDFS_results/MDFS_replicate_%d.rds"
rf_fpattern<- "RF_models/RF_replicate_%d.rds" 

library(matrixStats)
library(randomForest)
library(pROC)

tables_trials<-list()
CM1_trials<-list()
CMall_trials<- list()
RF1_trials<-list()
RFall_trials<-list()
yoob_RF1pred_RFallpred<-list()

for (B in 1:N_TRIALS) {

B_tnames<-colnames(agp_sets[[B]]$abundance_sampled)

stopifnot(all(consensus_AllRelevant %in% B_tnames))
# what are indexes of columns of consensus_AllRelevant in data _B?
match( consensus_AllRelevant, B_tnames)-> AR_position

B_P<- agp_sets[[B]]$binary_taxa_sampled #presence
B_A<- agp_sets[[B]]$abundance_sampled #presence
B_y<- agp_sets[[B]]$y_sampled
B_idx<- bs_IDX_list[[B]]

B_Psampled<- B_P[ B_idx, ]
B_Asampled<- B_A[ B_idx, ]
B_ysampled<- B_y[ B_idx ]

B_oob<- setdiff( 1:nrow(B_A), unique(B_idx))

B_Poob<- B_P[ B_oob, ]
B_Aoob<- B_A[ B_oob, ]
B_yoob<- B_y[ B_oob ]

B_mdfs_fname<- sprintf( mdfs_fpattern, B)
B_fs_fname<- sprintf( fs_fpattern, B)
B_rf_fname<- sprintf( rf_fpattern, B)
B_MDFS<- readRDS( B_mdfs_fname)
B_FS<- readRDS( B_fs_fname)
B_RF<- readRDS(B_rf_fname)$RF_all
B_RF1<- readRDS(B_rf_fname)$RF_1D

stopifnot( all(rownames(B_RF$importance)==consensus_AllRelevant))

B_df<- data.frame(
name= consensus_AllRelevant,
mean= colMeans(B_Asampled[,consensus_AllRelevant]),
freq= colMeans(B_Psampled[,consensus_AllRelevant]),
mean_d= colMeans(B_Asampled[B_ysampled==1,consensus_AllRelevant]),
freq_d= colMeans(B_Psampled[B_ysampled==1,consensus_AllRelevant]),
mean_h= colMeans(B_Asampled[B_ysampled==0,consensus_AllRelevant]),
freq_h= colMeans(B_Psampled[B_ysampled==0,consensus_AllRelevant]),
rel_1D= 1.*(consensus_AllRelevant %in% B_FS$rel_1D),
rel_2D= 1.*(consensus_AllRelevant %in% B_FS$rel_2D),
IG_1D= B_MDFS[[1]]$statistic[ AR_position ],
IG_2D= B_MDFS[[2]]$statistic[ AR_position ],
RF= B_RF$importance[,3],
RF_d= B_RF$importance[,'1'],
RF_h= B_RF$importance[,'0']
)

tables_trials[[B]]<- B_df

B_RFall_prob<- predict( B_RF, newdata= B_Poob, type="prob")[,2]
B_RF1_prob<- predict( B_RF1, newdata= B_Poob, type="prob")[,2]
yoob_RF1pred_RFallpred[[B]] = cbind( y_oob = B_yoob,
				   rf1_pred_oob= B_RF1_prob,
				   rfa_pred_oob= B_RFall_prob)


B_RFall_cl<-  B_RFall_prob>0.5
B_RF1_cl<- B_RF1_prob>0.5

cm_all<- table(Predicted=B_RFall_cl,Actual=B_yoob)
cm_1<- table(Predicted=B_RF1_cl,Actual=B_yoob)

CM1_trials[[B]]<- cm_1
CMall_trials[[B]]<- cm_all

odds_ratio<- function(CM) {(CM[2,2]*CM[1,1])/(CM[1,2]*CM[2,1]) } 



RFall_score<- c(auc= auc(response=B_yoob, predictor= B_RFall_prob),
		OR=odds_ratio(cm_all),
		acc= sum( B_yoob == B_RFall_cl)/length(B_yoob)
)
RF1_score<- c(auc= auc(response=B_yoob, predictor= B_RF1_prob),
	      OR=odds_ratio(cm_1),
		acc= sum( B_yoob == B_RF1_cl)/length(B_yoob)
	      )

RF1_trials[[B]]<- RF1_score
RFall_trials[[B]]<- RFall_score

}

mean_tab<- tables_trials[[1]]
mean_tab[,2:ncol(mean_tab)] = Reduce("+",
				     lapply(tables_trials,
				function(x) x[,2:ncol(x)])
				    )
mean_CM_all<- Reduce("+", CMall_trials)
mean_CM_1D<- Reduce("+", CM1_trials)

dont_divide<-  colnames(mean_tab) %in% c("name", "rel_1D","rel_2D")



mean_tab[, !dont_divide] = mean_tab[, !dont_divide]/N_TRIALS

write.csv(mean_tab,"resampling_mean_table.csv")
saveRDS(yoob_RF1pred_RFallpred, "list_of_OOB_rf1pr_rfALLpr.rds")

print(length(agp_sets[[1]]$y_sampled))
print("mean stats for RF on 1D:")
print(colMeans(do.call(rbind,RF1_trials)) )
print(colSds(do.call(rbind,RF1_trials)) )
print(mean_CM_1D)
print(sum(mean_CM_1D))

print("mean stats for RF on 1,2D and partners:")
print(colMeans(do.call(rbind,RFall_trials)) )
print(colSds(do.call(rbind,RFall_trials)) )
print(mean_CM_all)
print(sum(mean_CM_all))






