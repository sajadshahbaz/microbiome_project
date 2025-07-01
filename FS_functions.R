
#collection of functions
#that as an input take (at least)
#X (samples x features)
#y (samples x 1)
#significance lvl
#p.adjust.method 

#and return:
# relevant.variables - indexes of features in X found relevant wrt y
# p.value - corresponding p-values
# statistic - value of statistic used for each feature to determine p-value
# adjusted.p.value - p-values after multiple test correction
# rel_set - names of features found relevant 

U_test_FS<- function(X, y, lvl=0.05, p.adjust.method="holm"){
data=X
decision=y
ut<- lapply(1:ncol(data),
		function(j){
		res<-wilcox.test(  data[decision==1,j], data[decision==0,j] )
		list(p.value=res$p.value, statistic= res$statistic)
		}
	   ) 
pv<- sapply(ut, function(x) x$p.value)
pva<- p.adjust(pv, p.adjust.method)
stat<- sapply(ut, function(x) x$statistic)
return(list(relevant.variables=which(pva<lvl),
	    statistic=stat,
	    p.value=pv,
	    adjusted.p.value=pva,
	    rel_set= colnames(data)[which(pva<lvl)]
	   ))
}

#helper function to binarize taxa abundance table
# columns -taxons, rows - samples

binarize_taxa<- function(taxonomy, thresholds){
  stopifnot(length(thresholds)==ncol(taxonomy))
  taxonomy<- as.matrix(taxonomy)
tb<- t( t(taxonomy)>thresholds)*1. #exploits recycling mechanism of R:
			      #column [[i]] of `taxonomy` gets binarized by meds[[i]]
return(tb)}

# helper function to filter columns of binary table
filter_sparse_cols<- function(data, thr) {
  #including data==0 matters only for few common taxa in "zero" based binarization
  pmin(colSums(data==0),colSums(data==1))-> minority_count
  data[, minority_count >= thr]
  
}


#helper function for MDFS 2D:
# given input data and result of a feature selection, 
# function finds the best fitting synergistic partner to each
# variable in interesting.vars.

#get_interaction_partner<- function(dataset, decision, interesting.vars, ...){
##
#if (!(all(interesting.vars %in% 1:ncol(dataset))))
#        stop("interesting.vars must be in 1:ncol(dataset)")
#
#if (is.null(colnames(dataset)))
#        stop("datatset should have column names (taxa names)")
#
#best.tuples <- ComputeInterestingTuplesDiscrete(dataset,
#                                        decision,
#                                        interesting.vars =interesting.vars,
#					dimensions=2, ...
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

get_all_interaction_partners<- function(dataset,
					decision,
					interesting.vars,
					MDFS_1D_result,
					synergy_thr=NULL, 
					...){

if (!(all(interesting.vars %in% 1:ncol(dataset))))
        stop("interesting.vars must be in 1:ncol(dataset)")

if (is.null(colnames(dataset)))
        stop("datatset should have column names (taxa names)")

best.tuples <- ComputeInterestingTuplesDiscrete(dataset,
                                        decision,
                                        interesting.vars =interesting.vars,
					dimensions=2,
					...
                                          )

best.tuples<- best.tuples[ best.tuples$Var %in% interesting.vars, ]
  if(is.null(synergy_thr)) synergy_thr = min(MDFS_1D_result$statistic[MDFS_1D_result$relevant.variables])
best.tuples$IG_1D<- MDFS_1D_result$statistic[ best.tuples$Var ]
best.tuples$IG_increase= best.tuples$IG - best.tuples$IG_1D 
best.tuples<- best.tuples [ best.tuples$IG_increase > synergy_thr, ]
data.frame( rel_IDX= best.tuples$Var,
	    partner_IDX= ifelse(best.tuples$Tuple.1==best.tuples$Var,
				best.tuples$Tuple.2,
				best.tuples$Tuple.1))-> res_df
res_df$rel_name<- colnames(dataset)[res_df$rel_IDX]
res_df$partner_name<- colnames(dataset)[res_df$partner_IDX]
res_df$rel_IG<- best.tuples$IG
res_df$rel_IG_1D<- best.tuples$IG_1D
res_df$IG_increase<- best.tuples$IG_increase
res_df
}



#helper function to helper function:
#using output of get_interaction_partner and 1D MDFS,
# function adds additional column to the partner data frame 
# that allow one to decide if synergy is relevant
# IG_increase - how much bigger statistic for each relevant feature grows in 2D 
#               in comparison with 1D
# synergy - TRUE, if IG_increase > synergy_thr (which if left null, is by default 
#						set to minimum value of stat in 1D test
#						    of 1D relevant features)

#find_true_partners<- function(partner_data_frame, MDFS_1D_result,
#			      synergy_thr= NULL
#					       ){
#  stopifnot(all(partner_data_frame$rel_IDX %in% 1:length(MDFS_1D_result$statistic)))
#  stopifnot(all(partner_data_frame$partner_IDX %in% 1:length(MDFS_1D_result$statistic)))
#
#
#  stat_1D_rel2D<- MDFS_1D_result$statistic [ partner_data_frame$rel_IDX ]
#  if(is.null(synergy_thr)) synergy_thr = min(MDFS_1D_result$statistic[MDFS_1D_result$relevant.variables])
#  partner_data_frame$IG_increase= partner_data_frame$rel_IG - stat_1D_rel2D
#  partner_data_frame$synergy= (partner_data_frame$IG_increase > synergy_thr)
#  partner_data_frame
#}

#helper to helper of helper...
#given partner data frame with additional columns as returned from find_true_partners()
# and output of MDFS 1D,
# this function returns a vector of feature names that were found as truly synergistic
# partners, but not relevant in 2D or 1D MDFS analysis by themselves.

#only_partners<- function( annotated_partner_df, MDFS_1D_result  ) {
#	stopifnot(all(c("IG_increase", "synergy") %in% colnames(annotated_partner_df)))
#	true_synergistic<-annotated_partner_df$partner_name [ 
#		annotated_partner_df$synergy ] |> unlist() |> as.character()
#	setdiff(true_synergistic, annotated_partner_df$rel_name)
#}




fit_interactions<- function( partner_data_df, decision,bin_thresholds, bin_data){
	stopifnot(!is.null(names(bin_thresholds)))
	stopifnot(all(colnames(bin_data) %in% names(bin_thresholds)))
	stopifnot( all( union( partner_data_df$rel_name, partner_data_df$partner_name) %in% names(bin_thresholds)))
	pairs<- partner_data_df[, c("rel_name", "partner_name"), drop=FALSE]
	if (nrow(pairs)){
		tholds_pairs<- do.call(rbind, apply(pairs,1, function(PAIR) bin_thresholds[c(PAIR[[1]], PAIR[[2]])],simplify=FALSE ) ) |> as.data.frame()
		res= cbind(pairs, tholds_pairs) |> as.data.frame()
		colnames(res)<-c("rel_name","partner_name","rel_thr","partner_thr")
		inter_mat<-matrix(nrow=nrow(res), ncol= 4)
		colnames(inter_mat)<- c("r0p0","r1p0","r0p1","r1p1")
		for (j in 1:nrow(res))
		{
			rj=bin_data[, res[j,"rel_name"] ]
			pj=bin_data[, res[j,"partner_name"] ]
			inter_mat[j,"r0p0"] = sum((rj==0)&(pj==0)&(decision==1))/sum((rj==0)&(pj==0))
			inter_mat[j,"r1p0"] = sum((rj==1)&(pj==0)&(decision==1))/sum((rj==1)&(pj==0))
			inter_mat[j,"r0p1"] = sum((rj==0)&(pj==1)&(decision==1))/sum((rj==0)&(pj==1))
			inter_mat[j,"r1p1"] = sum((rj==1)&(pj==1)&(decision==1))/sum((rj==1)&(pj==1))
			inter_mat[j,]= rank(inter_mat[j,],ties.method="first")
		}
		inter_mat<- as.data.frame(inter_mat)
		res= cbind(res, inter_mat)
		return(res)	
	} else {
		return(NA)
	}
}


generate_interaction_variables<- function(X, interaction_data){

	sufficient_subset<-union(interaction_data$rel_name, interaction_data$partner_name)
	X_sub<- X[, sufficient_subset, drop=FALSE]
	interactions<-matrix(nrow=nrow(X), ncol= nrow(interaction_data))
	for (j in 1:nrow(interaction_data)){
		labels_j<-interaction_data[j,c("r0p0","r1p0","r0p1","r1p1")]
		names(labels_j)<- c("r0p0","r1p0","r0p1","r1p1")
		rj<- X_sub[, interaction_data[j,"rel_name"] ]
		pj<- X_sub[, interaction_data[j,"partner_name"] ]
		rt<- interaction_data[j,"rel_thr"]
		pt<- interaction_data[j,"partner_thr"]
		interactions[(rj<= rt)&(pj<=pt) ,j]=labels_j[["r0p0"]]
		interactions[(rj> rt)&(pj<=pt) ,j]=labels_j[["r1p0"]]
		interactions[(rj<= rt)&(pj>pt) ,j]=labels_j[["r0p1"]]
		interactions[(rj> rt)&(pj>pt) ,j]=labels_j[["r1p1"]]
	}
	colnames(interactions)<- sapply(1:nrow(interaction_data), function(i) 
						paste0(interaction_data$rel_name[[i]], 
						       "__(x)__", 
						       interaction_data$partner_name[[i]] ))
	return(interactions)
}

#given the result of find_true_partners()
#and binarization thresholds (named vector, one entry per variable, its names must include all partner and relevant variable names
# in annotated_partner_df), function returns a data.frame with one row per each synergisitc pair in the data
# columns 1,2 are names of relevant variable and its partner, while 3,4 are corresponding binarization thresholds (X_b = 1 if X > thr else X_b=0)

#fit_interactions<- function( annotated_partner_df, bin_thresholds){
#	stopifnot( all(union(annotated_partner_df$rel_name, annotated_partner_df$partner_name) %in% names(bin_thresholds)) )
#	pairs<- annotated_partner_df[ annotated_partner_df$synergy, c("rel_name","partner_name"), drop=FALSE ]
#	if (nrow(pairs)){
#		tholds_pairs<- do.call(rbind, apply(pairs,1, function(PAIR) bin_thresholds[c(PAIR[[1]], PAIR[[2]])],simplify=FALSE ) ) |> as.data.frame()
#		res=cbind(pairs, tholds_pairs)	
#		colnames(res)<-c("rel_name","partner_name","rel_thr","partner_thr")
#		return(res) 
#	} else {
#		return(NA)
#	}
#}

# given output of fit_interactions, this function creates a matrix of discrete predictors based on interaction of all pairs present in fitted_interaction_data (as rows)
# names of synthetic variables are REL_NAME__(x)__PARTNER_NAME

#generate_interaction_variables<- function(X, fitted_interaction_data){
#
#	sufficient_subset<-union(fitted_interaction_data$rel_name, fitted_interaction_data$partner_name)
#	X_sub<- X[, sufficient_subset, drop=FALSE]
#	#construct a map from names of variables to their binarization thresholds
#	sapply(colnames(X_sub), function(v_name) if (v_name %in% fitted_interaction_data$rel_name) {
#	       							fitted_interaction_data$rel_thr[ which(fitted_interaction_data$rel_name == v_name)[[1]] ]
#						} else {
#						    fitted_interaction_data$partner_thr[ which(fitted_interaction_data$partner_name== v_name)[[1]] ]
#						}
#	)-> X_sub_tholds
#	#binarize 
#	X_sub_b<- binarize_taxa( X_sub,X_sub_tholds)
#	# use binarized variables to build interactions
#	interactions<- matrix(nrow=nrow(X),ncol=nrow(fitted_interaction_data))	
#	for (i in 1:nrow(fitted_interaction_data))
#		interactions[,i]= interaction( X_sub_b[, fitted_interaction_data$rel_name[[i]] ],
#					       X_sub_b[, fitted_interaction_data$partner_name[[i]] ]
#					     ) |> as.integer()
#	colnames(interactions)<- sapply(1:nrow(fitted_interaction_data), function(i) 
#						paste0(fitted_interaction_data$rel_name[[i]], 
#						       "__(x)__", 
#						       fitted_interaction_data$partner_name[[i]] ))
#	return(interactions)
#}





#function performs 1D & 2D analysis jointly
#for 2 discretization methods for features in X (based on 0, or median of each col in X)
#
# mc - minimum binarized feature class size, if it is less than mc, feature is filtered out.
#
#output of each of the 4 analysis variants
#is stored in a corresponding sublist:
#res_1Dm, res_2Dm, res_1D0, res_2D0 (two for median based and two for zero based)
#
#each sublist contains the fields mentioned at the beginning of the file.
#
#
######moreover, each of the 2D sublists contains additional components:
######	- partner_df: data.frame with one row per each feature in rel_set
######		      where information on the best fitting, synergistic partner is given
######	- partner_set:names of features that were found as synergistic partners,
######		      but not found relevant by themselves

MDFS_FS<- function(X, y, lvl=0.05,
		   p.adjust.method="holm",
		   seed=NULL,
		   mc=30){
   attributes(y)<-NULL
   tholds_m= colMedians(X|> as.matrix()) 
   binarize_taxa(X,tholds_m)|> filter_sparse_cols(thr=mc) -> Xm
   tholds_0= rep(0,ncol(X)) 
   binarize_taxa(X,tholds_0)|> filter_sparse_cols(thr=mc) -> X0
   names(tholds_m)<-names(tholds_0)<- colnames(X)
   to_keep<- c("statistic","p.value","adjusted.p.value","relevant.variables")
   MDFS.discrete(data=Xm,decision = y,dimensions = 1,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=seed)[to_keep]-> res_1Dm
   MDFS.discrete(data=Xm,decision = y,dimensions = 2,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=seed)[to_keep]-> res_2Dm
   MDFS.discrete(data=X0,decision = y,dimensions = 1,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=seed)[to_keep]-> res_1D0
   MDFS.discrete(data=X0,decision = y,dimensions = 2,
                      p.adjust.method =p.adjust.method, 
		      level=0.05,seed=seed)[to_keep]-> res_2D0
  
   res_1Dm$rel_set= colnames(Xm)[res_1Dm$relevant.variables]
   res_1Dm$feature_names= colnames(Xm)
   res_2Dm$rel_set= colnames(Xm)[res_2Dm$relevant.variables]
   res_2Dm$feature_names= colnames(Xm)
   res_1D0$rel_set= colnames(X0)[res_1D0$relevant.variables]
   res_1D0$feature_names= colnames(X0)
   res_2D0$rel_set= colnames(X0)[res_2D0$relevant.variables]
   res_2D0$feature_names= colnames(X0)

   res_2Dm$partner_df = get_all_interaction_partners(dataset=Xm,
					decision=y,
					interesting.vars=res_2Dm$relevant.variables,
					MDFS_1D_result=res_1Dm
					)
   res_2Dm$partner_set<- setdiff( res_2Dm$partner_df$partner_name, res_2Dm$partner_df$rel_name)
   res_2Dm$interaction_data<-fit_interactions(partner_data_df=res_2Dm$partner_df,
					      decision=y,
					      bin_thresholds=tholds_m,
					      bin_data=Xm)


   res_2D0$partner_df = get_all_interaction_partners(dataset=X0,
					decision=y,
					interesting.vars=res_2D0$relevant.variables,
					MDFS_1D_result=res_1D0
					)
   res_2D0$partner_set<- setdiff( res_2D0$partner_df$partner_name, res_2D0$partner_df$rel_name)
   res_2D0$interaction_data<-fit_interactions(partner_data_df=res_2D0$partner_df,
					      decision=y,
					      bin_thresholds=tholds_0,
					      bin_data=X0)


#
#   res_2Dm$partner_df= get_interaction_partner(dataset=Xm, decision=y,
#					       	interesting.vars= res_2Dm$relevant.variables)
#   res_2Dm$partner_df= find_true_partners( partner_data_frame= res_2Dm$partner_df, 
#					   MDFS_1D_result= res_1Dm)
#   res_2Dm$partner_set=only_partners(annotated_partner_df= res_2Dm$partner_df, res_1Dm)
#
#   res_2D0$partner_df= get_interaction_partner(dataset=X0, decision=y,
#					       	interesting.vars= res_2D0$relevant.variables)
#   res_2D0$partner_df= find_true_partners( partner_data_frame= res_2D0$partner_df, 
#					   MDFS_1D_result= res_1D0)
#   res_2D0$partner_set=only_partners(annotated_partner_df= res_2D0$partner_df, res_1D0)
#
#   fit_interactions( annotated_partner_df=res_2Dm$partner_df, bin_thresholds= tholds_m) -> m2D_interaction_data
#   fit_interactions( annotated_partner_df=res_2D0$partner_df, bin_thresholds= tholds_0) -> z2D_interaction_data
#   res_2Dm$interaction_data<- m2D_interaction_data
#   res_2D0$interaction_data<- z2D_interaction_data
#
   RESULT<-list(
     res_1Dm=res_1Dm,
     res_2Dm=res_2Dm,
     res_1D0=res_1D0,
     res_2D0=res_2D0
   )
   
}
