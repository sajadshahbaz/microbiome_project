
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

binarize_taxa<- function(taxonomy, binarize_on_median=TRUE){
  taxonomy<- as.matrix(taxonomy)
meds= if (binarize_on_median) colMedians(taxonomy) else rep(0,ncol(taxonomy))
tb<- t( t(taxonomy)>meds )*1. #exploits recycling mechanism of R:
			      #column [[i]] of `taxonomy` gets binarized by meds[[i]]
return(tb)}

# helper function to filter columns of binary table
filter_sparse_cols<- function(data, thr) {
  #including data==0 matters only for few common taxa in "zero" based binarization
  pmin(colSums(data==0),colSums(data==1))-> minority_count
  data[, minority_count >= thr]
  
}



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
#====================================================================================
#====================================================================================
#====================================================================================
### @TO DO: estabilish a way of determining the synergy threshold from one run!
#	    right now, it is only possible to include below only "post-hoc",
#		after we estabilish robust set of 1D and 2D relevant features:
#
######moreover, each of the 2D sublists contains additional components:
######	- partner_df: data.frame with one row per each feature in rel_set
######		      where information on the best fitting, synergistic partner is given
######	- partner_set:names of features that were found as synergistic partners,
######		      but not found relevant by themselves
#====================================================================================
#====================================================================================
#====================================================================================

MDFS_FS<- function(X, y, lvl=0.05,
		   p.adjust.method="holm",
		   mc=30){
   attributes(y)<-NULL
   binarize_taxa(X,binarize_on_median = TRUE)|> filter_sparse_cols(thr=mc) -> Xm
   
   binarize_taxa(X,binarize_on_median = FALSE)|> filter_sparse_cols(thr=mc) -> X0
   to_keep<- c("statistic","p.value","adjusted.p.value","relevant.variables")
   MDFS.discrete(data=Xm,decision = y,dimensions = 1,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=123)[to_keep]-> res_1Dm
   MDFS.discrete(data=Xm,decision = y,dimensions = 2,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=123)[to_keep]-> res_2Dm
   MDFS.discrete(data=X0,decision = y,dimensions = 1,
                      p.adjust.method = p.adjust.method,
		      level=0.05,seed=123)[to_keep]-> res_1D0
   MDFS.discrete(data=X0,decision = y,dimensions = 2,
                      p.adjust.method =p.adjust.method, 
		      level=0.05,seed=123)[to_keep]-> res_2D0
  
   res_1Dm$rel_set= colnames(Xm)[res_1Dm$relevant.variables]
   res_2Dm$rel_set= colnames(Xm)[res_2Dm$relevant.variables]
   res_1D0$rel_set= colnames(X0)[res_1D0$relevant.variables]
   res_2D0$rel_set= colnames(X0)[res_2D0$relevant.variables]

#   res_2Dm$partner_df= get_interaction_partner(dataset=Xm, decision=y,
#					       	interesting.vars= res_2Dm$relevant.variables)
#   res_2Dm$partner_df= find_true_partners( partner_data_frame= res_2Dm$partner_df, 
#					   MDFS_1D_result= res_1Dm)
#   res_2Dm$partner_set=only_partners(annotated_partner_df= res_2Dm$partner_df, res_1Dm)

#   res_2D0$partner_df= get_interaction_partner(dataset=X0, decision=y,
#					       	interesting.vars= res_2D0$relevant.variables)
#   res_2D0$partner_df= find_true_partners( partner_data_frame= res_2D0$partner_df, 
#					   MDFS_1D_result= res_1D0)
#   res_2D0$partner_set=only_partners(annotated_partner_df= res_2D0$partner_df, res_1D0)

   RESULT<-list(
     res_1Dm=res_1Dm,
     res_2Dm=res_2Dm,
     res_1D0=res_1D0,
     res_2D0=res_2D0
   )
   
}
#helper function for MDFS 2D:
# given input data and result of a feature selection, 
# function finds the best fitting synergistic partner to each
# variable in interesting.vars.

get_interaction_partner<- function(dataset, decision, interesting.vars, ...){

if (!(all(interesting.vars %in% 1:ncol(dataset))))
        stop("interesting.vars must be in 1:ncol(dataset)")

if (is.null(colnames(dataset)))
        stop("datatset should have column names (taxa names)")

best.tuples <- ComputeInterestingTuplesDiscrete(dataset,
                                        decision,
                                        interesting.vars =interesting.vars,
					dimensions=2, ...
                                          )
do.call(rbind,
       lapply(interesting.vars, function(V){

               V_subset<- best.tuples$Var == V
               V_which.max<- which.max( best.tuples[ V_subset,]$IG )
              m_t1<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.1"]]
              m_t2<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.2"]]
              p_IDX<- if( V== m_t1) m_t2 else m_t1
              data.frame( rel_IDX= V, partner_IDX= p_IDX, rel_name = colnames(dataset)[[V]],                                                                             partner_name= colnames(dataset)[[p_IDX]],
                 rel_IG= best.tuples[ V_subset, ][V_which.max,][["IG"]] )
              }
             )
       )
}

#helper function to helper function:
#using output of get_interaction_partner and 1D MDFS,
# function adds additional column to the partner data frame 
# that allow one to decide if synergy is relevant
# IG_increase - how much bigger statistic for each relevant feature grows in 2D 
#               in comparison with 1D
# synergy - TRUE, if IG_increase > synergy_thr (which if left null, is by default 
#						set to minimum value of stat in 1D test
#						    of 2D relevan features)

find_true_partners<- function(partner_data_frame, MDFS_1D_result,
			      synergy_thr= NULL
					       ){
  stopifnot(all(partner_data_frame$rel_IDX %in% 1:length(MDFS_1D_result$statistic)))
  stopifnot(all(partner_data_frame$partner_IDX %in% 1:length(MDFS_1D_result$statistic)))


  stat_1D_rel2D<- MDFS_1D_result$statistic [ partner_data_frame$rel_IDX ]
  if(is.null(synergy_thr)) synergy_thr = min(stat_1D_rel2D)
  partner_data_frame$IG_increase= partner_data_frame$rel_IG - stat_1D_rel2D
  partner_data_frame$synergy= (partner_data_frame$IG_increase > synergy_thr)
  partner_data_frame
}

#helper to helper of helper...
#given partner data frame with additional columns as returned from find_true_partners()
# and output of MDFS 1D,
# this function returns a vector of feature names that were found as truly synergistic
# partners, but not relevant in 2D or 1D MDFS analysis by themselves.

only_partners<- function( annotated_partner_df, MDFS_1D_result  ) {
	stopifnot(all(c("IG_increase", "synergy") %in% colnames(annotated_partner_df)))
	annotated_partner_df$partner_name [ 
		annotated_partner_df$synergy ] |> unlist() |> as.character()
}
