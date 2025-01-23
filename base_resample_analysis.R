setwd("agp_projects/source_vs_disease_interactions/")
# union of all 30 runs of FS and mean IG
load("rel_MDFS_group1.rdata")

Reduce(union, rel_group[[3]])-> union_
library(magrittr)
# TO DO: this is not used right now i think.
source("useful_set_algebra.R")

Reduce("+",lapply(rel_group[[3]], function(SET)
 union_ %in% SET) 
      )-> rel_freq
Reduce("+",lapply(rel_group[[1]], function(SET)
  union_ %in% SET) 
)-> rel_freq_1D
Reduce("+",lapply(rel_group[[2]], function(SET)
  union_ %in% SET) 
)-> rel_freq_2D


names(rel_freq_1D)<-names(rel_freq_2D)<-names(rel_freq)<- union_

meanIG_per_dim<-list()


for (dimension in 1:2){
  dim_fs<-MDFS_group[[dimension]]
  lapply(dim_fs, function(run) 
		run$statistic[ run$relevant.variables ]
	)-> IG_of_rel_per_run 
  lapply(seq_along(dim_fs), function(i_run)
		rel_group[[dimension]][[i_run]]) -> taxa_per_run


  meanIG_per_dim[[dimension]]<- rep(0, length(union_))
  names(meanIG_per_dim[[dimension]])<-union_
  
  for (i_run in seq_along(IG_of_rel_per_run))
	meanIG_per_dim[[dimension]][ taxa_per_run[[i_run]] ]= meanIG_per_dim[[dimension]][ taxa_per_run[[i_run]] ]   +
							  IG_of_rel_per_run[[i_run]] 

  meanIG_per_dim[[dimension]]= meanIG_per_dim[[dimension]] / rel_freq 
  
}
# 30 runs of sampling from agp data once per host:
load("randomdataset.RData")
agp_sets<- randomdataset; rm(randomdataset)

library(MDFS)

# extracts from one run
get_interaction_partner<- function(dataset, decision, interesting.vars, ...){

if (!(all(interesting.vars %in% 1:ncol(dataset))))
	stop("interesting.vars must be in 1:ncol(dataset)")

if (is.null(colnames(dataset))) 
	stop("datatset should have column names (taxa names)")

best.tuples <- ComputeInterestingTuples(dataset,
                                        decision, 
                                        interesting.vars =interesting.vars,
					...
                                          )
do.call(rbind,
       lapply(interesting.vars, function(V){

	       V_subset<- best.tuples$Var == V
	       V_which.max<- which.max( best.tuples[ V_subset,]$IG )
	      m_t1<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.1"]]
	      m_t2<-best.tuples[ V_subset, ][V_which.max, ][["Tuple.2"]]
	      p_IDX<- if( V== m_t1) m_t2 else m_t1
	      data.frame( rel_IDX= V, partner_IDX= p_IDX, rel_name = colnames(dataset)[[V]],
		 partner_name= colnames(dataset)[[p_IDX]],
		 rel_IG= best.tuples[ V_subset, ][V_which.max,][["IG"]] )
	      }
	     )
       )
}

partners_per_run<- list()

for (i_set in seq_along(agp_sets))
{
partners_per_run[[i_set]] <- 
get_interaction_partner(dataset=agp_sets[[i_set]]$binary_taxa_sampled,
			 decision=agp_sets[[i_set]]$y_sampled,
			 interesting.vars=MDFS_group[[2]][[i_set]]$relevant.variables , 
			range=0, dimensions=2)
}


### find most frequent partner of each taxa according to agp sets

do.call(rbind, partners_per_run)-> partners_all_runs
aggregate(partners_all_runs$rel_IG,
	  list(rel=partners_all_runs$rel_name,
	      partner=partners_all_runs$partner_name),
	  mean)-> partner_IG_tabulation

per_taxa_partner_scores<- as.list(rep(NA, length(union_)))
lapply(union_, function(taxon)
{
	taxon_subset<- partner_IG_tabulation$rel == taxon
	if (any(taxon_subset))
	{
		taxon_partners<-partner_IG_tabulation[ taxon_subset, ,drop=FALSE]
	 	taxon_partners$partner[[ which.max(taxon_partners$x) ]]	
	} else NA

}	) %>% unlist() -> most_freq_partners 

names(most_freq_partners)<- union_

rel_and_partners<- unname(union(union_, most_freq_partners[!is.na(most_freq_partners)]))


# healthy and disease status stats
# and abundance/presence interactions inside aech group & total
source("pairwise_t_test.R")



analyze_presence_abundance<- function( PRESENCE,
					 ABUNDANCE,
					 Y,
					 IG_per_dim,
					 AB_multiplier=10000,
					 lvl=0.01,
					 report_progress=TRUE){

pairwise_means_sderrs(abundance = ABUNDANCE, presence = PRESENCE, 
                      disease = Y, report_progress = report_progress) -> interaction_components

difference_in_differences(interaction_components, lvl=lvl) -> interaction_matrices

get_sign_matrices(interaction_matrices) -> sign_matrices

summarize_abp_IG_HvsD(abundance = ABUNDANCE*AB_multiplier, presence = PRESENCE,disease = Y,
                      IG_per_dim =IG_per_dim )-> abp_summary

#summary_signs_ab_freq<-summarize_HvsD(sign_matrices = sign_matrices,abundance = ABUNDANCE*AB_multiplier,presence = PRESENCE,disease = Y,
#                                      IG_per_dim = IG_per_dim) 

list(means_SDerrs= interaction_components,
     Z_matrices= interaction_matrices,
     sign_matrices= sign_matrices,
     abp_IG_summary=abp_summary)
    
}

PA_tables<- list()

for (i_set in seq_along(agp_sets)){

P= as.data.frame(matrix( nrow= nrow(agp_sets[[i_set]]$binary_taxa_sampled),
           ncol= length(rel_and_partners)))
A= as.data.frame(matrix( nrow= nrow(agp_sets[[i_set]]$binary_taxa_sampled),
           ncol= length(rel_and_partners)))
print(dim(P))
print(dim(A))
colnames(P)<-colnames(A)<- rel_and_partners
P[,]<-A[,]<-0
present_in_set_i<- intersect(colnames(agp_sets[[i_set]]$binary_taxa_sampled),
                             rel_and_partners)
for (taxon in present_in_set_i)
{P[,taxon]=agp_sets[[i_set]]$binary_taxa_sampled[,taxon]
 A[,taxon]=agp_sets[[i_set]]$abundance_sampled[,taxon]
} 
  
print(dim(P))
print(dim(A))
PA_tables[[i_set]]<- analyze_presence_abundance(PRESENCE= P,
					        ABUNDANCE= A,
Y= agp_sets[[i_set]]$y_sampled,
IG_per_dim= meanIG_per_dim
)
PA_tables[[i_set]]$P=P
PA_tables[[i_set]]$A=A
}

array_mean<- function(array_list, columns2exclude){
  
  mean_array<- array_list[[1]]
  num_cols<- lapply(array_list, function(x) x[,!(colnames(x) %in% columns2exclude), drop=FALSE ] )
  Reduce("+", num_cols)/ length(array_list) -> mean_num_cols
  mean_array[,!( colnames(mean_array) %in% columns2exclude)] = mean_num_cols
  mean_array
}
lapply(PA_tables, function(x) x$abp_IG_summary)-> summaries
array_mean(summaries, c("name","rel_1D","rel_2D","IG_1D","IG_2D"))-> mean_PA
mean_PA[  ,c("rel_2D","rel_1D")] = 0
mean_PA[ match(union_, mean_PA$name)    ,"rel_1D" ] = rel_freq_1D 
mean_PA[ match(union_, mean_PA$name)    ,"rel_2D" ] = rel_freq_2D 
mean_PA$partner_name= lapply(1:nrow(mean_PA),
                             function(i)
                             {
                               if (mean_PA$rel_2D[[i]]>0)
                                 most_freq_partners[[ mean_PA$name[[i]] ]]
                               else NA
                             }
) %>% unlist()

#TO DO
# programatically do the filtering that we have done in xls to produce list of names below:


left_names<-"sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__Neisseria
sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae
sk__Bacteria;k__;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Propionibacteriales;f__Propionibacteriaceae;g__Cutibacterium
sk__Bacteria;k__;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;g__Chryseobacterium
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Corynebacteriales;f__;g__Lawsonella
sk__Bacteria;k__;p__Firmicutes;c__Tissierellia;o__Tissierellales;f__Peptoniphilaceae
sk__Bacteria;k__;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;g__Bergeyella
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Xanthomonas
sk__Bacteria;k__;p__candidate_division_CPR1
sk__Bacteria;k__;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium
sk__Bacteria;k__;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Aliivibrio
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae
sk__Bacteria;k__;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella
sk__Bacteria;k__;p__Bacteroidetes;c__Chitinophagia;o__Chitinophagales;f__Chitinophagaceae
sk__Bacteria;k__;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;g__Capnocytophaga
sk__Bacteria;k__;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingomonas
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae
sk__Bacteria;k__;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Dolosigranulum
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Citrobacter
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Propionibacteriales;f__Nocardioidaceae;g__Nocardioides
sk__Bacteria;k__;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Leptotrichia
sk__Bacteria;k__;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae
sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__;g__Paucibacter
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales
sk__Bacteria;k__;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Robinsoniella
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales
sk__Bacteria;k__;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Carnobacterium
sk__Bacteria;k__;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Peptostreptococcaceae;g__Peptostreptococcus
sk__Bacteria;k__;p__Firmicutes;c__Bacilli;o__Bacillales;f__;g__Gemella
sk__Bacteria;k__;p__Verrucomicrobia;c__Spartobacteria;o__Chthoniobacterales;f__Chthoniobacteraceae;g__Candidatus_Udaeobacter
sk__Bacteria;k__;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae
sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae;g__Rothia
sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Methylophilaceae;g__Methylobacillus
sk__Bacteria;k__;p__Tenericutes;c__Mollicutes;o__Mycoplasmatales;f__Mycoplasmataceae;g__Mycoplasma
sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Rhodocyclales;f__Azonexaceae
sk__Bacteria;k__;p__Actinobacteria;c__Rubrobacteria;o__Gaiellales
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas
sk__Bacteria;k__;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Oscillospira
sk__Bacteria;k__;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Holdemanella
sk__Bacteria;k__;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__gamma_proteobacterium_symbiont_of_Gonopsis_affinis
"

strsplit(left_names, "\n") %>% unlist() -> left_names

summary_reduced<- mean_PA[ match(left_names, mean_PA$name), ]

View(summary_reduced)

summary_reduced$IG_1D[ summary_reduced$rel_1D < 30 ] = -999
summary_reduced$IG_2D[ summary_reduced$rel_2D < 30 ] = -999


write.csv(summary_reduced, "strict_union_base_result_missingIG.csv")

save.image("after_base_resapling_analysis.rda")
