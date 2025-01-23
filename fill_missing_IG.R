
table_to_fill<- read.delim("strict_union_base_result_missingIG.csv", sep=",")
table_to_fill<- table_to_fill[,2:ncol(table_to_fill)]

table_to_fill$name[ table_to_fill[,"IG_1D"]==-999 ] -> missing_1D
table_to_fill$name[ table_to_fill[,"IG_2D"]==-999 ] -> missing_2D

ig1_missing<- rep(0, length(missing_1D))
ig2_missing<- rep(0, length(missing_2D))
for (i in 1:30){

	ith_ig<-readRDS(sprintf("missing_IG_agp%d.rds",i))
stopifnot(length(ith_ig[[1]]$IG)==length(ig1_missing))
stopifnot(length(ith_ig[[2]]$IG)==length(ig2_missing))
ig1_missing= ig1_missing + ith_ig[[1]]$IG
ig2_missing= ig2_missing + ith_ig[[2]]$IG

}
ig1_missing= ig1_missing/30
ig2_missing= ig2_missing/30

rownames(table_to_fill)<- table_to_fill$name
table_to_fill[ missing_1D ,"IG_1D"] = ig1_missing
table_to_fill[ missing_2D ,"IG_2D"] = ig2_missing

write.csv(table_to_fill, "filled_IG_table.csv")
