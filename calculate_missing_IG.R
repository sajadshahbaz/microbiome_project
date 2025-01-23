
args=commandArgs(trailingOnly=TRUE)
i<- as.integer(args[[1]])

table_to_fill<- read.delim("strict_union_base_result_missingIG.csv", sep=",")
table_to_fill<- table_to_fill[,2:ncol(table_to_fill)]

table_to_fill$name[ table_to_fill[,"IG_1D"]==-999 ] -> missing_1D
table_to_fill$name[ table_to_fill[,"IG_2D"]==-999 ] -> missing_2D

load("randomdataset.RData"); agp_sets<- randomdataset; rm(randomdataset)

library(MDFS)

ComputeMaxInfoGains(data=agp_sets[[i]]$binary_taxa_sampled[, missing_1D], 
		decision= agp_sets[[i]]$y_sampled,dimensions=1, divisions=1, discretizations=1, range=0)-> IG1D


ComputeMaxInfoGains(data=agp_sets[[i]]$binary_taxa_sampled,
		decision= agp_sets[[i]]$y_sampled,dimensions=2, divisions=1, discretizations=1, range=0)-> IG2D
rownames(IG2D)<- colnames(agp_sets[[i]]$binary_taxa_sampled)
IG2D<- IG2D[ missing_2D,, drop=FALSE]
print(dim(IG2D))
print(length(IG2D))
rownames(IG2D)<-NULL
print(IG2D)


IG1D$names<- missing_1D
IG2D$names<- missing_2D

saveRDS(list(IG1D, IG2D ), 
	sprintf("missing_IG_agp%d.rds",i))

#IG1Dall$names<- colnames(agp_sets[[i]]$binary_taxa_sampled)
#IG2Dall$names<- colnames(agp_sets[[i]]$binary_taxa_sampled)

#saveRDS(IG1Dall, "allIG1Dagain.rds")
#saveRDS(IG2Dall, "allIG2Dagain.rds")
