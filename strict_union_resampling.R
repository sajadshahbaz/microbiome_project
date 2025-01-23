N_TRIALS=30
args=commandArgs(trailingOnly = TRUE)
B<-as.integer(args[[1]])
stopifnot(B %in% 1:N_TRIALS)
load("randomdataset.RData"); set_B<- randomdataset[[B]]; rm(randomdataset)
readRDS("Bootstrap_IDX.rds")-> bs_IDX_list; IDX_B<- bs_IDX_list[[B]]
set.seed(32432)
session_seeds<- sample.int(323423,N_TRIALS)
set.seed(session_seeds[[B]])
#FS
library(MDFS)
  MDFS_1D<- MDFS(
    data = set_B[["binary_taxa_sampled"]][IDX_B,],
    decision = set_B[["y_sampled"]][IDX_B],
    dimensions = 1,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed= session_seeds[[B]], 
    level = 0.05
  )
  rel_1D<- colnames(set_B[["binary_taxa_sampled"]])[MDFS_1D$relevant.variables]
  # Perform MDFS with 2D
  MDFS_2D<- MDFS(
    data = set_B[["binary_taxa_sampled"]][IDX_B,],
    decision = set_B[["y_sampled"]][IDX_B],
    dimensions = 2,
    divisions = 1,
    discretizations = 1,
    range = 0,
    seed = session_seeds[[B]],
    level = 0.05
  )
  rel_2D<- colnames(set_B[["binary_taxa_sampled"]])[MDFS_2D$relevant.variables]
  rel_1or2D<-  union( rel_1D, rel_2D)  
  FS<- list(rel_1D= rel_1D,
            rel_2D= rel_2D,
            rel_1or2D=rel_1or2D)
  saveRDS( FS, sprintf("MDFS_results/FS_replicate_%d.rds", B ) )
  saveRDS( list(MDFS_1D=MDFS_1D,
                MDFS_2D=MDFS_2D),
            sprintf("MDFS_results/MDFS_replicate_%d.rds", B ) )
#RF
  library(randomForest)
  #library(pROC)
read.delim("filled_IG_table.csv", sep=",")-> ref_table
ref_table$name[ ref_table$rel_1D==30 ] -> consensus_1D
ref_table$name -> consensus_AllRelevant
unique(IDX_B)-> uqIDX_B
set.seed(session_seeds[[B]])
randomForest(x = set_B[["binary_taxa_sampled"]][uqIDX_B, consensus_1D],
             y =  as.factor(set_B[["y_sampled"]][uqIDX_B]),
             importance = TRUE,
             keep.forest = TRUE)-> RF_1D
randomForest(x = set_B[["binary_taxa_sampled"]][uqIDX_B,consensus_AllRelevant],
             y =  as.factor(set_B[["y_sampled"]][uqIDX_B]),
             importance = TRUE,
             keep.forest = TRUE)-> RF_all
saveRDS(list(RF_1D=RF_1D,
             RF_all=RF_all),
        sprintf("RF_models/RF_replicate_%d.rds", B)
)



