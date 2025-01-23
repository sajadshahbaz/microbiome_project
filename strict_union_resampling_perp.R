

load("randomdataset.RData")
agp_sets<- randomdataset; rm(randomdataset)
N_TRIALS=30

keep_frac=1

n1<-  sum(agp_sets[[1]]$y_sampled==1)
n0<-  sum(agp_sets[[1]]$y_sampled==0)

n1k<- ceiling(n1*keep_frac)
n0k<- ceiling(n0*keep_frac)

set.seed(32523)

B_idx<- lapply(1:N_TRIALS,
                       function(B)
                       {
                         y<-agp_sets[[B]]$y_sampled
                         
                         idx_1= which(y==1)
                         idx_0= which(y==0)
                         c( sample(idx_1, n1k, replace=TRUE),
                        sample(idx_0, n0k, replace=TRUE)
                         )
                       }
                      )

saveRDS(B_idx, "Bootstrap_IDX.rds")
