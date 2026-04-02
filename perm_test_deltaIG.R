
library(magrittr)
library(matrixStats)
library(MDFS)
metadata.filename <- "metadata.csv"
taxa.filename <- "taxa.tsv"
tdata <- as.data.frame(t(read.delim(taxa.filename, header = TRUE, sep = "\t")))
colnames(tdata) <- tdata[1,]
tdata <- tdata[-1,]

metadata<-read.csv(metadata.filename, header = TRUE, sep = ",")
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]


c("sex","height.cm","weight.kg","age.years")-> confounding_factors
confounding_factors %in% colnames(metadata)
metadata$sex[!(metadata$sex %in% c("male","female"))]<- NA
metadata$height[metadata$height.cm<38]<-NA
metadata$weight[metadata$weight.kg<1]<-NA
metadata$age_years[metadata$age.years<=0]<-NA
metadata$height[metadata$height.cm>220]<-NA
metadata$weight[metadata$weight.kg<=1]<-NA
metadata$weight[metadata$weight.kg>200]<-NA
to_remove_mask <- (metadata[,confounding_factors] %>% is.na() %>% rowSums()) > 0
metadata<-metadata[!to_remove_mask, ]

tdata<- tdata[, grep('sk__Bacteria', colnames(tdata)) ]
print("inital number of taxa samples")
print(nrow(tdata))
#2. keep stool samples
metadata$env_material-> sample_source
feces_rownames<- rownames(metadata)[ sample_source=="feces" ]
tdata<-tdata[ intersect(rownames(tdata),feces_rownames),] 
print("number of taxa samples from stool")
print(nrow(tdata))

  t_rownames<- rownames(tdata)
  tdata<- lapply(tdata, as.numeric) %>% as.data.frame()
  rownames(tdata)<- t_rownames
  lib_sizes<- rowSums(tdata) 
  t_props<-( tdata
            )/lib_sizes #due to R recycling, each row gets divided by its sum
  shannon_div<- -rowSums((t_props+1e-10)*log(t_props+1e-10))
  simpson_index<- rowSums ( 1 - t_props^2  )
  inverse_simp<- 1/rowSums( t_props^2 )
  size_order<-order(lib_sizes)

  taxa_min_lib_size= 3636 
  crit_rank= which.min(abs(lib_sizes[size_order] - taxa_min_lib_size))



  tdata<- tdata[ rowSums(tdata) >= taxa_min_lib_size, ]
  tdata<- tdata/rowSums(tdata)


ids<-metadata$host.subject.id
names(ids)<-rownames(metadata)

donors<-unique(ids)

common_names<- intersect(rownames(metadata), rownames(tdata))
tdata<- tdata[common_names,]
metadata<-metadata[common_names,]
print("# of samples left after filtering by alpha diversity (with repeating donors)")
print(length(common_names))
print("# of unique donors")
print(length(unique(metadata$host.subject.id)))
source("subsample_nonrepeating_donors.R")
y<- 1- metadata$allergic.to.i.have.no.food.allergies.that.i.know.of 
names(y)<- rownames(metadata)
host_id<- metadata$host.subject.id
names(host_id)<- rownames(metadata)
subsample<-list()
n_repeats=2000
if (!file.exists("agp_subsamples2000_perm.rds")){
for (i in 1:n_repeats){
  set.seed(i)
  p_idx<-sample(1:length(y), length(y), replace=FALSE)
  y_p<-y
  y_p[]<- y[p_idx]
  print(sprintf("len(y_p)=%d,len(y)=%d,nrow(tdata)=%d",length(y_p),length(y),nrow(tdata)))
  subsample_donors(taxonomy=tdata, disease=y_p,
                   host_id=host_id,
                   subsampleSize = 0.63,
                   print_summary = FALSE
                   )-> subsample[[i]]
  subsample[[i]]$p_idx<-p_idx
if (i%% 50==0) message(sprintf("%d/%d done",i,n_repeats))
}
  invisible()
saveRDS(subsample, "agp_subsamples2000_perm.rds") 
} else {
  subsample<-readRDS("agp_subsamples2000_perm.rds")
}

source("FS_functions.R")
lop= 0.37 #leave out percentage

set.seed(123)
seeds= sample.int(234234,n_repeats)

# IG2D( Z, Y)= max_j H(Y| X_j) - H(Y| X_j, Z)
# IG1D(X_j | Z)= H( Y| Z) - H( Y| X_j,Z)
# H( Y| X_j,Z) =H( Y| Z) - IG1D(X_j | Z)

#IG1D(Z)= H(Y) - H(Y|Z)
#H(Y|Z) = H(Y) - IG1D(Z)
# so also H(Y| X_j) = H(Y) - IG1D(X_j)

#therefore
# IG2D(Z,Y)= max_j{ H(Y) - IG1D(X_j)  - H(Y|Z) + IG1D(X_j | Z) } 
# IG2D(Z,Y)= max_j{IG1D(Z) - IG1D(X_j) + IG1D(X_j | Z) }
# DELTA IG(Z):= IG2D(Z) - IG1D(Z) = max_j { IG1D(X_j | Z) - IG1D(X_j) }
DELTA_IGs<- function(X, y,
		     V_subset,
		   lvl=0.05,
		   p.adjust.method="holm",
		   seed=NULL,
		   mc=30){
   attributes(y)<-NULL
   tholds_m= colMedians(X|> as.matrix()) 
   binarize_taxa(X,tholds_m)|> filter_sparse_cols(thr=mc) -> Xm
   tholds_0= rep(0,ncol(X)) 
   binarize_taxa(X,tholds_0)|> filter_sparse_cols(thr=mc) -> X0
   names(tholds_m)<-names(tholds_0)<- colnames(X)
   common<- intersect(colnames(X0), colnames(Xm))
   X0<-X0[,common]
   Xm<-Xm[,common]
   to_keep<- c("statistic")

   MDFS.discrete(data=Xm, decision=y,dimensions = 1, 
       pc.xi = 0.25,seed=123)-> mdfs_1dm
   MDFS.discrete(data=X0, decision=y,dimensions = 1, 
       pc.xi = 0.25,seed=123)-> mdfs_1d0

   IG1D0<-mdfs_1d0$statistic
   IG1Dm<-mdfs_1dm$statistic
   delta_IG<- rep(NA, length(V_subset))
   names(delta_IG)<- V_subset
   i=1
   for (v in V_subset){
      IG1D_given_v0<- MDFS.discrete.confounded(data= X0, decision=y,
		       confounders= X0[, v], dimensions=1)$statistic
      IG1D_given_vm<- MDFS.discrete.confounded(data= Xm, decision=y,
		      confounders= Xm[, v], dimensions=1)$statistic
      delta_IG[[i]]=max(pmax(IG1D_given_v0 - IG1D0, IG1D_given_vm - IG1Dm))
      i=i+1}

   delta_IG
   
}
deltaIG_pattern="features/perm_agp_deltaIGs_set_j=%d_lop=%.2f.rds"
plausible_synergizers<- readRDS("plausible_relevant2D_and_synergistic.rds")
delta_igs<-list()
t_start<-Sys.time()
for (j in 1:2000){
  y_p<-y
  p_idx<- subsample[[j]]$p_idx #permutation index
  y_p[]<- y[p_idx] #permuted decision
  d_j_fname= sprintf(deltaIG_pattern,j,lop)
  if(!file.exists(d_j_fname)){
    message(d_j_fname)
      set.seed(seeds[[j]])
      DELTA_IGs(X=tdata[ subsample[[j]]$keep, ],
		y= y_p[ subsample[[j]]$keep ],
		V_subset=plausible_synergizers,
		   lvl=0.05,
		   p.adjust.method="holm",
		   seed=seeds[[j]],
		   mc=30) -> d_j
      saveRDS(d_j, d_j_fname)
      delta_igs[[length(delta_igs)+1]]=d_j
  
  message(sprintf("%d /%d done",j,n_repeats))
  message(Sys.time() - t_start)
  }
}
