
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
n_repeats=500
if (!file.exists("agp_subsamples63_perm.rds")){
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
if (i%% 50==0) message(sprintf("%d/500 done",i))
}
  invisible()
saveRDS(subsample, "agp_subsamples63_perm.rds") 
} else {
  subsample<-readRDS("agp_subsamples63_perm.rds")
}

source("FS_functions.R")
features_file_string="features/perm_agp_features_selected_lop=%.2f.rds"
lop= 0.37 #leave out percentage

#if feature selection was not performed yet, set up the lists
fs_fname=sprintf(features_file_string,lop)
if (!file.exists(fs_fname)) {
features<- list(md1D=list(),
                md12D=list(),
                md12Dp=list(),
                u=list()
                )

} else { # otherwise load the saved result
  features<-readRDS(fs_fname)
}

set.seed(123)
seeds= sample.int(234234,n_repeats)

fs_j_pattern="features/perm_agp_fs_set_j=%d_lop=%.2f.rds"
 full_MDFSres<- full_Ures<-list()
for (j in 1:n_repeats){
  y_p<-y
  p_idx<- subsample[[j]]$p_idx #permutation index
  y_p[]<- y[p_idx] #permuted decision
  fs_j_fname= sprintf(fs_j_pattern,j,lop)
  if(!file.exists(fs_j_fname)){
    message(fs_j_fname)
      set.seed(seeds[[j]])
      MDFS_FS(X= tdata[ subsample[[j]]$keep, ], 
              y= y_p[ subsample[[j]]$keep ], #note permuted decision
	      lvl=0.05, p.adjust.method = "holm",
              seed= seeds[[j]],
              mc=30 #minimum allowed class size of binarized taxa abundance
              )-> all_mdfs_j
      U_test_FS(X=tdata[ subsample[[j]]$keep, ],
                y = y_p[ subsample[[j]]$keep ], #permuted.
                lvl = 0.05,
                p.adjust.method = "holm",
                min_presence=30)->u_j
      fs_j=list(mdfs=all_mdfs_j, u=u_j)
      saveRDS(fs_j, fs_j_fname)
      } else {
        fs_j=readRDS(fs_j_fname)
        all_mdfs_j=fs_j$mdfs
        u_j=fs_j$u
  }
  full_MDFSres[[j]]=all_mdfs_j
  full_Ures[[j]]=u_j
  if  (!file.exists(features_file_string)) {
            features$u[[j]]=u_j$rel_set
            features$md1D[[j]]=union(all_mdfs_j$res_1Dm$rel_set, all_mdfs_j$res_1D0$rel_set)
            features$md12D[[j]]=Reduce(union, list(all_mdfs_j$res_2Dm$rel_set, 
                                                  all_mdfs_j$res_2D0$rel_set,
                                                  features$md1D[[j]]))
            features$md12Dp[[j]]=Reduce(union, list(all_mdfs_j$res_2Dm$partner_set, 
                                                  all_mdfs_j$res_2D0$partner_set,
                                                  features$md12D[[j]]))
          }
  
  message(sprintf("%d /%d done",j,n_repeats))
  
}
