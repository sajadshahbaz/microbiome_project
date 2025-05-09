

draw_once_per_host<- function( taxonomy, disease,host_id,
                               min_minority_class_taxa=20,
                               min_minority_class_dis=30,
                               print_summary=TRUE,
                               report_progess=FALSE,
			       binarize_on_median=TRUE	) {
  if (!is.null(dim(disease)))
    stop("disease must be a 1-d atomic vector")
  
  if (!is.null(dim(host_id)))
    stop("host_id must be a 1-d atomic vector")
  
  if (is.null(dim(taxonomy))) 
    stop("taxonomy must have dim of length 2")
  
  if (length(dim(taxonomy))<2)
    stop("taxonomy must have dim of length 2")
  
 #names check 
  
 if (is.null(names(host_id)))
   stop(" host_id vector must have names of samples")
  

 if (is.null(names(disease)))
   stop(" disease vector must have names of samples")
  
 if (is.null(rownames(taxonomy)))
   stop(" taxonomy must have names of samples as rownames")
  

common_names<- Reduce(intersect, list(rownames(taxonomy), names(disease), names(host_id)))
 
taxonomy<- taxonomy[common_names,]
disease<-disease[common_names]
host_id<- host_id[common_names]
 
preproc_summary=list()

# sampling process:

## find repeating host ids
if (report_progess) print("working on sampling non repeating hosts...")
#sampling_start<-Sys.time()
table(host_id)-> hid_counts
names(hid_counts[hid_counts>1])-> hosts_w_repeats

## which samples come from hosts with more than one sample?
repeats_id_mask<- host_id %in% hosts_w_repeats
preproc_summary[["# of hosts with more than one sample"]]= length(hosts_w_repeats)
preproc_summary[["# of samples from hosts with more than one sample"]]= sum(repeats_id_mask)

# put aside data to sample from:
taxonomy_rep_hosts<- taxonomy[  repeats_id_mask,,drop=FALSE]
disease_rep_hosts<- disease[  repeats_id_mask]
host_id_rep_hosts<- host_id[  repeats_id_mask]


# prepare one sample from each repeating host_id
unlist(lapply(hosts_w_repeats, function(ID)
{
   samples_from_host<-  which(host_id_rep_hosts == ID)
   sample(samples_from_host,size = 1)
}))-> sample_IDX_per_host
#print("time for sampling")
#print(Sys.time()-sampling_start)
stopifnot(length(unique(sample_IDX_per_host))==length(sample_IDX_per_host))

## concatenate data without repeats with sampled data

taxonomy<- rbind( taxonomy[!repeats_id_mask,,drop=FALSE],
                  taxonomy_rep_hosts[sample_IDX_per_host,,drop=FALSE]
)
disease<- c( disease[!repeats_id_mask],
             disease_rep_hosts[sample_IDX_per_host])
host_id<- c(host_id[!repeats_id_mask],
            host_id_rep_hosts[sample_IDX_per_host])
stopifnot(length(host_id)==nrow(taxonomy))
stopifnot(length(host_id)==length(disease))

if  ( min(sum(disease==1),sum(disease==0))< min_minority_class_dis)
  stop(" after sampling, minority class size of disease is too small")
preproc_summary[["# of disease samples after sampling"]]= sum(disease)


## filter out too sparse taxa again
t_0c<- colSums(taxonomy==0)
taxonomy<- taxonomy[, !(t_0c < min_minority_class_taxa), drop=FALSE ]
if (ncol(taxonomy)==0)
  stop("after sampling, no taxa are left with given source_type 
       and info on host_id and disease which are not too sparse")
preproc_summary[["# of taxa left with more than min_minority_class_taxa nonzeros after sampling"]]= ncol(taxonomy)

if (report_progess) print("   ...done")
#bin_start<-Sys.time()
taxonomy<- as.matrix(taxonomy)

if (report_progess) print("working on binarization")
meds= if (binarize_on_median) colMedians(taxonomy) else rep(0,ncol(taxonomy))
tb<- t( t(taxonomy)>meds )*1. #exploit recycling mechanism.
			      #column [[i]] of `taxonomy` gets binarized by meds[[i]]
#tb2<-tb
#print("time for binarization:")
#print(Sys.time()-bin_start)
#for (i in 1:ncol(taxonomy)) tb2[,i]= (taxonomy[,i]>meds[[i]])*1.
#tb2<-as.data.frame(tb2)

#minority_time<- Sys.time()
t_mc<- pmin(colSums(tb==0),colSums(tb==1))
#print("time for minority class calc:")
#print(Sys.time()- minority_time)
tb<- tb[, t_mc >= min_minority_class_taxa, drop=FALSE ]
#tb2<-tb2[, t_mc >= min_minority_class_taxa, drop=FALSE ]
taxonomy<- taxonomy[, t_mc >= min_minority_class_taxa, drop=FALSE ]

taxonomy<-as.data.frame(taxonomy)
tb<-as.data.frame(tb)
preproc_summary[["# of binary taxa left with sufficient minority class size "]]= ncol(taxonomy)
if (ncol(tb)==0)
  stop("after sampling and median binarization
        no taxa are left with sufficient minority class size")
preproc_summary[[" final sample size"]]= nrow(taxonomy)

if (report_progess) print("    ...done")
if (print_summary)
  print(preproc_summary)
 list( taxa_abundance= taxonomy,
       taxa_binary= tb,
#	tb2=tb2,
       disease_status=disease,
       host_id=host_id,
       summary_data=preproc_summary
 )
}

