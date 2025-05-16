

subsample_donors<- function( taxonomy, #all samples, including repeating donors
			      	disease,host_id,  #1d vectors, must have names equal to rownames(taxonomy)
				subsampleSize=0.9, # in terms of fraction
                               min_minority_class_dis=30,
                               print_summary=TRUE
			       ) {

   stopifnot( (subsampleSize<1) && (subsampleSize>0) )
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

stopifnot(length(disease)==nrow(taxonomy))
stopifnot(length(disease)==length(host_id))
stopifnot(all(names(disease)==rownames(taxonomy)))
stopifnot(all(names(disease)==names(host_id)))

 
preproc_summary=list()
preproc_summary[["# of high quality samples, including repeating donors:"]]=nrow(taxonomy)

# sampling process:

# shuffle rows
permutation<-sample(1:nrow(taxonomy), nrow(taxonomy), replace=FALSE)

sample_IDs<- rownames(taxonomy)
sample_IDs<- sample_IDs[ permutation ]

disease<-disease[sample_IDs]
host_id<- host_id[sample_IDs]

#remove duplicates

dups<- duplicated(host_id)
preproc_summary[["# of donors with more than one sample:"]]=length(unique(host_id[ dups ]))
disease<-disease[!dups]
host_id<- host_id[!dups]
preproc_summary[["# of unique donors:"]]= length(host_id)

n_keep=floor(subsampleSize*length(host_id))
keep_idx= 1: n_keep
sample_IDs<- names(host_id)
stopifnot(length(sample_IDs)==length(disease))
stopifnot(all(names(disease)==sample_IDs))
keep_sampleIDs<-sample_IDs[keep_idx]
leaveOut_sampleIDs<- sample_IDs[-keep_idx]



if  ( min(sum(disease[keep_sampleIDs]==1),sum(disease[keep_sampleIDs]==0))< min_minority_class_dis)
  stop(" after sampling, minority class size of disease is too small")


if (print_summary)
  print(preproc_summary)
 list( keep=keep_sampleIDs,
       leave=leaveOut_sampleIDs,
       summary=preproc_summary)
}

