
get_minority_class_size<- function(COL) {
                                          tab<-table(COL); l1<- length(tab)==1; (1-l1)*min(tab)
                                        }

draw_once_per_host<- function( taxonomy, disease, source_material, host_id,
                               source_type="feces",
                               taxa_to_skip_IDX=c(1:9, 1022:1028), #archaea.eucaryotes
                               taxa_min_lib_size=3636,
                               min_minority_class_taxa=20,
                               min_minority_class_dis=30,
                               plot_alpha_diversity=FALSE,
                               alpha_div_plot2pdf=FALSE,
                               print_summary=TRUE,
                               report_progess=FALSE) {
  if (!is.null(dim(disease)))
    stop("disease must be a 1-d atomic vector")
  
  if (!is.null(dim(source_material)))
    stop("source_material must be a 1-d atomic vector")
  
  if (!is.null(dim(host_id)))
    stop("host_id must be a 1-d atomic vector")
  
  if (is.null(dim(taxonomy))) 
    stop("taxonomy must have dim of length 2")
  
  if (length(dim(taxonomy))<2)
    stop("taxonomy must have dim of length 2")
  
 #names check 
  
 if (is.null(names(host_id)))
   stop(" host_id vector must have names of samples")
  
 if (is.null(names(source_material)))
   stop(" source_material vector must have names of samples")

 if (is.null(names(disease)))
   stop(" disease vector must have names of samples")
  
 if (is.null(rownames(taxonomy)))
   stop(" taxonomy must have names of samples as rownames")
  
 sources_available<-unique(source_material)
 if (!(source_type %in% sources_available))
   stop("given source_type is not a member of source_material")
 
   source_type_names<- names( source_material[ source_material== source_type ] )
   
 Reduce( intersect, list(names(host_id), names(disease),source_type_names, rownames(taxonomy)))-> common_names
  if (!length(common_names))
    stop("names of host_id, disease, taxonomy and source_material with given source_type must have nonempty intersection.")

 
   
    
 
 if (length(unique(disease[!is.na(disease)]))!=2)
   stop("disease must take values either 1 or 0")
 
 if (!all(disease[!is.na(disease)]%in% c(0,1)) )
   stop("disease must take values either 1 or 0")
 
 
 
preproc_summary=list()
   
#### prepare data types
 taxonomy<- taxonomy[,-taxa_to_skip_IDX, drop=FALSE] 
 if (ncol(taxonomy)<2)
   stop(" too many taxa_to_skip_IDX ")
 for(i in 1:ncol(taxonomy)) taxonomy[,i] <- as.numeric(taxonomy[,i])
  
preproc_summary[["# NAs in disease replaced by 0"]]=sum(is.na(disease))
 disease[is.na(disease)]<-0
 dn<- names(disease)
 disease<-as.numeric(disease)
 names(disease)<-dn

 
#### prepare taxa
preproc_summary[["initial # of taxa samples"]]=nrow(taxonomy)
st_names<- intersect(rownames(taxonomy), source_type_names)  
taxonomy<- taxonomy[ st_names, , drop=FALSE]
if (nrow(taxonomy)<max(min_minority_class_taxa,min_minority_class_dis))
  stop("too little taxa left with given source_type for given minority class minimums")
preproc_summary[["# of taxa samples with given source_type"]]=nrow(taxonomy)

preproc_summary[["minimum lib size of taxa samples"]]=taxa_min_lib_size
if (plot_alpha_diversity)
  {
  if (report_progess) print("working on alpha diversity...")
  lib_sizes<- rowSums(taxonomy) #+1e-10
  t_props<-( taxonomy#+1e-10
            )/lib_sizes #due to R recycling, each row gets divided by its sum
  shannon_div<- -rowSums((t_props+1e-10)*log(t_props+1e-10))
  simpson_index<- rowSums ( 1 - t_props^2  )
  inverse_simp<- 1/rowSums( t_props^2 )
  size_order<-order(lib_sizes)
  crit_rank= which.min(abs(lib_sizes[size_order] - taxa_min_lib_size))
  preproc_summary[["rank of min sample taxa lib size kept"]]=crit_rank
  if (alpha_div_plot2pdf)
    pdf("alpha_diversity.pdf")
  par(mfrow=c(3,1))
  plot(shannon_div[size_order],type='l')
  abline(v = crit_rank, col="red",lwd=2) 
  plot(simpson_index[size_order],type='l')
  abline(v = crit_rank, col="red",lwd=2) 
  plot(inverse_simp[size_order],type='l')
  abline(v = crit_rank, col="red",lwd=2) 
  if (alpha_div_plot2pdf)
    dev.off()
  if (report_progess) print("     ...done")
  }
taxonomy<- taxonomy[ rowSums(taxonomy) >= taxa_min_lib_size, , drop=FALSE ]
if (nrow(taxonomy)<max(min_minority_class_taxa,min_minority_class_dis))
  stop("too little taxa left with given source_type after min lib size filter for given minority class minimums")
preproc_summary[["# taxa samples kept due to minimum lib size fliter"]]=nrow(taxonomy)

taxonomy<- taxonomy/rowSums(taxonomy)

# keep only samples which have all info needed
 Reduce(intersect,list(names(host_id), names(disease),source_type_names, rownames(taxonomy)))-> common_names
 
 if (!length(common_names))
   stop(" after filtering samples by minimum library size, no samples from given source_type are left
        with info on disease and host_id")

taxonomy <- taxonomy[common_names,, drop=FALSE]
host_id<-host_id[common_names]
disease<-disease[common_names]

if (nrow(taxonomy)<max(min_minority_class_taxa,min_minority_class_dis))
  stop("too little samples left with given source_type after min lib size filter for given minority class minimums with 
       full info on host_id and disease data")
stopifnot(length(disease)==length(host_id))
stopifnot(length(host_id)==nrow(taxonomy))
preproc_summary[["# final samples available with repeats in host_id"]]=nrow(taxonomy)
length(unique(host_id))-> n_unique_hosts
preproc_summary[["# unique hosts"]]= n_unique_hosts


# check disease minority class size

dis_mc<- get_minority_class_size(disease)
print(length(disease))
print(length(names(disease)))
print(table(disease))
print(dis_mc)
if (dis_mc< min_minority_class_dis)
  stop(" after all sample filtering done, disease does not contain at least min_minority_class_dis samples")
preproc_summary[["# of disease samples with repeats in host_id"]]= sum(disease)

saveRDS(taxonomy,"taxonomy_with_repeatsANDsparseVars.rds")


# remove too sparse taxa

taxonomy<- as.data.frame(taxonomy)
t_0c<- colSums(taxonomy==0)
taxonomy<- taxonomy[, !(t_0c < min_minority_class_taxa), drop=FALSE ]
if (ncol(taxonomy)==0)
  stop("no taxa are left with given source_type and info on host_id and disease which are not too sparse")
preproc_summary[["# of taxa left with more than min_minority_class_taxa nonzeros before sampling"]]= ncol(taxonomy)

saveRDS(taxonomy, "taxonomy_with_repeats.rds")
saveRDS(host_id, "host_id_with_repeats.rds")
saveRDS(disease, "disease_vector_with_repeats.rds")

# sampling process:

## find repeating host ids
if (report_progess) print("working on sampling non repeating hosts...")
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

if (get_minority_class_size(disease)< min_minority_class_dis)
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
taxonomy<- as.matrix(taxonomy)

if (report_progess) print("working on binarization")
tb<-0*taxonomy
for(i in 1:ncol(taxonomy)) tb[,i]<-(taxonomy[,i]>median(taxonomy[,i]))

taxonomy<-as.data.frame(taxonomy)
tb<-as.data.frame(tb)

unlist(lapply(tb,get_minority_class_size))-> t_mc
tb<- tb[, t_mc >= min_minority_class_taxa, drop=FALSE ]
taxonomy<- taxonomy[, t_mc >= min_minority_class_taxa, drop=FALSE ]

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
       disease_status=disease,
       host_id=host_id,
       summary_data=preproc_summary
 )
}

