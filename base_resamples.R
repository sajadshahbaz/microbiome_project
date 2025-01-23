
source("draw_once_per_host.R")


metadata.filename <- "metadata.csv"
taxa.filename <- "taxa.tsv"


#Read microbiome data
tdata <- as.data.frame(t(read.delim(taxa.filename, header = TRUE, sep = "\t")))
colnames(tdata) <- tdata[1,]
tdata <- tdata[-1,]



library(magrittr)
metadata<-read.csv(metadata.filename, header = TRUE, sep = ",")
rownames(metadata)<-metadata[,1]
metadata<-metadata[,-1]
#metadata<-metadata[rownames(tdata),]
c("sex","height.cm","weight.kg","age.years")-> confounding_factors
confounding_factors %in% colnames(metadata)
metadata$sex[!(metadata$sex %in% c("male","female"))]<- NA
metadata$height[metadata$height.cm<38]<-NA
metadata$weight[metadata$weight.kg<1]<-NA
metadata$age_years[metadata$age.years<=0]<-NA
metadata$height[metadata$height.cm>220]<-NA
metadata$weight[metadata$weight.kg<=1]<-NA
metadata$weight[metadata$weight.kg>200]<-NA
to_remove_mask <- (metadata[,confounding_factors] %>% is.na() %>% rowSums()) > 1

metadata<-metadata[!to_remove_mask, ]
#Disease variables 
disease <- metadata[, c(6:9, 22:29, 33, 35, 36, 48:50, 54, 67:75, 98:110, 113:115,
                        134, 202:206, 212, 216, 218:229, 234:240, 264, 266, 289:293, 298, 300)]
#Binarize disease data (1- disease, 0- no disease)
disease[disease %in% c("Self-diagnosed",
                       "Diagnosed by an alternative medicine practitioner",
                       "Diagnosed by a medical professional (doctor, physician assistant)",
                       "Unspecified",
                       "Surgery only",
                       "Radiation therapy",
                       "Chemotherapy",
                       "Quite worried",
                       "A little worried",
                       "Type II diabetes",
                       "Type I diabetes",
                       "Gestational diabetes")]<-1
# disease must be 1d vector, binary but can contain NAs. NAs will be replaced by 0
y= 1- disease$allergic.to.i.have.no.food.allergies.that.i.know.of
#length(y2)
#length(y)
### TO DO: check if we use y from above instead of one from .rds below we get same resuls!
#y=readRDS("ANY_FOOD_ALLERGY_DECISION.rds")*1.
metadata$host.subject.id -> h_id
metadata$env_material-> sauce

#disease vector, host_id vector and source_information vectors all must have sample names set as its names
names(y)<- rownames(disease)
names(h_id)<-rownames(metadata)
names(sauce)<-rownames(metadata)
# (taxonomy table must also have sample names as rownames)
#generate sample from agp data wthout repeat of host id

randomdataset<-list()

for (i in 1:30){
current_seed <- i  # Use 'i' as the seed value
  set.seed(current_seed)  # Set the seed for reproducibility
draw_once_per_host(taxonomy = t_data, disease = y, source_material = sauce, host_id = h_id,
                   taxa_to_skip_IDX=c(1:9, 1022:1028), #archaea.eucaryotes,change if they are 
                   #in different places in 'taxonomy'
                   source_type ="feces" , #which source to keep
                   taxa_min_lib_size = 3636 , #minimum total unnorm. abundance count per sample to keep
                   min_minority_class_taxa = 20 , #for binarizing taxa
                   min_minority_class_dis = 30,  #for disease.
                   # if at any point of filtering, the minority class size of any taxa gets below < min_minority_class_taxa
                   # or minority class size of disease gets below < min_minority_class_dis,
                   # execution stops and informative message is printed
                   plot_alpha_diversity = (i==1),  #set true if u want a plot of alpha diversity (takes more time then)
                   alpha_div_plot2pdf = TRUE,   #set true if plot is to be saved in pdf instead of displayed
                   print_summary =TRUE, # prints sizes of the dataset and disease status distribution after each filter step
                   report_progess = TRUE
)-> agp_sample 

#TO DO: below shuffling, taking things out of agp_sample and putting them into data_list is completely 
#        unnecessary.... we just have to save agp_sample each time
#        but doing this, other scripts use randomdataset.RData, so names might not match.

binary_taxa_sampled<- agp_sample$taxa_binary
abundance_sampled<-agp_sample$taxa_abundance
y_sampled<- agp_sample$disease_status
h_id_sampled<- agp_sample$host_id # host id of each sample, for reference
agp_sample$summary_data-> preproc_summary # summary of the preprocessing process for reference
preproc_summary$seed_used <- current_seed  #useful to add it to the list
data_list <- list(
    agp_sample, binary_taxa_sampled, abundance_sampled, 
    y_sampled, h_id_sampled, preproc_summary, preproc_summary$seed_used
  )
  
  # Correctly name the datasets
  datasetnames <- c(
    "agp_sample", "binary_taxa_sampled", "abundance_sampled", 
    "y_sampled", "h_id_sampled", "preproc_summary", "seed_used"
  )
  
  names(data_list) <- datasetnames

randomdataset[[i]]<-data_list
}
save(randomdataset, file = "randomdataset.RData")
