taxonomy_with_repeatsANDsparseVars <- readRDS("D:/OwnCloud/Microbiome/alergia/taxonomy_with_repeatsANDsparseVars.rds")
meta_rds <- readRDS("D:/OwnCloud/Microbiome/alergia/meta_rds.rds")
disease_vector_with_repeats <- readRDS("D:/OwnCloud/Microbiome/alergia/disease_vector_with_repeats.rds") 
library(matrixStats)

#This is not very strict!
meta_rds$sex[!(meta_rds$sex %in% c("male","female"))]<- NA
meta_rds$height[meta_rds$height.cm<38]<-NA
meta_rds$weight[meta_rds$weight.kg<1]<-NA
meta_rds$age_years[meta_rds$age.years<=0]<-NA
meta_rds$height[meta_rds$height.cm>220]<-NA
meta_rds$weight[meta_rds$weight.kg<=1]<-NA
meta_rds$weight[meta_rds$weight.kg>200]<-NA

ids<-meta_rds$host.subject.id
names(ids)<-rownames(meta_rds)
ids<-ids[names(ids) %in% rownames(taxonomy_with_repeatsANDsparseVars)]

donors<-unique(ids)

taxonomy.mean<-matrix(NA,length(donors),ncol(taxonomy_with_repeatsANDsparseVars))
rownames(taxonomy.mean)<-donors
colnames(taxonomy.mean)<-colnames(taxonomy_with_repeatsANDsparseVars)
meta.mean<-as.data.frame(matrix(NA,length(donors),6))
rownames(meta.mean)<-donors
colnames(meta.mean)<-c('age','sex','weight.kg','height.cm','BMI_cat','age_cat')
decision.mean<-numeric(length(donors))
names(decision.mean)<-donors

for (i in 1:length(donors)) {
 donor<-donors[i]
 samples<-names(ids)[ids==donor]
 
 #decision
 decs<-disease_vector_with_repeats[samples]
 decision.mean[i]<-floor(median(decs))
 
 #taxa
 taxa<-taxonomy_with_repeatsANDsparseVars[samples,,drop=F]
 taxonomy.mean[i,]<-colMedians(as.matrix(taxa))
 
 #metadata
 #age
 ages<-meta_rds[samples,'age.years']
 ages<-ages[!is.na(ages)]
 if (length(ages)>0) meta.mean[i,1]<-round(median(ages))
 #sex
 sexes<-meta_rds[samples,'sex']
 sexes<-sexes[!is.na(sexes)]
 if (length(sexes)>0) meta.mean[i,2]<-names(table(sexes))[which.max(table(sexes))]
 #weight
 weights<-as.numeric(meta_rds[samples,'weight.kg'])
 weights<-weights[!is.na(weights)]
 if (length(weights)>0) meta.mean[i,3]<-median(weights)
 #height
 heights<-as.numeric(meta_rds[samples,'height.cm'])
 heights<-heights[!is.na(heights)]
 if (length(heights)>0) meta.mean[i,4]<-median(heights)
 
 
 if(!(i%%10)) print(i)
}

mask.meta<-rowSums(is.na(meta.mean[,1:4]))==0

meta.mean.final<-meta.mean[mask.meta,]
meta.mean.final[meta.mean.final[1]<15,6]<-'0-14'
meta.mean.final[meta.mean.final[1]>=15 & meta.mean.final[1]<64,6]<-'15-64'
meta.mean.final[meta.mean.final[1]>=64,6]<-'65+'
bmi<-meta.mean.final[,'weight.kg']/meta.mean.final[,'height.cm']^2*10000
meta.mean.final[bmi<18.5,5]<-'Underweight'
meta.mean.final[bmi>=18.5 & bmi<25,5]<-'Normal'
meta.mean.final[bmi>=25 & bmi<30,5]<-'Overweight'
meta.mean.final[bmi>=30,5]<-'Obese'


decision.mean.final<-decision.mean[mask.meta]
taxonomy.mean.final<-taxonomy.mean[mask.meta,]
taxonomy.mean.final<-taxonomy.mean.final[,colSums(taxonomy.mean.final>0)>20]
taxonomy.mean.final<-as.data.frame(taxonomy.mean.final/rowSums(taxonomy.mean.final)*10000)

binary.mean.final<-taxonomy.mean.final*0
for (i in 1:ncol(binary.mean.final)) binary.mean.final[,i]<-taxonomy.mean.final[,i]>median(taxonomy.mean.final[,i])

save(meta.mean.final,decision.mean.final,taxonomy.mean.final,binary.mean.final,file='mean_data.RData')


#Table 1 
table(decision.mean)

table(meta.mean.final[decision.mean==0,2])
table(meta.mean.final[decision.mean==1,2])

fisher.test(rbind(table(meta.mean.final[decision.mean==0,2]),
            table(meta.mean.final[decision.mean==1,2])))
            

table(meta.mean.final[decision.mean==0,6])
table(meta.mean.final[decision.mean==1,6])

fisher.test(rbind(table(meta.mean.final[decision.mean==0,6]),
                  table(meta.mean.final[decision.mean==1,6])))#,simulate.p.value = T)


table(meta.mean.final[decision.mean==0,5])
table(meta.mean.final[decision.mean==1,5])

fisher.test(rbind(table(meta.mean.final[decision.mean==0,5]),
                  table(meta.mean.final[decision.mean==1,5])))

