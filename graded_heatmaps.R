
#load results of base analysis
load("after_base_resampling_analysis.rda")

library(matrixStats)

reduce_sign_matrices_to_name_subset<- function(PA_tabs, chosen_names){
  
reduced_sign_matr<- lapply(PA_tabs, function (x) x$sign_matrices)
for (i in 1:length(reduced_sign_matr))
{
  colnames(reduced_sign_matr[[i]]$total)<-rownames(reduced_sign_matr[[i]]$total)<- reduced_sign_matr[[i]]$taxa_names
  colnames(reduced_sign_matr[[i]]$d)<-rownames(reduced_sign_matr[[i]]$d)<- reduced_sign_matr[[i]]$taxa_names
  colnames(reduced_sign_matr[[i]]$h)<-rownames(reduced_sign_matr[[i]]$h)<- reduced_sign_matr[[i]]$taxa_names
  reduced_sign_matr[[i]]$total<- reduced_sign_matr[[i]]$total[ chosen_names,
                                                               chosen_names]
  reduced_sign_matr[[i]]$d<- reduced_sign_matr[[i]]$d[ chosen_names,
                                                               chosen_names]
  reduced_sign_matr[[i]]$h<- reduced_sign_matr[[i]]$h[ chosen_names,
                                                               chosen_names]
}
reduced_sign_matr
  
}
rd_signs<- reduce_sign_matrices_to_name_subset(PA_tabs = PA_tables, chosen_names = left_names)


consensus_heatmap<- function(matrices_list, type="total",min_freq=30){
sign_v<-do.call(rbind, lapply(matrices_list, function(x) 
                  as.vector(x[[type]] )
                ))
 n_p=colSums(sign_v==1)  
 n_n=colSums(sign_v==-1)   
 n_0=-n_p -n_n + nrow(sign_v)
 print(max(n_0))
 rbind(n_p,n_n,n_0)-> sign_counts
 colMaxs(sign_counts)->n_mostfreq
 print(max(n_mostfreq))
 apply(sign_counts,2,function(COL)  c(1,-1,0)[[which.max(COL)]])-> consensus_signs
 consensus_signs[ n_mostfreq < 30 ] = 0
 map_<- matrix(nrow=nrow( matrices_list[[1]][[type]]),
               ncol=ncol(matrices_list[[1]][[type]])
              )
 map_[,]= consensus_signs
 rownames(map_)<-colnames(map_)<- rownames( matrices_list[[1]][[type]])
 map_
}


interaction_heatmap<- function(map_matrix, main="interaction heatmap"){
  rownames(map_matrix)<-colnames(map_matrix)<-NULL
  heatmap(map_matrix, Rowv=NA,  Colv=NA, scale="none", revC=TRUE, col=c("red","white","blue"), main=main)
}
map_tc<- consensus_heatmap(rd_signs, type = "total", min_freq = 30)
table(map_tc)
interaction_heatmap(map_tc)

map_hc<- consensus_heatmap(rd_signs, type = "h", min_freq = 30)
table(map_hc)
interaction_heatmap(map_hc)

map_dc<- consensus_heatmap(rd_signs, type = "d", min_freq = 30)
table(map_dc)
interaction_heatmap(map_dc)
#not the final plots (base ones)
#pdf("abundance_heatmap_total_finalOrder_base_resamples.pdf")
#interaction_heatmap(map_tc)
#dev.off()
#
#
#pdf("abundance_heatmap_dis_finalOrder_base_resamples.pdf")
#interaction_heatmap(map_dc)
#dev.off()
#
#
#pdf("abundance_heatmap_hlth_finalOrder_base_resamples.pdf")
#interaction_heatmap(map_hc)
#dev.off()

stopifnot( all(rownames(map_dc)==summary_reduced$name))
stopifnot( all(rownames(map_hc)==summary_reduced$name))
stopifnot( all(rownames(map_tc)==summary_reduced$name))

#compare with conensus over resampling
read.delim("filled_IG_table.csv", sep=",")-> ref_table
stopifnot( all(ref_table$name %in% rownames(map_tc)) )
resampling_hmaps<- readRDS("sign_matrices_consensus_over_resampling27.rds")
for (i in seq_along(resampling_hmaps)) rownames(resampling_hmaps[[i]])<-
                                      colnames(resampling_hmaps[[i]])<- ref_table$name
resampling_hmaps$t<- resampling_hmaps$t[left_names, left_names]
resampling_hmaps$d<- resampling_hmaps$d[left_names, left_names]
resampling_hmaps$h<- resampling_hmaps$h[left_names, left_names]
interaction_heatmap(map_tc)
interaction_heatmap(resampling_hmaps$t)
interaction_heatmap(map_dc)
interaction_heatmap(resampling_hmaps$d)

graded_sign_matrix<- function(base_map, resmapled_map){
 graded<- base_map + resmapled_map
  print(table(graded))
  return(graded)
}



gmap_t<-graded_sign_matrix(map_tc, resampling_hmaps$t)
gmap_h<-graded_sign_matrix(map_hc, resampling_hmaps$h)
gmap_d<-graded_sign_matrix(map_dc, resampling_hmaps$d)

write.csv(map_tc,"heatmap_total.csv")
write.csv(map_hc,"heatmap_healthy.csv")
write.csv(map_dc,"heatmap_disease.csv")
write.csv(resampling_hmaps$t,"heatmap_total_resample.csv")
write.csv(resampling_hmaps$h,"heatmap_healthy_resample.scv")
write.csv(resampling_hmaps$d,"heatmap_disease_resample.csv")



graded_interaction_heatmap<- function(gmap_matrix, main="interaction heatmap"){
  rownames(gmap_matrix)<-colnames(gmap_matrix)<-NULL
  heatmap(gmap_matrix, Rowv=NA,  Colv=NA, scale="none", revC=TRUE, col=c(rgb(1,0,0), rgb(1,0.6,0.6), 
                                                                        "white",
                                                                        rgb(0.6,0.6,1),rgb(0,0,1)),
          main=main)
}

#PUBLICATION PLOTS
pdf("abundance_map_graded_total27.pdf")
graded_interaction_heatmap(gmap_t, "total data")
dev.off()
pdf("abundance_map_graded_healthy27.pdf", height=8)
par(cex.main=3)
graded_interaction_heatmap(gmap_h, "healthy")
dev.off()
pdf("abundance_map_graded_disease27.pdf", height=8)
par(cex.main=3)
graded_interaction_heatmap(gmap_d, "allergic")
dev.off()

