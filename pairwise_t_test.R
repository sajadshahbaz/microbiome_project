
# TO DO: leave only function used in base_resample_analysis.R

### DOES NOT WORK :(

### variance formula (mean of squares - square of means) is too unstable numerically 

pairwise_interaction_by_t_test<- function(abundancy,presence){
 A<-abundancy
 P<-presence
 nP<- 1- P #<- absence
 
 t(P)-> t.P #one row per one taxa
 t(A)-> t.A #same
 
 N0_ij<- t.P %*% nP # number of times i was present and j absent
 N1_ij<- t.P %*% P  # same but j also present
 small_N<- (N0_ij<2) | (N1_ij<2) # these ones we skip
 
 t.Pxt.A<- t.P * t.A # helper, zero out abundance if i is absent
  
 A0_ij<- ( t.Pxt.A %*% nP )/N0_ij #mean abundance of i when i is present and j absent
 A1_ij<- ( t.Pxt.A %*% P )/N1_ij #mean abundance of i when i is present and j present
 A0_ij[small_N]=NA
 A1_ij[small_N]=0
 
 
 t.Pxt.SA <- t.P * t.A^2 #zeroed out squares of abundances if i is absent 
 SA0_ij <- t.Pxt.SA %*% nP #sums of squares of ab. of i when i is present and j absent
 SA1_ij <- t.Pxt.SA %*% P  #sums of squares of ab. of i when i is present and j present
 
 
 
 V0_ij<- (1/(N0_ij-1))*( SA0_ij - N0_ij*A0_ij^2) # unbiased variance of ab. of i when i is present and j absent
 V1_ij<- (1/(N1_ij-1))*( SA1_ij - N1_ij*A1_ij^2) # unbiased variance of ab. of i when i is present and j present
 
 V0_ij[small_N]<-V1_ij[small_N]<-0
 
 # "essentially constant data" like in base R t.test (two-tailed)
 stderr1<- sqrt(V0_ij/N0_ij)
 stderr2<- sqrt(V1_ij/N1_ij)
 stderr12<-sqrt(stderr1^2 + stderr2^2)
 is_constant<- (small_N) | (stderr12 < 10 * .Machine$double.eps * max( abs(A0_ij),abs(A1_ij)))
 
 
 T_ij<- (A0_ij - A1_ij)/sqrt( (V0_ij/N0_ij) + (V1_ij/N1_ij) )  #statistic
 
 T_ij[ is_constant ] = NA
 
 
 DF_ij<- (((V0_ij/N0_ij) + (V1_ij/N1_ij))^2)/
   ( ((V0_ij^2)/((N0_ij^2)*(N0_ij-1)) ) + ((V1_ij^2)/((N1_ij^2)*(N1_ij-1))  ) ) #degrees of freedom for each test
  
 DF_ij [ is_constant ] = NA
 
 PV_ij=2*pt(abs(T_ij),DF_ij, lower.tail = FALSE) #p values
 
 SGN<- (PV_ij<0.05)*ifelse(T_ij>0,1,-1)  #sign of interaction
 SGN[ is_constant ] = 0  #small group sizes / zero variance
 
 return( list(statistic=T_ij,
              p.value=PV_ij,
              sign=SGN)
 )
}


# Function to perform t-test and determine significant increase or decrease
perform_t_test <- function(x, y) {
  # Check if there are enough observations to perform the t-test
  if (length(x) < 2 || length(y) < 2) {
    return(0)
  }
  
  if ( (all(x==0) && all(y==0) ) )
    return(0)
  
  
  t_test_result <- t.test(x, y)
  #  if (is.na(t_test_result$p.value))
  #    {print(table(x)); print(table(y)) }
  # Check if the p-value is less than 0.05
  if (t_test_result$p.value < 0.05) {
    # Check if the mean of y is significantly different from the mean of x
    if (t_test_result$statistic > 0) {
      return(1)
    } else {
      return(-1)
    }
  } else {
    return(0)
  }
}




pairwise_means_sderrs<- function(abundance, presence, disease, report_progress=TRUE, 
                                 min_size=5) {
  mn= min_size -1
  stopifnot(ncol(abundance)==ncol(presence))
  stopifnot(nrow(abundance)==nrow(presence))
  fa= abundance
  fp= presence
  da= abundance[disease==1,]
  dp= presence[disease==1,]
  ha= abundance[disease==0,]
  hp= presence[disease==0,]
  
  mu_ij1_f<-mu_ij0_f<-v_ij1_f<-v_ij0_f<-n_ij1_f<-n_ij0_f<-matrix(ncol=ncol(da),nrow=ncol(da)) #full data
  
  mu_ij1_h<-mu_ij0_h<-mu_ij1_d<-mu_ij0_d<-v_ij1_h<-v_ij0_h<-v_ij1_d<-v_ij0_d<-mu_ij1_f #class specific
  n_ij1_h<-n_ij0_h<-n_ij1_d<-n_ij0_d<-mu_ij1_h
  
  for (i in 1:ncol(da)) {
    if (report_progress) if (i%%10 == 0 ) print(sprintf("%d/%d",i,ncol(da)))
    for (j in 1:ncol(da)) {
      
      f_i1j1<- (fp[,i]==1)& (fp[,j]==1)
      f_i1j0<- (fp[,i]==1)& (fp[,j]==0)
      
      n_ij1_f[i,j]<- sum(f_i1j1)
      n_ij0_f[i,j]<- sum(f_i1j0)
      
      d_i1j1<- (dp[,i]==1)& (dp[,j]==1)
      d_i1j0<- (dp[,i]==1)& (dp[,j]==0)
      h_i1j1<- (hp[,i]==1)& (hp[,j]==1)
      h_i1j0<- (hp[,i]==1)& (hp[,j]==0)
      
      n_ij1_h[i,j]<- sum(h_i1j1)
      n_ij0_h[i,j]<- sum(h_i1j0)
      n_ij1_d[i,j]<- sum(d_i1j1)
      n_ij0_d[i,j]<- sum(d_i1j0)
      
      if (n_ij0_f[i,j]>mn) {
        mu_ij0_f[i,j]= mean(fa[f_i1j0,i])
        v_ij0_f[i,j]= var(fa[f_i1j0,i])
      }
      if (n_ij1_f[i,j]>mn) {
        mu_ij1_f[i,j]= mean(fa[f_i1j1,i])
        v_ij1_f[i,j]= var(fa[f_i1j1,i])
      }
      if (n_ij1_d[i,j]>mn) {
        mu_ij1_d[i,j]= mean(da[d_i1j1,i])
        v_ij1_d[i,j]= var(da[d_i1j1,i])
      }
      if (n_ij0_d[i,j]>mn) {
        mu_ij0_d[i,j]= mean(da[d_i1j0,i])
        v_ij0_d[i,j]= var(da[d_i1j0,i])
      }
      if (n_ij0_h[i,j]>mn) {
        mu_ij0_h[i,j]= mean(ha[h_i1j0,i])
        v_ij0_h[i,j]= var(ha[h_i1j0,i])
      }
      if (n_ij1_h[i,j]>mn) {
        mu_ij1_h[i,j]= mean(ha[h_i1j1,i])
        v_ij1_h[i,j]= var(ha[h_i1j1,i])
      }
    }
  }
  list(n1=n_ij1_f,
       mu1=mu_ij1_f,
       v1=v_ij1_f,
       
       n0=n_ij0_f,
       mu0=mu_ij0_f,
       v0=v_ij0_f,
      
       n1_h=n_ij1_h,
       mu1_h=mu_ij1_h,
       v1_h=v_ij1_h,
       
       n1_d=n_ij1_d,
       mu1_d=mu_ij1_d,
       v1_d=v_ij1_d,
       
       n0_h=n_ij0_h,
       mu0_h=mu_ij0_h,
       v0_h=v_ij0_h,
       
       n0_d=n_ij0_d,
       mu0_d=mu_ij0_d,
       v0_d=v_ij0_d
  )-> RESULT
  RESULT<- lapply(RESULT, t)
  RESULT$taxa_names= colnames(abundance)
  return(RESULT)
}

replace_NA<- function(arr,with_what=0) { arr[is.na(arr)]=with_what; arr } 

difference_in_differences<- function(pairwise_means_sderrs_result,lvl=0.01){
 pairwise_means_sderrs_result-> .r
  t_change=(.r$mu1 - .r$mu0) #total
  d_change=(.r$mu1_d - .r$mu0_d)
  h_change=(.r$mu1_h - .r$mu0_h)
  
  t_change_var=(.r$v1)/(.r$n1) +
    (.r$v0)/(.r$n0)
  d_change_var=(.r$v1_d)/(.r$n1_d) +
    (.r$v0_d)/(.r$n0_d) 
  h_change_var=(.r$v1_h)/(.r$n1_h) +
    (.r$v0_h)/(.r$n0_h) 
  diff_change=  d_change- h_change 
  
  diff_change_var= (.r$v1_d)/(.r$n1_d) +
                (.r$v0_d)/(.r$n0_d) +
                (.r$v1_h)/(.r$n1_h) +
                (.r$v0_h)/(.r$n0_h) 
  
  2*pnorm(abs(t_change),mean = 0, sd= sqrt(t_change_var), lower.tail=FALSE)-> p.value_t
  2*pnorm(abs(d_change),mean = 0, sd= sqrt(d_change_var), lower.tail=FALSE)-> p.value_dis
  2*pnorm(abs(h_change),mean = 0, sd= sqrt(h_change_var), lower.tail=FALSE)-> p.value_h
  2*pnorm(abs(diff_change),mean = 0, sd= sqrt(diff_change_var), lower.tail=FALSE)-> p.value_diff
  
  
  t_sign<- (p.value_t<lvl)*ifelse(t_change>0,1,-1)
  diff_sign<- (p.value_diff<lvl)*ifelse(diff_change>0,1,-1)
  dis_sign<- (p.value_dis<lvl)*ifelse(d_change>0,1,-1)
  h_sign<- (p.value_h<lvl)*ifelse(h_change>0,1,-1)
  
  t_sign<-replace_NA(t_sign)
  diff_sign<-replace_NA(diff_sign)
  h_sign<-replace_NA(h_sign)
  dis_sign<-replace_NA(dis_sign)
  
  
  list(ji_effect_DvsHt=t_change,
       var_ji_effect_DvsHt=t_change_var,
       signt=t_sign,
       ji_effect_DvsH=diff_change,
       var_ji_effect_DvsH=diff_change_var,
       sign=diff_sign,
       ji_effect_DvsHd=d_change, 
       var_ji_effect_DvsHd=d_change_var,
       signd=dis_sign,
       ji_effect_DvsHh=h_change, 
       var_ji_effect_DvsHh=h_change_var,
       signh=h_sign,
       taxa_names=.r$taxa_names
       )
}

get_sign_matrices<- function(pairwise_interaction_matrix_list){
 pairwise_interaction_matrix_list-> .r
  list(total= .r$signt,
       d=.r$signd,
       h=.r$signh,
       taxa_names=.r$taxa_names)
}

summarize_signs<- function(sign_matrix, taxa_names){
  rbind(
  data.frame(name="overall",
             mean=NA,
             freq=NA,
             n_plus=sum(sign_matrix>0),
             n_minus=sum(sign_matrix<0)
  ),
  data.frame(name=taxa_names,
             n_plus=rowSums(sign_matrix>0),
             n_minus=rowSums(sign_matrix<0)
              )
  )
  
}
summarize_abp<- function(abundance,presence){
  
  if (ncol(abundance)!=ncol(presence))
    stop("presence and abundance must have same number of columns")
  data.frame(name=colnames(abundance),
             mean=unname(colMeans(abundance)),
             freq=unname(colMeans(presence))
              )
}

summarize_abp_IG_HvsD<- function(
                          abundance,
                          presence,
                          disease,
                          IG_per_dim
                          ){
  summarize_abp(
                  abundance,
                  presence)-> total
  summarize_abp(
                  abundance[disease==1,],
                  presence[disease==1,])-> d
  summarize_abp(
                  abundance[disease==0,],
                  presence[disease==0,])-> h
  colnames(d)<- paste0(colnames(d),"_d") 
  colnames(h)<- paste0(colnames(h),"_h") 
  cbind(total, d[,2:ncol(d)], h[,2:ncol(h)])-> df
  rel_taxa_per_dim<- lapply(IG_per_dim,names)
  #print(rel_taxa_per_dim)
  df$rel_1D= df$name %in% rel_taxa_per_dim[[1]]
  df$rel_2D= df$name %in% rel_taxa_per_dim[[2]]
  df$IG_1D= rep(NA, nrow(df))
  df$IG_2D= rep(NA, nrow(df))
  for (i in 1:nrow(df))
  { name_i<- df$name[[i]]
    df$IG_1D[[i]]= ifelse(name_i %in% rel_taxa_per_dim[[1]], IG_per_dim[[1]][[name_i]],NA) 
    df$IG_2D[[i]]= ifelse(name_i %in% rel_taxa_per_dim[[2]], IG_per_dim[[2]][[name_i]],NA) 
  }
  df
}



sign_heatmap_display<- function(sign_matrix, pdf_file=NULL, color_code=c("red","white","blue")){
  stopifnot(all(sign_matrix %in% c(-1,0,1)))
  if (!is.null(pdf_file)) pdf(pdf_file)
  heatmap(sign_matrix, Rowv=NA,Colv=NA,scale="none",revC=TRUE, col=color_code)
  if (!is.null(pdf_file)) dev.off()
}




library(magrittr)


#t2sign<-heatmap_list$signd
#t2sign<- t1sign
#t2sign[,]=NA
#t2sign_h<-t2sign
#a_rv<- abundance_sampled[,rv]
#b_rv<- binary_taxa_sampled[,rv]
## Iterate through each taxa in di==1
#for (i in 1:ncol(a_rv)) {
#  if (i%%10 == 0 ) print(paste("i_di=",i))
#  for (j in 1:ncol(a_rv)) {
#    # Filter the dataset to include only rows where the disease is present (disease == 1)
#    filtered_present_di <- a_rv[ (y_sampled==1) &(b_rv[,i]==1) & (b_rv[,j]==1) ,]
#    filtered_absent_di <-a_rv[ (y_sampled==1) &(b_rv[,i]==1) &(b_rv[,j]==0) ,]
#    
#    filtered_present_h<- a_rv[ (y_sampled==0) &(b_rv[,i]==1) & (b_rv[,j]==1) ,]
#    filtered_absent_h<-a_rv[ (y_sampled==0) &(b_rv[,i]==1) &(b_rv[,j]==0) ,]
# 
#    t2sign[i, j] = perform_t_test(filtered_present_di[, i], filtered_absent_di[, i])
#    t2sign_h[i,j]= perform_t_test(filtered_present_h[, i], filtered_absent_h[, i])
#  }
#}
#
#
#dim(t2sign)
#dim(t1sign)
#all(t1sign==t2sign)
#sum(t1sign==t2sign)
#length(t1sign)

###TO DO:
##1) for each relevant V (1D or 2D): 
###  IG (1D & 2D), 
###  freq & abundance per disease status, 
###  interaction partner (2D)
###2) pairwise heatmap:
### - general
### - health
### - disease
### - permuted?
###3) permute taxa? abundance? presence? both? one taxa at the time, keep rest fixed?
###4) cluster patients
### - based on taxa - can we correlate clusters with external clinical variables?
### - based on clinical - can we correlate clusters with bacteria?


#library(umap)
#set.seed(123)
#umap(d = abundance_sampled, n_neigbors=30)$layout->UMAP
#par(mfrow=c(1,1))
#plot(UMAP, col=ifelse(y_sampled,"red","blue"), pch=16)
#metumap.defaults
