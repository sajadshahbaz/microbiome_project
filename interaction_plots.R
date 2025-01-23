load("D:/OwnCloud/Microbiome/alergia/mean_data.RData")

last.name<-function(x,sep=';',n=1) {
 return(unlist(lapply(strsplit(x,sep),FUN=function(x) tail(x,n))))
}

decision<-decision.mean.final

taxon2<-'sk__Bacteria;k__;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;g__Bergeyella'

for (taxon1 in c(
 'sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__Neisseria',
 'sk__Bacteria;k__;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae')
 ) {

 X1<-taxonomy.mean.final[,taxon1]
 short1<-last.name(taxon1)

 X2<-taxonomy.mean.final[,taxon2]
 short2<-last.name(taxon2)

 min1<-min(X1[X1>0])
 max1<-max(X1)
 min2<-min(X2[X2>0])
 max2<-max(X2)
 
 min1r<-10^floor(log10(min1))
 max1r<-10^ceiling(log10(max1))
 min2r<-10^floor(log10(min2))
 max2r<-10^ceiling(log10(max2))
 tx1<-10^(log10(min1r):log10(max1r))
 tx2<-10^(log10(min2r):log10(max2r))
 
 add1<-min(X1[X1>0])/2
 add2<-min(X2[X2>0])/2
 
 add1m<-add1/sqrt(10)
 add2m<-add2/sqrt(10)
 
 prop000<-sum(X1==0 & X2==0 & decision==0)/sum(X1==0 & X2==0)
 prop001<-sum(X1==0 & X2==0 & decision==1)/sum(X1==0 & X2==0)
 prop010<-sum(X1==0 & X2>0 & decision==0)/sum(X1==0 & X2>0)
 prop011<-sum(X1==0 & X2>0 & decision==1)/sum(X1==0 & X2>0)
 prop100<-sum(X1>0 & X2==0 & decision==0)/sum(X1>0 & X2==0)
 prop101<-sum(X1>0 & X2==0 & decision==1)/sum(X1>0 & X2==0)
 prop110<-sum(X1>0 & X2>0 & decision==0)/sum(X1>0 & X2>0)
 prop111<-sum(X1>0 & X2>0 & decision==1)/sum(X1>0 & X2>0)
 
 
 x1<-X1
 x2<-X2
 x1[X1==0]<-x1[X1==0]+exp(log(add1/9)+(log(add1)-log(add1/9))*runif(length(X1[X1==0]))^(2*(!decision[X1==0])+0.5*decision[X1==0]))
 x2[X2==0]<-x2[X2==0]+exp(log(add2/9)+(log(add2)-log(add2/9))*runif(length(X2[X2==0]))^(2*(!decision[X2==0])+0.5*decision[X2==0]))
 x1[X1>0]<-x1[X1>0]+exp(runif(length(X1[X1>0]),log(add1/9),log(add1)))
 x2[X2>0]<-x2[X2>0]+exp(runif(length(X2[X2>0]),log(add2/9),log(add2)))
 
 
 bg.h<-'cyan3'
 bg.d<-'red'
 tx.h<-'black'
 tx.d<-'black'
 par(mar=c(5,5,1,1))
 plot(x1,x2,col=1+decision,type='n',log='xy',cex=1,cex.lab=2.5,
      xlim=c(add1/20,max1r),ylim=c(add2/20,max2r*2),
      xlab=short1,ylab=short2,axes=F)
 axis(1,at=c(0,tx1)+add1m,labels = c(0,tx1),cex.axis=2)
 axis(2,at=c(0,tx2)+add2m,labels = c(0,tx2),cex.axis=2)
 rect(add1/9,add2/9,add1,max2r,col='lightgray',border=NA)
 rect(add1/9,add2/9,max1r,add2,col='lightgray',border=NA)
 
 rect(add1,add2/25,max1r,add2/10,col=bg.d,border=NA)
 rect(add1,add2/25,max1r^prop100*add1^prop101,add2/10,col=bg.h,border=NA)
 text(sqrt(add1*max1r^prop100*add1^prop101),add2/10/sqrt(2.5),sum(X1>0 & X2==0 & decision==0),col=tx.h,cex=2.5)
 text(sqrt(max1r*max1r^prop100*add1^prop101),add2/10/sqrt(2.5),sum(X1>0 & X2==0 & decision==1),col=tx.d,cex=2.5)
 
 rect(add1/10,add2/25,add1,add2/10,col=bg.d,border=NA)
 rect(add1/10,add2/25,add1^prop000*(add1/10)^(prop001),add2/10,col=bg.h,border=NA)
 text(sqrt(add1/10*add1^prop000*(add1/10)^prop001),add2/10/sqrt(2.5),sum(X1==0 & X2==0 & decision==0),col=tx.h,cex=2.5)
 text(sqrt(add1*add1^prop000*(add1/10)^prop001),add2/10/sqrt(2.5),sum(X1==0 & X2==0 & decision==1),col=tx.d,cex=2.5)
 
 rect(add1/25,add2,add1/10,max2r,col=bg.d,border=NA)
 rect(add1/25,add2,add1/10,max2r^prop010*add2^prop011,col=bg.h,border=NA)
 text(add1/10/sqrt(2.5),sqrt(add2*add2^prop011*max2r^prop010),sum(X1==0 & X2>0 & decision==0),col=tx.h,cex=2.5,srt=90)
 text(add1/10/sqrt(2.5),sqrt(max2r*add2^prop011*max2r^prop010),sum(X1==0 & X2>0 & decision==1),col=tx.d,cex=2.5,srt=90)
 
 rect(add1,max2r,max1r,max2r*2.5,col=bg.d,border=NA)
 rect(add1,max2r,max1r^prop110*add1^prop111,max2r*2.5,col=bg.h,border=NA)
 text(sqrt(add1*max1r^prop110*add1^prop111),max2r*sqrt(2.5),sum(X1>0 & X2>0 & decision==0),col=tx.h,cex=2.5)
 text(sqrt(max1r*max1r^prop110*add1^prop111),max2r*sqrt(2.5),sum(X1>0 & X2>0 & decision==1),col=tx.d,cex=2.5)
 
 points(x1,x2,col=c(bg.h,bg.d)[1+decision],pch=19,cex=0.66)
 
 OR_01<-(sum(X1>0 & X2==0 & decision==1)/sum(X1>0 & X2==0 & decision==0))/
  (sum(!(X1>0 & X2==0) & decision==1)/sum(!(X1>0 & X2==0) & decision==0))
 
 print(paste('Odds ratio for',short1,'and',short2,OR_01))
} 
