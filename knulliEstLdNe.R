## estimating Ne for ABC analyses based on interchromosomal LD
## source ld ne functions
source("ldNe.R")

## read data
library(data.table)
L<-64650
N<-138

## T. knulli genetics and depth
g<-as.matrix(fread("../genotypes_rw/Entropy/G_tknulli.txt",sep=",",header=FALSE))
dat<-read.table("/uufs/chpc.utah.edu/common/home/u6000989/projects/timema_confiers/timemaRw/RedwoodDatCombined.csv",header=TRUE,sep=",")
ph<-dat[1:138,]

snps<-read.table("TknulliSnps.txt",header=FALSE)
lg11<-which(snps[,1]==500)
## drop LG11=scaffold 500
## round to nearest integer for LD
subSnps<-snps[-lg11,]
subG<-round(g[,-lg11],0)

## ids for each
pops<-vector("list",3)
pops[[1]]<-which(ph$Host=="C" & ph$Population=="BCTURN") ## 37 
pops[[2]]<-which(ph$Host=="RW" & ph$Population=="BCE") ## 24
pops[[3]]<-which(ph$Host=="C" & ph$Population=="BCE") ## 68

CC<-20000 ## number of independent pairs

ne<-matrix(NA,nrow=3,ncol=1000)
for(i in 1:1000){
	## compute delta for CC pairs of independent loci
	r2<-matrix(NA,nrow=3,ncol=CC)
	Ns<-dim(subSnps)[1]
	for(k in 1:3){
		j<-0
		while(j<CC){
			sp<-sample(1:Ns,2,replace=FALSE)
			sc<-subSnps[sp,1]
			if(sc[1] != sc[2]){## verify independent
				j<-j+1
				r2[k,j]<-compDelta(subG[pops[[k]],sp[1]],subG[pops[[k]],sp[2]])
				if(is.na(r2[k,j]) | is.finite(r2[k,j])==FALSE){
					j<-j-1
				}
				else{
					r2[k,j]<-r2[k,j]^2
				}
			}
		}
	}

	mnLd<-apply(r2,1,mean)
	ne[1,i]<-NeLargeS(S=37,r2=mnLd[1])
	ne[2,i]<-NeSmallS(S=24,r2=mnLd[2])
	ne[3,i]<-NeLargeS(S=68,r2=mnLd[3])
}

apply(ne,1,summary)
#            [,1]      [,2]     [,3]
#Min.    19.19155  65.63491 63.15598
#1st Qu. 22.15803  96.37296 72.01282
#Median  23.14748 113.32965 75.34265
#Mean    23.22062 124.98650 75.68692
#3rd Qu. 24.15100 137.49077 79.07493
#Max.    29.06530 522.44498 95.09986


save(list=ls(),file="ne.rdat")
