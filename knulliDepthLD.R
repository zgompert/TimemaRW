## looking at coverage and LD on chromosome 11


## orindation of genotype data, and formatting for gemma
library(data.table)
L<-64650
N<-138

## T. knulli genetics and depth
dep<-as.matrix(fread("../genotypes_rw/depthKnulli.txt",header=FALSE))
g<-as.matrix(fread("../genotypes_rw/Entropy/G_tknulli.txt",sep=",",header=FALSE))
dat<-read.table("/uufs/chpc.utah.edu/common/home/u6000989/projects/timema_confiers/timemaRw/RedwoodDatCombined.csv",header=TRUE,sep=",")
ph<-dat[1:138,]

snps<-read.table("TknulliSnps.txt",header=FALSE)
lg11<-which(snps[,1]==500)
lg4<-which(snps[,1]==6886)

## pc1 is cluster
ko<-kmeans(pco$x[,1],centers=3,iter.max=100,nstart=50)

par(mfrow=c(3,1))
mns<-matrix(NA,nrow=3,ncol=length(lg11))
for(i in 1:3){
        k<-which(ko$cluster==i)
        plot(apply(dep[lg11,k],1,mean),pch=19)
        mns[i,]<-apply(dep[lg11,k],1,mean)
}

pos<-snps[lg11,2]
## 100 windows of ~25 SNPs each on average
wins<-seq(1,max(pos),length.out=100)
Nw<-length(wins)
wmns<-matrix(NA,nrow=3,ncol=Nw-1)
for(i in 1:3){
        for(j in 1:(Nw-1)){
                xx<-which(pos > wins[j] & pos <= wins[j+1])
                wmns[i,j]<-mean(mns[i,xx])
        }
}

plot(wins[-Nw],wmns[1,],pch=19)
points(wins[-Nw],wmns[3,],pch=19,col="cadetblue")
segments(wins[-Nw],wmns[3,],wins[-Nw],wmns[1,])
lb<-c(12696313,15226220)
ub<-c(42020457,44744971)
abline(v=c(lb,ub),col="red")

## differential coverage at one window within the region but not right at bounds
## windows are ~500kb

## LD
LD<-vector("list",3)
for(i in 1:3){
    k<-which(ko$cluster==i & dat$Population[1:138]!="BCTURN") ## drop BCTURN for LD
    G<-g[k,lg11]    
    LD[[i]]<-cor(G,G)
    }

library(RColorBrewer)    
bnds<-c(-1,0.005,0.01,0.05,0.1,0.5,1)
bnds<-c(-1,0.1,0.5,1)
cs<-brewer.pal(n=6,"Reds")[c(1,5,6)]

image(LD[[1]]^2,col=cs,breaks=bnds)
image(LD[[3]]^2,col=cs,breaks=bnds)
image(LD[[2]]^2,col=cs,breaks=bnds)
   
## 100 kb windows
wins100kb<-seq(1,max(pos),100000)
N100<-length(wins100kb)
pos<-snps[lg11,2]
winLd<-matrix(NA,nrow=3,ncol=N100-1)
for(i in 1:3){
    for(j in 1:(N100-1)){
      xx<-which(pos > wins100kb[j] & pos <= wins100kb[j+1])
      subM<-LD[[i]][xx,xx]^2
      winLd[i,j]<-mean(subM[upper.tri(subM)],na.rm=TRUE)
    }
}   

pdf("LDplot.pdf",width=7,height=6)
par(mar=c(5,5,1,1))
mids<-(wins100kb[-1]+wins100kb[-N100])/2
plot(mids,winLd[1,],type='h',ylim=c(-1,1),axes=FALSE,xlab="Position (bp on Chrom. 11)",ylab="LD",cex.lab=1.5)
points(mids,-1*winLd[3,],type='h',col="cadetblue")
abline(v=c(lb,ub),col="red",lty=3)
axis(1)
axis(2,at=c(-1,-.5,0,.5,1),c(1,.5,0,.5,1))
box()
dev.off()

pdf("LDdifplot.pdf",width=7,height=6)
par(mar=c(5,5,1,1))
mids<-(wins100kb[-1]+wins100kb[-N100])/2
plot(mids,winLd[1,]-winLd[3,],type='l',xlab="Position (bp on Chrom. 11)",ylab="LD difference",cex.lab=1.5)
abline(h=0,lty=2)
abline(v=c(lb,ub),col="red",lty=3)
dev.off()


#################################################      
   
## LG 4 for comp.
LD<-vector("list",3)
for(i in 1:3){
    k<-which(ko$cluster==i & dat$Population[1:138]!="BCTURN") ## drop BCTURN for LD
    G<-g[k,lg4]    
    LD[[i]]<-cor(G,G)
    }

library(RColorBrewer)    
bnds<-c(-1,0.005,0.01,0.05,0.1,0.5,1)
bnds<-c(-1,0.1,0.5,1)
cs<-brewer.pal(n=6,"Reds")[c(1,5,6)]

image(LD[[1]]^2,col=cs,breaks=bnds)
image(LD[[3]]^2,col=cs,breaks=bnds)
image(LD[[2]]^2,col=cs,breaks=bnds)
      
