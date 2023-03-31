## R script to test SV/fusion genotype effect on performance
library(data.table)
L<-64650
N<-138

## T. knulli genetics and depth
g<-as.matrix(fread("../genotypes_rw/Entropy/G_tknulli.txt",sep=",",header=FALSE))
dat<-read.table("/uufs/chpc.utah.edu/common/home/u6000989/projects/timema_confiers/timemaRw/RedwoodDatCombined.csv",header=TRUE,sep=",")
ph<-dat[1:138,]

snps<-read.table("TknulliSnps.txt",header=FALSE)
P<-apply(g[dat$Population[1:138]!="BCTURN",],2,mean)/2
lg11<-which(snps[,1]==500)
lg11Com<-which(snps[,1]==500 & P > 0.05 & P < 0.95)


pc<-prcomp(g[dat$Population[1:138]!="BCTURN",lg11],center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3,nstart=10,iter.max=100)
gen<-ko$cluster ## these are sorted
plot(pc$x[,1],pc$x[,2])
text(pc$x[,1],pc$x[,2],gen)
save(list=ls(),file="ld.rdat")## can load to get clusters ordered right

gbce<-g[dat$Population[1:138]!="BCTURN",lg11Com]
#gbce<-g[dat$Population[1:138]!="BCTURN",lg11]

ld1<-cor(gbce[gen==1,])^2
ld2<-cor(gbce[gen==2,])^2
ld3<-cor(gbce[gen==3,])^2
ldX<-cor(gbce[gen!=2,])^2
ld1[upper.tri(ld1)]<-NA
ld1[diag(ld1)]<-NA
ld2[upper.tri(ld2)]<-NA
ld2[diag(ld2)]<-NA
ld3[upper.tri(ld3)]<-NA
ld3[diag(ld3)]<-NA
ldX[upper.tri(ldX)]<-NA
ldX[diag(ldX)]<-NA


brks<-c(-.1,seq(0.1,0.5,0.1),1.1)
cs<-c(brewer.pal("YlOrRd",n=6))
#brks<-c(-.1,seq(0.1,0.9,0.1),1.1)
#cs<-c(brewer.pal("YlOrRd",n=9),"black")
#brks<-c(-.1,0.01,0.05,0.1,0.2,0.5,1.1)
#cs<-brewer.pal("YlOrRd",n=6)
pdf("ClusterLDHet.pdf",width=4,height=4)

par(mar=c(1.5,1.5,0.5,0.5))
image(ld2,axes=FALSE,col=cs,breaks=brks)
legend(.1,.95,fill=cs,c(seq(0.0,0.4,0.1),"0.5+"),bty='n',title=expression(r^2))
dev.off()

pdf("ClusterLD.pdf",width=8,height=8)

par(mfrow=c(2,2))
par(mar=c(1.5,1.5,2.5,2.5))
image(ld1,axes=FALSE,col=cs,breaks=brks)
legend(.1,.95,fill=cs,c(seq(0.0,0.4,0.1),"0.5+"),bty='n',title=expression(r^2))
title(main="(A) RW/RW homozygote",cex.main=cm)
#legend(.1,.95,fill=cs,seq(0.05,0.95,0.1))
image(ld2,axes=FALSE,col=cs,breaks=brks)
title(main="(B) RW/C heterozygote",cex.main=cm)
image(ld3,axes=FALSE,col=cs,breaks=brks)
title(main="(C) C/C homozygote",cex.main=cm)
image(ldX,axes=FALSE,col=cs,breaks=brks)
title(main="(D) Rw/RW + C/C homozygotes",cex.main=cm)
dev.off()

gbce<-g[dat$Population[1:138]!="BCTURN",lg11]

h1<-apply(round(gbce[gen==1,],0)==1,2,mean)
h2<-apply(round(gbce[gen==2,],0)==1,2,mean)
h3<-apply(round(gbce[gen==3,],0)==1,2,mean)
q1<-quantile(h1,probs=seq(0,1,0.01))
q2<-quantile(h2,probs=seq(0,1,0.01))
q3<-quantile(h3,probs=seq(0,1,0.01))

perf<-which(snps[lg11,2] >= 13093370 & snps[lg11,2] <= 43606674) 
mean(h1[perf])
#[1] 0.176713
mean(h2[perf])
#[1] 0.2492733
mean(h3[perf])
#[1] 0.04608586

cl<-1.4;cm<-1.4

pdf("ClusterHet.pdf",width=8,height=8)

par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
plot(snps[lg11,2],h1,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.18",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(A) RW/RW homozygote",cex.main=cm)
plot(snps[lg11,2],h2,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.25",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(B) RW/C heterozygote",cex.main=cm)
plot(snps[lg11,2],h3,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.05",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(C) C/C homozygote",cex.main=cm)
plot(snps[lg11,2],h2-(h1+h3)/2,pch=19,col=alpha("black",.3),xlab="Position (bp)",ylab="Difference",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
title(main="(D) Difference in heterozygostiy",cex.main=cm)
abline(h=0)
d<-h2-(h1+h3)/2
wdx<-50:2500
wd<-rep(NA,length(wdx))
for(i in 1:length(wdx)){
	wd[i]<-mean(d[(wdx[i]-49):(wdx[i]+49)])
}
wdxP<-snps[lg11,2][50:2500]
lines(wdxP,wd,lwd=4,col="firebrick")
mean(wd)
#[1] 0.09767095
dev.off()

library(scales)
pdf("ClusterLdHet.pdf",width=12,height=8)

par(mfrow=c(2,3))
par(mar=c(1.5,1.5,2.5,2.5))
image(ld1,axes=FALSE,col=cs,breaks=brks)
legend(.1,.95,fill=cs,c(seq(0.0,0.4,0.1),"0.5+"),bty='n',title=expression(r^2))
title(main="(D) LD RW/RW",cex.main=cm)
#legend(.1,.95,fill=cs,seq(0.05,0.95,0.1))
image(ld3,axes=FALSE,col=cs,breaks=brks)
title(main="(E) LD C/C",cex.main=cm)
image(ldX,axes=FALSE,col=cs,breaks=brks)
title(main="(F) LD RW/RW + C/C",cex.main=cm)
par(mar=c(4.5,5.5,2.5,1.5))
plot(snps[lg11,2],h1,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.18",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(G) Heterozygosity RW/RW",cex.main=cm)
plot(snps[lg11,2],h3,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.05",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(H) Heterozygosity C/C",cex.main=cm)
plot(snps[lg11,2],h2,pch=19,col=alpha("black",.3),ylim=c(0,1),xlab="Position (bp)",ylab="Heterozygosity",cex.lab=cl)
abline(v=c(13093370,43606674),lty=2,col="darkgray",lwd=2)
mtext("0.25",side=3,line=-2,adj=0.05,cex=1.1)
title(main="(I) Heterozygosity RW/C",cex.main=cm)
dev.off()











load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)

## pca RW without BCTURN
pc<-prcomp(t(g_RW_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted


## lm fits, RW only 15 d weight significant... just a test to verify same result and thus clusters as
## lm analysis
summary(lm(gemma_phRW_sub[,1] ~ gen))

#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.018729   0.008554   2.190   0.0343 *
#gen         -0.008547   0.003880  -2.203   0.0333 *
#Residual standard error: 0.01932 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1058,	Adjusted R-squared:  0.08402 
#F-statistic: 4.853 on 1 and 41 DF,  p-value: 0.03328

#brks<-c(-.1,0.01,0.05,0.2,0.5,0.75,1.1)
brks<-c(-.1,0.01,0.05,0.1,0.2,0.5,1.1)
cs<-brewer.pal("YlOrRd",n=6)

ld1<-cor(t(g_RW_sub[a,gen==1]))^2
ld2<-cor(t(g_RW_sub[a,gen==2]))^2
ld3<-cor(t(g_RW_sub[a,gen==3]))^2
ldX<-cor(t(g_RW_sub[a,gen!=2]))^2

ld1[upper.tri(ld1)]<-NA
ld1[diag(ld1)]<-NA
ld2[upper.tri(ld2)]<-NA
ld2[diag(ld2)]<-NA
ld3[upper.tri(ld3)]<-NA
ld3[diag(ld3)]<-NA
ldX[upper.tri(ldX)]<-NA
ldX[diag(ldX)]<-NA

image(ld1,axes=FALSE,col=cs,breaks=brks)
image(ldX,axes=FALSE,col=cs,breaks=brks)

p1<-apply(g_RW_sub[a,gen==1],1,mean)/2
p2<-apply(g_RW_sub[a,gen==2],1,mean)/2
p3<-apply(g_RW_sub[a,gen==3],1,mean)/2
pX<-apply(g_RW_sub[a,gen!=2],1,mean)/2
