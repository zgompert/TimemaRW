## compute date based on ML estimates

dest<-read.table("dd_bce_cc_bce_rw_im_old.txt",header=FALSE)

Ttime<-dest[8,6]
Theta<-dest[8,2]
L<-157212 ## number of sites with mean cov > 2 and present 80%, including variable loci, from samtools depth on knulli only data, modest number of sites called variable but filtered out excluded
mu<-(2.9e-9 + 2.8e-9)/2 ## direct estimates from Heliconius and Drosophila, see https://doi.org/10.1093/molbev/msw226

## compute reference (ancestral) population size
## theta = 4 Nref * muL
Nref<-Theta/(4*mu*L)
#[1] 62037
## compute time in generations = years
RealTime<-Nref * Ttime * 2
#1,716,356
## migration rates, M = 2 Nref m
Nref
#[1] 62037
dest[8,7]/(2*Nref)
#[1] 2.96146e-06
dest[8,8]/(2*Nref)
#[1] 8.857139e-07

## for actual data, of SNPs that pass depth and coverage filter (18% of total) about 2/3 (12% of total) fail a different filter whereas 1/3 (6% of total) pass all filters
## here we try to correct for this to get the actual expected number of quality invariant sites
prFail<-0.120138
prPass<-0.06204065
Lc<-L * (prPass)/(prPass+prFail)
Nref<-Theta/(4*mu*Lc)
#182,167.9
## compute time in generations = years
RealTime<-Nref * Ttime * 2
#5,039,976


## jackknife
bdest<-read.table("jackVcf/mlEstimates.txt",header=FALSE)
Ttime<-bdest[,6]
Theta<-bdest[,2]
Nref<-Theta/(4*mu*Lc*.99) ## each jackknife interval is 99% of the total
RealTime<-Nref * Ttime * 2
hist(RealTime)

library(scales)
pdf("dadiPerformTime.pdf",width=4,height=5)
par(mar=c(.5,4.5,.5,.5))
boxplot(RealTime,ylim=c(0,7e6),pch=NA,col="white",xlab="",ylab="Time (MYA)",cex.lab=1.3,axes=FALSE)
axis(2,at=1e6 * 0:7,0:7)
box()
points(jitter(rep(1,100)),RealTime,pch=19,col=alpha("darkgray",.6))
points(1,5039976,pch=19,cex=1.5)
dev.off()

quantile(RealTime,probs=c(.05,.5))
#     5%     50%
#1930561 3414957
