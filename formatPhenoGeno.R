## orindation of genotype data, and formatting for gemma
library(data.table)
L<-64650
N<-138

## T. knulli phenotypes/populations ##
g<-as.matrix(fread("../genotypes_rw/Entropy/G_tknulli.txt",sep=",",header=FALSE))

o<-prcomp(g,center=TRUE,scale=FALSE)
plot(o)
plot(o$x[,1],o$x[,2],pch=19)
## first PC = 22%, 2 clear clusters

dat<-read.table("/uufs/chpc.utah.edu/common/home/u6000989/projects/timema_confiers/timemaRw/RedwoodDatCombined.csv",header=TRUE,sep=",")
ph<-dat[1:138,]

snps<-read.table("TknulliSnps.txt",header=FALSE)
rev(sort(table(snps[,1])))
##   29  6886   775  6852  6839    30  6895  1305   934   500   813  6840  1113 
##15995  9947  5970  5856  4514  4112  3849  3610  3610  2557  2528  1326    16 
##  792     6  6880  1004   121    25  1740   656  2762   697   498   267    48 
##   11    10     8     8     8     8     7     7     6     6     6     6     6 
## There are 12 large scaffolds with many SNPs, these are my LGs, which have not yet been matched up with other LGs

LGs<-as.numeric(names(which(table(snps[,1]) > 100)))

cs<-c("cadetblue","darkolivegreen","chocolate","darkgoldenrod1")
sy<-c(3,4)
#sy<-c(3,NA,4)
pdf("knulliPca.pdf",width=6,height=6)
par(mar=c(5,5,.5,.5))
plot(o$x[,1],o$x[,2],pch=sy[as.numeric(as.factor(ph$Host))],col=cs[as.numeric(as.factor(ph$Population))],xlab="PC1",ylab="PC2",cex.lab=1.4,cex.axis=1.1)
legend(-80,70,c("BCE","BCSH","BCTURN","BCXD"),fill=cs,bty='n')
legend(-80,40,c("C","RW"),pch=c(3,4),bty='n')
dev.off()
## BCTURN clearly different than others, PC2 pulls out BCE on C (versus everyone else)


pdf("TknulLgPCA.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4,5,2,1))

for(i in 1:12){
	a<-which(snps[,1]==LGs[i])
	olg<-prcomp(g[,a],center=TRUE,scale=FALSE)
	pcs<-summary(olg)
	pc1<-round(pcs$importance[2,1]*100,1)
	pc2<-round(pcs$importance[2,2]*100,1)
	plot(olg$x[,1],olg$x[,2],pch=sy[as.numeric(ph$Host)],col=cs[as.numeric(ph$Population)],xlab=paste("PC1 (",pc1,"%)",sep=""),ylab=paste("PC2 (",pc2,"%)",sep=""),cex.lab=1.4,cex.axis=1.1)
	title(main=LGs[i])
}
dev.off()
## now scaffold 500 stands out, which is likely melanic 11 and green striped 12033, but let's see


## pcas along scaffold 500
pdf("TknulScaf500PCAsV2.pdf",width=9,height=12)
par(mfrow=c(4,3))
par(mar=c(4,5,2,1))

a<-which(snps[,1]==500)

out<-prcomp(g[,a],center=TRUE,scale=FALSE)

## get CC BCE homozygotes
bce_cc<-which(out$x[,1] < -6 & out$x[,2] > 5)
## get RWRW BEC homozygotes
bce_rw<-which(out$x[,1] > 10 & out$x[,2] > -5)
## get CC BCTURN homozygotes
bcturn_cc<-which(out$x[,1] < -7 & out$x[,2] < -5)

write.table(ph[bce_cc,1],file="IDS_bce_cc.txt",quote=FALSE,row.names=F,col.names=FALSE)
write.table(ph[bce_rw,1],file="IDS_bce_rw.txt",quote=FALSE,row.names=F,col.names=FALSE)
write.table(ph[bcturn_cc,1],file="IDS_bcturn_cc.txt",quote=FALSE,row.names=F,col.names=FALSE)
# 33 IDS_bce_cc.txt
# 25 IDS_bce_rw.txt
# 21 IDS_bcturn_cc.txt

sy<-c(3,4)
cs<-c("cadetblue","darkolivegreen","chocolate","darkgoldenrod1")
for(i in 1:25){
	aa<-a[((i-1)*100+1):(i*100)]
	slb<-min(snps[aa,2]);sub<-max(snps[aa,2])
	olg<-prcomp(g[,aa],center=TRUE,scale=FALSE)
	pcs<-summary(olg)
	pc1<-round(pcs$importance[2,1]*100,1)
	pc2<-round(pcs$importance[2,2]*100,1)
	plot(olg$x[,1],olg$x[,2],pch=sy[as.numeric(as.factor(ph$Host))],col=cs[as.numeric(as.factor(ph$Population))],xlab=paste("PC1 (",pc1,"%)",sep=""),ylab=paste("PC2 (",pc2,"%)",sep=""),cex.lab=1.4,cex.axis=1.1)
	title(main=paste("Scaffold 500, window",i))
	mtext(paste(slb,sub),3,line=-1.5)
}
dev.off()
## ends of chromosome show geographic pattern, but rest (most) shows clusters suggestive of major SV with freq. differences based on host
## bounds = 12,696,313-15,226,220 ; 42,020,467-44,744,971

## subsets by host treatment and w vs. w/o BCTURN

## first for all
C<-which(ph$Treatment.March.19=="C")
RW<-which(ph$Treatment.March.19=="RW")

phC<-ph[C,]
phRW<-ph[RW,]

## traits are:
## 1. 15 day wgt control sex and stage
## 2. 21 day wgt control sex and stage
## 3. Survival
## 4. 21-15 day wgt change control sex and stage
## 5. 21-15 day wgt change control sex

## use genetic sex instead
sex<-read.table("../genotypes_rw/GeneticSex.txt",header=FALSE)
sexC<-as.factor(sex[C,])
sexRW<-as.factor(sex[RW,])
gemma_phC<-matrix(NA,nrow=length(C),ncol=5)
lmo<-lm(phC$Weight.at.Day.15 ~ phC$Stage.at.Day.15 * sexC)
gemma_phC[is.na(phC$Weight.at.Day.15)==FALSE,1]<-lmo$residuals
lmo<-lm(phC$Weight.at.day.21 ~ phC$Stage.at.Day.21 * sexC)
gemma_phC[is.na(phC$Weight.at.day.21)==FALSE,2]<-lmo$residuals
gemma_phC[,3]<-as.numeric(is.na(phC$Survival==TRUE))
lmo<-lm(phC$Weight.at.day.21-phC$Weight.at.Day.15 ~ phC$Stage.at.Day.15 * sexC)
gemma_phC[is.na(phC$Weight.at.day.21)==FALSE,4]<-lmo$residuals
lmo<-lm(phC$Weight.at.day.21-phC$Weight.at.Day.15 ~ sexC)
gemma_phC[is.na(phC$Weight.at.day.21)==FALSE,5]<-lmo$residuals

gemma_phRW<-matrix(NA,nrow=length(RW),ncol=5)
lmo<-lm(phRW$Weight.at.Day.15 ~ phRW$Stage.at.Day.15 * sexRW)
gemma_phRW[is.na(phRW$Weight.at.Day.15)==FALSE,1]<-lmo$residuals
lmo<-lm(phRW$Weight.at.day.21 ~ phRW$Stage.at.Day.21 * sexRW)
gemma_phRW[is.na(phRW$Weight.at.day.21)==FALSE,2]<-lmo$residuals
gemma_phRW[,3]<-as.numeric(is.na(phRW$Survival==TRUE))
lmo<-lm(phRW$Weight.at.day.21-phRW$Weight.at.Day.15 ~ phRW$Stage.at.Day.15 * sexRW)
gemma_phRW[is.na(phRW$Weight.at.day.21)==FALSE,4]<-lmo$residuals
lmo<-lm(phRW$Weight.at.day.21-phRW$Weight.at.Day.15 ~ sexRW)
gemma_phRW[is.na(phRW$Weight.at.day.21)==FALSE,5]<-lmo$residuals

## now subset
keep<-which(phC$Population != 'BCTURN')
gemma_phC_sub<-gemma_phC[keep,]
keep<-which(phRW$Population != 'BCTURN')
gemma_phRW_sub<-gemma_phRW[keep,]

## write out
write.table(round(gemma_phC,5),file="pheno_C.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(gemma_phC_sub,5),file="pheno_C_sub.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(gemma_phRW,5),file="pheno_RW.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(gemma_phRW_sub,5),file="pheno_RW_sub.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## now subset genotype by ind. similarly
g_C<-t(g[C,])
g_RW<-t(g[RW,])
keep<-which(phC$Population != 'BCTURN')
g_C_sub<-g_C[,keep]
keep<-which(phRW$Population != 'BCTURN')
g_RW_sub<-g_RW[,keep]

write.table(round(g_C,5),file="geno_C.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(g_C_sub,5),file="geno_C_sub.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(g_RW,5),file="geno_RW.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(g_RW_sub,5),file="geno_RW_sub.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


save(list=ls(),file="gp.rdat")

