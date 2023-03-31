## annotation summary for perform locus
## perform = scaffold 500, 13,093,370 to 43,606,674 

library(scales)

lb<-13093370 
ub<-43606674
nano_lb<-9706606
nano_ub<-48357002

gtf<-read.table("clean_braker.gtf",header=FALSE)

sum(gtf[,3]=="gene") ## number of "genes"
#[1] 36055

table(gtf[,3])

#        CDS        exon        gene      intron start_codon  stop_codon 
#     154239      154236       36055      116983       37012       37009 
# transcript 
#      36892 


gtf_gene_perform<-gtf[which((gtf[,1]=="ScFnSIn_500_HRSCAF958") & (gtf[,5] >= lb) & (gtf[,4] <= ub) & (gtf[,3]=="gene")),1:5] 
dim(gtf_gene_perform)
#[1] 876   5
# = 876 genes identified within Perform
summary(gtf_gene_perform[,5]-gtf_gene_perform[,4]) ## gene size distribution
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    200     956    3618    6884    8997   75805 


gtf_cds_perform<-gtf[which((gtf[,1]=="ScFnSIn_500_HRSCAF958") & (gtf[,5] >= lb) & (gtf[,4] <= ub) & (gtf[,3]=="CDS")),1:5] 
dim(gtf_cds_perform)
#[1] 3781    5
# i.e., 4.31621 CDS per gene

gtf_gene_lg11<-gtf[which((gtf[,1]=="ScFnSIn_500_HRSCAF958")  & (gtf[,3]=="gene")),1:5] 
dim(gtf_gene_lg11)
#[1] 2153    5

pdf("PerformGeneDensity.pdf",width=5,height=4.5)
par(mar=c(4.5,4.5,.5,.5))
plot((gtf_gene_lg11[,4]+gtf_gene_lg11[,5])/2,log10(gtf_gene_lg11[,5]-gtf_gene_lg11[,4]),pch=19,col=alpha("black",.3),xlab="Position (bp)",ylab="Gene length (log10 bp)",cex.lab=1.3)
abline(v=c(lb,ub),col="firebrick",lwd=2)
dev.off()

bnds<-seq(1,54070627,by=500000)
gd<-rep(NA,108)
for(i in 1:108){
	gd[i]<-sum(gtf_gene_lg11[,5] >= bnds[i] & gtf_gene_lg11[,4] <= bnds[i+1])
}
mids<-(bnds[-109]+bnds[-1])/2
pdf("PerformGeneDensityV2.pdf",width=5,height=4.5)
par(mar=c(4.5,4.5,.5,.5))
plot(mids,gd,xlab="Position (bp)",ylab="Number of genes",cex.lab=1.3,type='l')
abline(v=c(lb,ub),col="firebrick",lwd=2)
dev.off()


## many of the timema chromosomes have fewer genes in the center
pdf("GeneDensityKnulliChroms.pdf",width=8,height=11)
par(mfrow=c(4,3))
par(mar=c(4.5,5,2.5,.5))
lgs<-table(gtf[,1])
chr<-c(9,1,10,11,6,12,8,4,5,13,2,7)
for(k in 1:13){
	j<-which(chr==k)
	if(length(j) > 0){
	gtf_gene_sub<-gtf[which((gtf[,1]==names(lgs)[j])  & (gtf[,3]=="gene")),1:5]
	bnds<-seq(1,max(gtf_gene_sub[,5]),by=500000)
	ll<-length(bnds)-1
	gd<-rep(NA,ll)
	for(i in 1:ll){
		gd[i]<-sum(gtf_gene_sub[,5] >= bnds[i] & gtf_gene_sub[,4] <= bnds[i+1])
	}
	mids<-(bnds[-(ll+1)]+bnds[-1])/2
	plot(mids,gd,xlab="Position (bp)",ylab="Number of genes",cex.lab=1.3,type='l')
	title(main=paste("Chromosome ",k,sep=""),cex.main=1.2)
	if(k==11){
		abline(v=c(lb,ub),col="firebrick",lwd=2)
	}}
}
dev.off()

## get the protein details
ano<-read.table("chrom11_annots_clean.txt",header=FALSE,fill=TRUE)
r1<-which(ano[,1] >= (lb-5000000) & ano[,2] <= (lb+5000000))
r3<-which(ano[,1] >= (ub-5000000) & ano[,2] <= (ub+5000000))
r2<-which(ano[,1] >= (lb+5000000) & ano[,2] <= (ub-5000000))
length(r1)
#[1] 317 
length(r2)
#[1] 382
length(r3)
#[1] 537

write.table(file="ch11BndLeft.txt",data.frame(names(table(ano[r1,3])),table(ano[r1,3])),quote=FALSE,row.names=FALSE)
write.table(file="ch11BndRight.txt",data.frame(names(table(ano[r3,3])),table(ano[r3,3])),quote=FALSE,row.names=FALSE)
write.table(file="ch11BndMid.txt",data.frame(names(table(ano[r2,3])),table(ano[r2,3])),quote=FALSE,row.names=FALSE)

