## generate bootstrap replicates of vcf file

## read vcf for perform
vcf<-read.table("/uufs/chpc.utah.edu/common/home/gompert-group3/data/timema_clines_rw_SV/variants_rw_plus/morefilter_filtered2x_tcr_rw_knulli_perform.recode.vcf",comment.char="#",header=FALSE)

Nsnps<-dim(vcf)[1] ## 1728

Njack<-100
bnds<-round(seq(from=1,to=1728,length.out=(Njack+1)),0)

for(i in 1:Njack){
	jack<-vcf[-c(bnds[i]:(bnds[i+1]-1)),] ## 18 SNPs dropped at a time
	write.table(jack,file=paste("jackrep",i,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

