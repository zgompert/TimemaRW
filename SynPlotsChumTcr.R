library(data.table)

## read synteny dat
dat<-fread("out_synteny_CrisChum.psl",header=FALSE)
dfdat<-as.data.frame(dat)



tab<-tapply(X=dfdat[,1],INDEX=list(qg=dfdat[,10],tg=dfdat[,14]),sum)
cris_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,52,4)])
chum_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,40,4)])

ntab<-tab
for(i in 1:10){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}

## order by cris chrom
chord<-c(12,6,10,9,8,13,1,11,7,5,2,3,4)
ntab<-ntab[,chord]

pdf("SynTcrisTchum.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash",ylab="T. cristinae",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,1:13,las=2)
axis(1,at=seq(0,10,length.out=10)/10,chum_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
chtab<-matrix(c(8483,43,
	14640,1392,
	42935,43,
	42912,43,
	18722,56,
	9928,1469,
	10660,1510,
	7748,113,
	16151,43,
	14160,1213,
	12033,48,
	12380,1403,
	14101,1308),nrow=13,ncol=2,byrow=TRUE)

#bnds<-c(13093370,43606674)

chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,926504516,135538230,129098588,79636963,208324082)
crisSize<-c(69933647,70603406,147781425,151833030,69672476,72868255,63432198,97122551,69324240,68842367,65852211,38417004,134256509)

chn<-c(1:13)
pdf("AlnPlotsChumCris.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tcr[tcr %in% tch]
	subd<-dfdat[cc,]
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(as.numeric(subd[1,16:17]),as.numeric(subd[1,12:13]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,ylab="T. chumash",xlab="T. cristinae")
	title(main=paste("Chrom.",chn[i]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,16:17],subd[j,12:13])
		}
		else{
			lines(crisSize[i]-subd[j,16:17],subd[j,12:13],col="cadetblue")
		}
	}
#	if(chn[i]==11){
 #       	abline(v=bnds,col="red",lwd=1.5)
 #       }

}
dev.off()
