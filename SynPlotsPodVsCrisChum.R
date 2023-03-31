library(data.table)

### podura vs cristinae ####

## read synteny dat
dat<-fread("out_synteny_CrisGSPod.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## 14 big ones for podura
xx<-table(dfdat[,10])
podCh<-names(xx)[xx>1000]
keep<-dfdat[,10] %in% podCh

subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
cr_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,52,4)])
tp_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="[_-]",fixed=FALSE))[seq(2,56,4)])

## normalize with respect to knulli
ntab<-tab
for(i in 1:14){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
## cris chroms in order of tab
ch<-c(7,11,12,13,10,2,9,5,4,3,8,1,6)
## chrom vs scaf
cbind(cr_sc,ch)


pdf("SynTcrisGSTpod.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura",ylab="T. cristinae",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,ch,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tp_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
chtab<-matrix(c(8483,12,
	8483,3,
	14640,229,
	42935,2,
	42935,12,
	42935,3,
	42912,8,
	18722,4,
	9928,1,
	10660,239,
	7748,29,
	16151,7,
	14160,43,
	12033,5,
	12380,6,
	14101,16),nrow=16,ncol=2,byrow=TRUE)

crisSize<-c(69933647,70603406,147781425,151833030,69672476,72868255,63432198,97122551,69324240,68842367,65852211,38417004,134256509)

chn<-c(1,1,2,3,3,3,4,5:13)
pdf("AlnPlotsPodCrisGs.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:16){
	cat(i,"\n")
	tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tpo<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tcr[tcr %in% tpo]
	subd<-dfdat[cc,]
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(as.numeric(subd[1,16:17]),as.numeric(subd[1,12:13]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,ylab="T. podura",xlab="T. cristinae")
	title(main=paste("Chrom.",chn[i]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,16:17],subd[j,12:13])
		}
		else{
			lines(crisSize[chn[i]]-subd[j,16:17],subd[j,12:13],col="cadetblue")
			#lines(xub-subd[j,16:17],subd[j,12:13],col="cadetblue")
		}
	}

}
dev.off()

### podura vs chumash ####

## read synteny dat
dat<-fread("out_synteny_CrisGSPod.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## 14 big ones for podura
xx<-table(dfdat[,10])
podCh<-names(xx)[xx>1000]
keep<-dfdat[,10] %in% podCh

subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
cr_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,52,4)])
tp_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="[_-]",fixed=FALSE))[seq(2,56,4)])

## normalize with respect to knulli
ntab<-tab
for(i in 1:14){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
## cris chroms in order of tab
ch<-c(7,11,12,13,10,2,9,5,4,3,8,1,6)
## chrom vs scaf
cbind(cr_sc,ch)


pdf("SynTcrisGSTpod.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura",ylab="T. cristinae",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,ch,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tp_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
chtab<-matrix(c(8483,12,
	8483,3,
	14640,229,
	42935,2,
	42935,12,
	42935,3,
	42912,8,
	18722,4,
	9928,1,
	10660,239,
	7748,29,
	16151,7,
	14160,43,
	12033,5,
	12380,6,
	14101,16),nrow=16,ncol=2,byrow=TRUE)

crisSize<-c(69933647,70603406,147781425,151833030,69672476,72868255,63432198,97122551,69324240,68842367,65852211,38417004,134256509)

chn<-c(1,1,2,3,3,3,4,5:13)
pdf("AlnPlotsPodCrisGs.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:16){
	cat(i,"\n")
	tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tpo<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tcr[tcr %in% tpo]
	subd<-dfdat[cc,]
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(as.numeric(subd[1,16:17]),as.numeric(subd[1,12:13]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,ylab="T. podura",xlab="T. cris")
	title(main=paste("Chrom.",chn[i]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,16:17],subd[j,12:13])
		}
		else{
			lines(crisSize[chn[i]]-subd[j,16:17],subd[j,12:13],col="cadetblue")
			#lines(xub-subd[j,16:17],subd[j,12:13],col="cadetblue")
		}
	}

}
dev.off()
