library(data.table)

## read synteny dat
dat<-fread("out_synteny_knulli.psl",header=FALSE)
dfdat<-as.data.frame(dat)

tab<-tapply(X=dfdat[,1],INDEX=list(qg=dfdat[,10],tg=dfdat[,14]),sum)
kn_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,36,3)])
tc_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,50,4)])

## normalize with respect to knulli
ntab<-tab
for(i in 1:12){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
ch<-c(7,11,12,13,10,2,9,5,4,3,8,1,6)
## chrom vs scaf
cbind(tc_sc,ch)
#      tc_sc ch
# [1,] 10660  7
# [2,] 12033 11
# [3,] 12380 12
# [4,] 14101 13
# [5,] 14160 10
# [6,] 14640  2
# [7,] 16151  9
# [8,] 18722  5
# [9,] 42912  4
#[10,] 42935  3
#[11,]  7748  8
#[12,]  8483  1
#[13,]  9928  6
pdf("SynTcrTknul.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. knulli",ylab="T. cristinae",cex.lab=1.4)
axis(2,at=seq(0,13,length.out=13)/13,ch,las=2)
axis(1,at=seq(0,12,length.out=12)/12,kn_sc,las=2)
box()
dev.off()

## get subset Tcr 11 (12033) vs Tknul 500
tcr11<-grep(x=dfdat[,14],pattern="12033") 
tkn500<-grep(x=dfdat[,10],pattern="500")
cc11x500<-tcr11[tcr11 %in% tkn500]

## aprox boudns of SV signal
lb<-c(12696313,15226220)
ub<-c(42020457,44744971)

subd<-dfdat[cc11x500,]

## 12/13 qstart qend
## 16/17 tstart tend
pdf("AlignTcr11Tknul500.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
xub<-54072659
yub<-65781494
plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. knulli",ylab="T. cristinae")

N<-dim(subd)[1]
for(i in 2:N){
	if(subd[i,9]=="++"){
		lines(subd[i,12:13],subd[i,16:17])
	}
	else{
		lines(subd[i,12:13],subd[i,16:17],col="cadetblue")
	}
}
abline(v=c(lb,ub),col="red",lwd=1.5)
dev.off()	

crisSize<-c(69933647,70603406,147781425,151833030,69672476,72868255,63432198,97122551,69324240,68842367,65852211,38417004,134256509)

## colinearity plots for all homologous chromsomes
chtab<-matrix(c(8483,29,
	14640,813,
	42935,29,
	42912,6886,
	18722,6895,
	9928,6839,
	10660,934,
	7748,6852,
	16151,1305,
	14160,30,
	12033,500,
	12380,6840,
	14101,775),nrow=13,ncol=2,byrow=TRUE)

bnds<-c(13093370,43606674)
pdf("AlnPlotsKnulTcr.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:13){
	tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tkn<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tcr[tcr %in% tkn]
	subd<-dfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. knulli",ylab="T. cristinae")
	title(main=paste("Chrom.",i),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],crisSize[i]-subd[j,16:17],col="cadetblue")
		}
	}
	if(i==11){
		abline(v=bnds,col="red",lwd=1.5)
	}
}
dev.off()
