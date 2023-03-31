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
pdf("AlnPlotsKnulNano.pdf",width=5,height=5)
par(mar=c(4.5,5.5,2.5,1.5))
i<-11
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

invs<-read.table("~/../gompert-group3/data/Tknulli_nanopore/alignment_BCEC-22-4/LG11_inversions.txt",header=FALSE)
invs_st<-invs[,2]
invs_len<-as.numeric(gsub(pattern="SVLEN=",x=invs[,10],replacement=""))
keep<-which(invs_len>1000000)

big_st<-invs_st[keep]
big_end<-invs_st[keep]+invs_len[keep]
for(i in 1:5){
lines(x=c(big_st[i],big_end[i]),y=rep(i*10000000,2),lwd=2.5,col="orange")}

dev.off()

pdf("AlnPlotsKnulNanoMain.pdf",width=5,height=5)
par(mar=c(4.5,5.5,2.5,1.5))
i<-11
tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
tkn<-grep(x=dfdat[,10],pattern=chtab[i,2])
cc<-tcr[tcr %in% tkn]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17])	

plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='l',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.3,xlab="LG11 T. knulli",ylab="LG11 T. cristinae")
title(main="T. knulli variants",cex.main=1.3)
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

invs<-read.table("~/../gompert-group3/data/Tknulli_nanopore/alignment_BCEC-22-4/LG11_inversions.txt",header=FALSE)
invs_st<-invs[,2]
invs_len<-as.numeric(gsub(pattern="SVLEN=",x=invs[,10],replacement=""))
keep<-which(invs_len>1000000)

big_st<-invs_st[keep]
big_end<-invs_st[keep]+invs_len[keep]
i<-1
polygon(x=c(big_st[i],big_end[i],big_end[i],big_st[i]),y=c(2.2e7,2.2e7,5.5e7,5.5e7),border=NA,col=alpha("blue",.25))

dev.off()

