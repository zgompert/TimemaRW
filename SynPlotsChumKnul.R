library(data.table)

## read synteny dat
dat<-fread("out_synteny_KnulChum.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## 10 big ones for chumash
xx<-table(dfdat[,10])
chCh<-names(xx)[xx>1000]
keep<-dfdat[,10] %in% chCh

subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
kn_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,36,3)])
tc_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,40,4)])

## normalize with respect to knulli
ntab<-tab
for(i in 1:10){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
## knulli chroms in order of tab
ch<-c(9,1,10,11,6,12,8,4,5,13,2,7)
## chrom vs scaf
cbind(kn_sc,ch)


pdf("SynTknulTchum.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. chumash",ylab="T. knulli",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,ch,las=2)
axis(1,at=seq(0,10,length.out=10)/10,tc_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
chtab<-matrix(c(29,43,
	813,1392,
	6886,43,
	6895,56,
	6839,1469,
	934,1510,
	6852,113,
	1305,43,
	30,1213,
	500,48,
	6840,1403,
	775,1308),nrow=12,ncol=2,byrow=TRUE)

bnds<-c(13093370,43606674)

chumSize<-c(926504516,142035733,926504516,149761397,151976569,125220160,174643345,926504516,135538230,129098588,79636963,208324082)
knulSize<-c(238587175,63944390,161977359,80749336,83614905,62264078,95486440,79647727,63251091,54078262,43249258,118658088)

chn<-c(1,2,4:13)
pdf("AlnPlotsChumKnul.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:12){
	tkn<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tch<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tkn[tkn %in% tch]
	subd<-dfdat[cc,]
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(as.numeric(subd[1,16:17]),as.numeric(subd[1,12:13]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,ylab="T. chumash",xlab="T. knulli")
	title(main=paste("Chrom.",chn[i]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,16:17],subd[j,12:13])
		}
		else{
			lines(knulSize[i]-subd[j,16:17],subd[j,12:13],col="cadetblue")
			#lines(xub-subd[j,16:17],subd[j,12:13],col="cadetblue")
		}
	}
	if(chn[i]==11){
        	abline(v=bnds,col="red",lwd=1.5)
        }

}
dev.off()
