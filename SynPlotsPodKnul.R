library(data.table)

## read synteny dat
dat<-fread("out_synteny_KnulPod.psl",header=FALSE)
dfdat<-as.data.frame(dat)
## 14 big ones for podure
xx<-table(dfdat[,10])
podCh<-names(xx)[xx>1000]
keep<-dfdat[,10] %in% podCh

subDfdat<-dfdat[keep,]


tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
kn_sc<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,36,3)])
tp_sc<-as.numeric(unlist(strsplit(x=rownames(tab),split="[_-]",fixed=FALSE))[seq(2,56,4)])

## normalize with respect to knulli
ntab<-tab
for(i in 1:14){
	ntab[i,]<-ntab[i,]/sum(ntab[i,])
}
## knulli chroms in order of tab
ch<-c(9,1,10,11,6,12,8,4,5,13,2,7)
## chrom vs scaf
cbind(kn_sc,ch)


pdf("SynTknulTpod.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura",ylab="T. knulli",cex.lab=1.4)
axis(2,at=seq(0,12,length.out=12)/12,ch,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tp_sc,las=2)
box()
dev.off()


## colinearity plots for all homologous chromsomes
chtab<-matrix(c(29,12,
	29,2,
	29,3,
	813,229,
	6886,8,
	6886,3,
	6886,12,
	6895,4,
	6839,1,
	934,239,
	6852,29,
	1305,7,
	30,43,
	500,5,
	6840,6,
	775,16),nrow=16,ncol=2,byrow=TRUE)

bnds<-c(13093370,43606674)

knulSize<-c(238587175,63944390,161977359,80749336,83614905,62264078,95486440,79647727,63251091,54078262,43249258,118658088)

chn<-c(1,1,1,2,4,4,4,5:13)
pdf("AlnPlotsPodKnul.pdf",width=10,height=10)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:16){
	tkn<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
	tpo<-grep(x=dfdat[,10],pattern=chtab[i,2])
	cc<-tkn[tkn %in% tpo]
	subd<-dfdat[cc,]
	yub<-max(subd[,13]);xub<-max(subd[,17])	

	plot(as.numeric(subd[1,16:17]),as.numeric(subd[1,12:13]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,ylab="T. podura",xlab="T. knulli")
	title(main=paste("Chrom.",chn[i]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 1:N){
		if(subd[j,9]=="++"){
			lines(subd[j,16:17],subd[j,12:13])
		}
		else{
			lines(knulSize[chn[i]]-subd[j,16:17],subd[j,12:13],col="cadetblue")
			#lines(xub-subd[j,16:17],subd[j,12:13],col="cadetblue")
		}
	}
	if(chn[i]==11){
        	abline(v=bnds,col="red",lwd=1.5)
        }

}
dev.off()
