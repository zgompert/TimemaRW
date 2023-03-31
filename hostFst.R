## Fst between parapatric/sympatric different host Timema

## R version 4.0.2

pfiles<-list.files(pattern="^p_") ## allele freq files
N<-length(pfiles)
P<-vector("list",N)
Lg<-vector("list",N)
Scaf<-vector("list",N)
Pos<-vector("list",N)
for(i in 1:N){
	dat<-read.table(pfiles[i],header=FALSE)
	Nl<-dim(dat)[1]
	Lg[[i]]<-as.numeric(unlist(strsplit(dat[,1],split=":"))[seq(1,by=3,length.out=Nl)])
	Scaf[[i]]<-as.numeric(unlist(strsplit(dat[,1],split=":"))[seq(2,by=3,length.out=Nl)])
	Pos[[i]]<-as.numeric(unlist(strsplit(dat[,1],split=":"))[seq(3,by=3,length.out=Nl)])
	P[[i]]<-as.numeric(dat[,3])
}

## pairs
spps<-matrix(c(1,2, ## cali SM M x Q
	3,6, ## knulli BCE RW x BCWP C
	4,5, ## knulli BCTRU P x C
	7,8, ## land BCBOG C x Q
	9,10, ## land BCSUM C x Q
	11,12), ## popp TBARN DF x RW
	nrow=6,ncol=2,byrow=TRUE)

fst<-function(p1,p2){
	pbar<-(p1+p2)/2
	Ht<-2*pbar*(1-pbar)
	Hs<-p1 * (1-p1) + p2 * (1-p2)
	Fst<-(Ht-Hs)/Ht
	return(Fst)
}

Np<-dim(spps)[1]
genFst<-vector("list",Np)
for(i in 1:Np){
	genFst[[i]]<-fst(P[[spps[i,1]]],P[[spps[i,2]]])
}

## means by Lg
LgFst<-vector("list",Np)
LgFstSd<-vector("list",Np)
LgFst10<-vector("list",Np)
LgFst90<-vector("list",Np)
for(i in 1:Np){
	LgFst[[i]]<-tapply(X=genFst[[i]],INDEX=Lg[[spps[i,1]]],mean,na.rm=TRUE)
	LgFstSd[[i]]<-tapply(X=genFst[[i]],INDEX=Lg[[spps[i,1]]],sd,na.rm=TRUE)
	LgFst10[[i]]<-tapply(X=genFst[[i]],INDEX=Lg[[spps[i,1]]],quantile,probs=.10,na.rm=TRUE)
	LgFst90[[i]]<-tapply(X=genFst[[i]],INDEX=Lg[[spps[i,1]]],quantile,probs=.90,na.rm=TRUE)
}


pdf("lgFstHostDif.pdf",width=8,height=10)
tit<-c("(a) T. californicum, MxQ","(b) T. knulli, CxRW","(c) T. knulli, CxP","(d) T. landelsensis, CxQ","(e) T. landelsensis CxQ","(f) T. poppensis, DFxRW")

par(mfrow=c(3,2))
par(mar=c(5,5,2.5,1))
for(i in 1:Np){
	dotchart(as.numeric(LgFst[[i]][-1]),xlim=c(0,.65),pch=19,labels=1:13,xlab="FST")
	title(main=tit[i])
	segments(LgFst[[i]][-1],1:13,LgFst90[[i]][-1],1:13)
}
dev.off()


