## make popgenplot
library(HiddenMarkov)
library(scales)
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)
g_comb_sub<-cbind(g_RW_sub[a,],g_C_sub[a,])
g_comb<-cbind(g_RW_sub,g_C_sub)
lg11_snps<-snps[a,]
host<-c(rep(1,49),rep(2,52))
## not right yet

bnds<-1:2451
eig<-rep(NA,2451)
pos<-rep(NA,2451)
for(i in 1:2451){
	x<-bnds[i]:(bnds[i]+99)
	o<-prcomp(t(g_comb_sub[x,]),center=TRUE,scale=TRUE)
	eig[i]<-o$sdev[1]
	pos[i]<-mean(lg11_snps[x,2])
}



delta<-c(0.5,0.5)
pii<-matrix(c(0.95,0.05,0.05,0.95),nrow=2,byrow=T)
sds<-rep(sd(eig)*.5,2)
mns<-as.numeric(quantile(eig,probs=c(.25,.75)))

init<-dthmm(eig,pii,delta,"norm",list(mean=mns,sd=sds),discrete=FALSE)
params<-BaumWelch(init,bwcontrol(maxiter = 500, tol = 1e-04, 
                prt = TRUE, posdiff = FALSE,converge = expression(diff < tol)))
fit<-Viterbi(params)

bnds<-summary(which(fit==2))[c(1,6)]
sv<-pos[bnds]
#[1] 13093370 43606674

pops<-as.character(dat$Population)[1:138]
Prw<-apply(g[pops=="BCE" & dat$Host[1:138]=="RW",]/2,2,mean)
Pc<-apply(g[pops=="BCE" & dat$Host[1:138]=="C",]/2,2,mean)
fst<-function(p1,p2){
	pbar<-(p1+p2)/2
	Ht<-2*pbar*(1-pbar)
	Hs<-p1 * (1-p1) + p2 * (1-p2)
	ff<-mean(Ht-Hs)/mean(Ht)
	return(ff)
}

L<-length(Pc)
fvec<-rep(NA,L)
for(i in 1:L){
	fvec[i]<-fst(Prw[i],Pc[i])
}
cmap<-matrix(c(1,29,
	2,813,
	4,6886,
	5,6895,
	6,6839,
	7,934,
	8,6852,
	9,1305,
	10,30,
	11,500,
	12,6840,
	13,775),nrow=12,ncol=2,byrow=TRUE)
smap<-rep(NA,L)
for(i in 1:12){
	smap[snps[,1]==cmap[i,2]]<-cmap[i,1]
}

ofvec<-fvec[order(smap)]
chrom<-smap[order(smap)]
ofvec.nona<-ofvec[-which(is.na(chrom)==TRUE)]
chrom.nona<-chrom[-which(is.na(chrom)==TRUE)]
cbnds<-c(which(chrom.nona[-1] != chrom.nona[-length(chrom.nona)]),length(chrom.nona))
mids<-tapply(1:length(chrom.nona),INDEX=chrom.nona,median)

pdf("fig_tknulli_popgen.pdf",width=9,height=6)
layout(matrix(c(1,1,1,2:4),byrow=TRUE,ncol=3,nrow=2),widths=c(3,3,3),height=c(3,3))
cs<-alpha(c("blue",NA,"indianred4"),c(.75,NA,.75))
cl<-1.5;ca<-1.1;cm<-1.5

par(mar=c(4.5,5.5,2.5,1.5))

plot(ofvec.nona,type='n',xlab="Chromosome",ylab=expression(paste(F[ST])),axes=FALSE,cex.lab=cl)
for(i in seq(1,11,2)){
	polygon(cbnds[c(i,i+1,i+1,i)],c(-.1,-.1,1.1,1.1),col=alpha("gray",.4),border=NA)
}
points(ofvec.nona,pch=19,col=alpha("firebrick2",.4))
ns<-length(ofvec.nona)
wofvec<-rep(NA,ns-99)
for(j in 1:length(wofvec)){
	wofvec[j]<-mean(ofvec.nona[j:(j+99)],na.rm=TRUE)
}
lines(c(1:(ns-99))+49,wofvec,col="firebrick4",lwd=1.5)

axis(1,at=mids,c(1,2,4:12,"X"),cex.axis=ca)
axis(2,cex.axis=ca)
title(main="(a) Genetic differentiation between host races",cex.main=cm)
box()

## (b)
plot(pos,eig,type='n',xlab="Position (bp)",ylab="Eigenvalue",cex.lab=cl,cex.axis=ca,ylim=c(2.2,6.5))
a<-1:(bnds[1])
b<-bnds[2]:length(eig)
lines(pos[a],eig[a],col="darkgray",lwd=1.6)
lines(pos[bnds[1]:bnds[2]],eig[bnds[1]:bnds[2]],col="red",lwd=1.6)
lines(pos[b],eig[b],col="darkgray",lwd=1.6)
mtext("Perform",side=3,line=-2,cex=1)
title(main="(b) Chromosome-11 eigen",cex.main=cm)

## (c)
cs<-c("#0000FFBF","#8B3A3ABF")

pops<-as.character(dat$Population[1:138])
pn<-as.numeric(as.factor(pops))
sym<-c(19,19,3,19)
a<-which(snps[,1]==500)
o<-prcomp(g[,-a],center=TRUE,scale=TRUE)
plot(o$x[,1],o$x[,2],pch=sym[pn],bg="darkgray",col=cs[as.numeric(as.factor(dat$Host[1:138]))],xlab="PC1 (12.6%)",ylab="PC2 (3.4%)",cex.lab=cl,cex.axis=ca)
title(main="(c) Genome-wide PCA",cex.main=cm)
legend(-150,170,c("BCTURN, C","BCE, C","BCE, RW"),pch=c(3,19,19),col=cs[c(1,1,2)],bty='n')

## (d)
x<-which(snps[,1]==500)[bnds[1]:bnds[2]]
o<-prcomp(g[,x],center=TRUE,scale=TRUE)
plot(-1 * o$x[,1],-1 * o$x[,2],pch=sym[pn],bg="darkgray",cex.lab=cl,col=cs[as.numeric(as.factor(dat$Host[1:138]))],xlab="PC1 (20.5%)",ylab="PC2 (8.6%)",cex.main=cm)
title(main="(d) Perform PCA",cex.main=cm)

dev.off()

## boundaries of SV

#$Pi
#             [,1]       [,2]
#[1,] 0.9986083900 0.00139161
#[2,] 0.0005776114 0.99942239

#$pm
#$pm$mean
#[1] 2.827240 4.638348

#$pm$sd
#[1] 0.2903859 0.4466655

## rw homozygotes based on PCA (note bounds are for a speficic polarization from one PCA)
rr<-which(o$x[,1] < -20 & o$x[,2] < 5) ## 25 inds.
cc<-which(o$x[,1] > 10 & o$x[,2] < 0) ## 33 inds ## THESE Match what I have for beast, those were homozygotes!
