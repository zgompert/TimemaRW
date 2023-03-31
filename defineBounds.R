## R script to test SV/fusion genotype effect on performance
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)
g_comb_sub<-cbind(g_RW_sub[a,],g_C_sub[a,])
lg11_snps<-snps[a,]

o<-prcomp(t(g_comb_sub),center=TRUE,scale=TRUE)

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

pdf("SVlimPCAeigen.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(pos,eig,type='n',xlab="Position (bp)",ylab="Eigenvalue",cex.lab=1.5,cex.axis=1.1)
a<-1:(bnds[1])
b<-bnds[2]:length(eig)
lines(pos[a],eig[a],col="darkgray",lwd=1.6)
lines(pos[bnds[1]:bnds[2]],eig[bnds[1]:bnds[2]],col="firebrick",lwd=1.6)
lines(pos[b],eig[b],col="darkgray",lwd=1.6)
dev.off()

## boundaries of SV
pos[bnds]
#[1] 13093370 43606674

#$Pi
#             [,1]       [,2]
#[1,] 0.9986083900 0.00139161
#[2,] 0.0005776114 0.99942239

#$pm
#$pm$mean
#[1] 2.827240 4.638348

#$pm$sd
#[1] 0.2903859 0.4466655


