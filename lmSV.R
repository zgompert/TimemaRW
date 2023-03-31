## R script to test SV/fusion genotype effect on performance
load("gp.rdat")

## scaffold 500
a<-which(snps[,1]==500)

## pca RW without BCTURN
pc<-prcomp(t(g_RW_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted


## lm fits, RW only 15 d weight significant
summary(lm(gemma_phRW_sub[,1] ~ gen))

#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.018729   0.008554   2.190   0.0343 *
#gen         -0.008547   0.003880  -2.203   0.0333 *
#Residual standard error: 0.01932 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1058,	Adjusted R-squared:  0.08402 
#F-statistic: 4.853 on 1 and 41 DF,  p-value: 0.03328

summary(lm(gemma_phRW_sub[,2] ~ gen))

#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.008068   0.011748   0.687    0.496
#gen         -0.004356   0.005329  -0.818    0.418
#Residual standard error: 0.02653 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.01604,	Adjusted R-squared:  -0.007959 
#F-statistic: 0.6683 on 1 and 41 DF,  p-value: 0.4184

summary(lm(gemma_phRW_sub[,3] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.93683    0.13360   7.012 7.83e-09 ***
#gen         -0.01849    0.05999  -0.308    0.759    
#Residual standard error: 0.3088 on 47 degrees of freedom
#Multiple R-squared:  0.002017,	Adjusted R-squared:  -0.01922 
#F-statistic: 0.09499 on 1 and 47 DF,  p-value: 0.7593
summary(glm(gemma_phRW_sub[,3] ~ gen,family="binomial"))
#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)   2.6137     1.5036   1.738   0.0821 .
#gen          -0.2046     0.6520  -0.314   0.7537  

#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 32.196  on 47  degrees of freedom
#AIC: 36.196
#Number of Fisher Scoring iterations: 5

sexRW_sub<-sexRW[which(phRW$Population != 'BCTURN')]
tapply(X=gemma_phRW_sub[,3],INDEX=list(gen,sexRW_sub),mean)
#     female  male
#1 1.0000000 1.000
#2 0.7857143 0.875
#3 0.9166667 1.000

## pca C without BCTURN
pc<-prcomp(t(g_C_sub[a,]),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted ... ran till they were

## lm fits, C 15 and 21 d weight and survival all significant
summary(lm(gemma_phC_sub[,1] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.041435   0.014971   2.768  0.00800 **
#gen         -0.018215   0.006749  -2.699  0.00958 **
#Residual standard error: 0.03677 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1317,	Adjusted R-squared:  0.1137 
#F-statistic: 7.284 on 1 and 48 DF,  p-value: 0.009578

summary(lm(gemma_phC_sub[,2] ~ gen))
#             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  0.045064   0.018150   2.483   0.0168 *
#gen         -0.019387   0.007982  -2.429   0.0192 *
#Residual standard error: 0.04067 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.1159,	Adjusted R-squared:  0.09625 
#F-statistic: 5.899 on 1 and 45 DF,  p-value: 0.01921

summary(lm(gemma_phC_sub[,3] ~ gen))
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.59451    0.12126   4.903 1.04e-05 ***
#gen          0.14099    0.05519   2.555   0.0137 *  
#---
#Residual standard error: 0.3064 on 50 degrees of freedom
#Multiple R-squared:  0.1154,	Adjusted R-squared:  0.09775 
#F-statistic: 6.526 on 1 and 50 DF,  p-value: 0.01373
summary(glm(gemma_phC_sub[,3] ~ gen,family="binomial"))

#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)  -0.8772     1.2138  -0.723   0.4699  
#gen           1.7119     0.7955   2.152   0.0314 *
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 30.603  on 50  degrees of freedom
#AIC: 34.603
#Number of Fisher Scoring iterations: 6

sexC_sub<-sexC[which(phC$Population != 'BCTURN')]
tapply(X=gemma_phC_sub[,3],INDEX=list(gen,sexC_sub),mean)
#     female      male
#1 0.8571429 0.5714286
#2 0.9375000 0.8000000
#3 1.0000000 1.0000000


## combined
g_comb_sub<-cbind(g_RW_sub[a,],g_C_sub[a,])
pc<-prcomp(t(g_comb_sub),center=TRUE,scale=FALSE)
ko<-kmeans(pc$x[,1],centers=3)
gen<-ko$cluster ## these are sorted
save(list=ls(),file="svLm.rdat")

N_rw<-dim(gemma_phRW_sub)[1]
N_c<-dim(gemma_phC_sub)[1]
host_trt<-rep(c(1,2),c(N_rw,N_c))
summary(lm(gemma_phRW_sub[,1] ~ gen[host_trt==1]))
#                    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.018729   0.008554   2.190   0.0343 *
#gen[host_trt == 1] -0.008547   0.003880  -2.203   0.0333 *
#Residual standard error: 0.01932 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.1058,	Adjusted R-squared:  0.08402 
#F-statistic: 4.853 on 1 and 41 DF,  p-value: 0.03328

summary(lm(gemma_phRW_sub[,2] ~ gen[host_trt==1]))
#                    Estimate Std. Error t value Pr(>|t|)
#(Intercept)         0.008068   0.011748   0.687    0.496
#gen[host_trt == 1] -0.004356   0.005329  -0.818    0.418
#Residual standard error: 0.02653 on 41 degrees of freedom
#  (6 observations deleted due to missingness)
#Multiple R-squared:  0.01604,	Adjusted R-squared:  -0.007959 
#F-statistic: 0.6683 on 1 and 41 DF,  p-value: 0.4184

summary(glm(gemma_phRW_sub[,3] ~ gen[host_trt==1],family="binomial"))
#                   Estimate Std. Error z value Pr(>|z|)  
#(Intercept)          2.6137     1.5036   1.738   0.0821 .
#gen[host_trt == 1]  -0.2046     0.6520  -0.314   0.7537  
#    Null deviance: 32.295  on 48  degrees of freedom
#Residual deviance: 32.196  on 47  degrees of freedom
#AIC: 36.196
#Number of Fisher Scoring iterations: 5

summary(lm(gemma_phC_sub[,1] ~ gen[host_trt==2]))
#                    Estimate Std. Error t value Pr(>|t|)   
#(Intercept)         0.041435   0.014971   2.768  0.00800 **
#gen[host_trt == 2] -0.018215   0.006749  -2.699  0.00958 **
#Residual standard error: 0.03677 on 48 degrees of freedom
#  (2 observations deleted due to missingness)
#Multiple R-squared:  0.1317,	Adjusted R-squared:  0.1137 
#F-statistic: 7.284 on 1 and 48 DF,  p-value: 0.009578

summary(lm(gemma_phC_sub[,2] ~ gen[host_trt==2]))
#                    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.045064   0.018150   2.483   0.0168 *
#gen[host_trt == 2] -0.019387   0.007982  -2.429   0.0192 *
#Residual standard error: 0.04067 on 45 degrees of freedom
#  (5 observations deleted due to missingness)
#Multiple R-squared:  0.1159,	Adjusted R-squared:  0.09625 
#F-statistic: 5.899 on 1 and 45 DF,  p-value: 0.01921

summary(glm(gemma_phC_sub[,3] ~ gen[host_trt==2],family="binomial"))
#                   Estimate Std. Error z value Pr(>|z|)  
#(Intercept)         -0.8772     1.2138  -0.723   0.4699  
#gen[host_trt == 2]   1.7119     0.7955   2.152   0.0314 *
#    Null deviance: 37.193  on 51  degrees of freedom
#Residual deviance: 30.603  on 50  degrees of freedom
#AIC: 34.603
#Number of Fisher Scoring iterations: 6

cm<-1.5;cl<-1.5

pdf("KnulliSexWeightSV.pdf",width=9,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5.5,2.5,1.5))
cs<-rep("gray10",length(sexC_sub))
a<-which(sexC_sub=="male")
cs[a]<-"red"
plot(gen[host_trt==2],gemma_phC_sub[,1],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexC_sub=="male")
o<-lm(gemma_phC_sub[a,1] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexC_sub=="female")
o<-lm(gemma_phC_sub[a,1] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="C, 15d weight",cex.main=cm)

plot(gen[host_trt==2],gemma_phC_sub[,2],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexC_sub=="male")
o<-lm(gemma_phC_sub[a,2] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexC_sub=="female")
o<-lm(gemma_phC_sub[a,2] ~ gen[host_trt==2][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="C, 21d weight",cex.main=cm)


cs<-rep("gray10",length(sexRW_sub))
a<-which(sexRW_sub=="male")
cs[a]<-"red"
plot(gen[host_trt==1],gemma_phRW_sub[,1],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexRW_sub=="male")
o<-lm(gemma_phRW_sub[a,1] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexRW_sub=="female")
o<-lm(gemma_phRW_sub[a,1] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="RW, 15d weight",cex.main=cm)

plot(gen[host_trt==1],gemma_phRW_sub[,2],pch=19,col=alpha(cs,.5),axes=FALSE,xlab="Genotype",ylab="Weight",cex.lab=cl)
a<-which(sexRW_sub=="male")
o<-lm(gemma_phRW_sub[a,2] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="red")
a<-which(sexRW_sub=="female")
o<-lm(gemma_phRW_sub[a,2] ~ gen[host_trt==1][a])
abline(o$coefficients,lwd=1.8,col="gray10")
axis(1,at=1:3,c("0","1","2"))
axis(2)
box()
title(main="RW, 21d weight",cex.main=cm)
dev.off()


##################################################
csC<-c("cadetblue2","cadetblue3","cadetblue4")
csRW<-c("firebrick2","firebrick3","firebrick4")
cl<-1.6;lx<-2.1;cm<-1.6;ca<-1.2;cx<-1.3

pdf("knulliPerformanceSV.pdf",width=12,height=12)
layout(matrix(c(1,1,1,2:7),nrow=3,ncol=3,byrow=TRUE),widths=c(4,4,4),height=c(4,4,4))

plot(c(0,1),type='n',xlab="",ylab="",axes=FALSE)
title("(a) Experimental design",cex.main=cm)



par(mar=c(4.5,5.5,2.5,1.5))
plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,1],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phC_sub[,1],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.031",side=3,adj=.1,line=-2)
title("(b) 15d weight, C",cex.main=cm)
box()

plot(jitter(gen[host_trt==2]-1),gemma_phC_sub[,2],pch=19,col=csC[gen[host_trt==2]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phC_sub[,2],INDEX=gen[host_trt==2],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.002",side=3,adj=.1,line=-2)
title("(c) 21d weight, C",cex.main=cm)
box()

surv<-tapply(INDEX=gen[host_trt==2]-1,gemma_phC_sub[,3],mean)
barplot(surv,col=csC,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca,cex=cx)
mtext("P = 0.031",side=3,adj=.5,line=-2)
title("(d) Survival, C",cex.main=cm)

plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,1],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 15-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phRW_sub[,1],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.034",side=3,adj=.1,line=-2)
title("(e) 15d weight, RW",cex.main=cm)
box()


plot(jitter(gen[host_trt==1]-1),gemma_phRW_sub[,2],pch=19,col=csRW[gen[host_trt==1]],axes=FALSE,xlab="Genotype",ylab="Resid. 21-d weight",cex.lab=cl,cex.axis=ca,cex=cx)
mns<-tapply(X=gemma_phRW_sub[,2],INDEX=gen[host_trt==1],mean,na.rm=TRUE)
lines(c(-.2,.2),rep(mns[1],2),lwd=lx)
lines(c(.8,1.2),rep(mns[2],2),lwd=lx)
lines(c(1.8,2.2),rep(mns[3],2),lwd=lx)
axis(1,at=c(0,1,2))
axis(2)
mtext("P = 0.418",side=3,adj=.1,line=-2)
title("(f) 21d weight, RW",cex.main=cm)
box()


surv<-tapply(INDEX=gen[host_trt==1]-1,gemma_phRW_sub[,3],mean)
barplot(surv,col=csRW,xlab="Genotype",ylab="Survival proportion",cex.lab=cl,cex.axis=ca,cex=cx)
mtext("P = 0.754",side=3,adj=.5,line=-2)
title("(g) Survival, RW",cex.main=cm)
dev.off()


###################################################
dat<-dat[1:138,]
host_source_RW<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli")])
host_source_C<-as.character(dat$Host[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli")])
host_source<-c(host_source_RW,host_source_C)

dat_sub_C<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "C" & dat$Species == "knulli"),]
dat_sub_RW<-dat[(dat$Population != "BCTURN" & dat$Treatment.March.19 == "RW" & dat$Species == "knulli"),]


###################################################

tapply(X=gen==1,INDEX=host_source,mean)
#        C        RW 
#0.1052632 0.6800000 
tapply(X=gen==2,INDEX=host_source,mean)
#        C        RW 
#0.4605263 0.3200000 
tapply(X=gen==3,INDEX=host_source,mean)
#        C        RW 
#0.4342105 0.0000000 

