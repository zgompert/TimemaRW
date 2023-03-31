## summarize ABC output
## running R 4.0.2
library("abc")
library("data.table")

## read abc sims, 3 sets of priors
sims<-vector("list",3)
ff<-list.files(pattern="cout") ## est, hi, wide
for(i in 1:3){
	sims[[i]]<-as.matrix(fread(ff[i],header=FALSE))
}

#               C        RW
#BCE    0.3529412 0.8333333
#BCSH          NA 1.0000000
#BCTURN 0.2432432        NA
#BCXD   0.1875000        NA
obs<-c(0.2432,0.8333,0.3529)

## model post. probs
mo<-vector("list",3)
for(i in 1:3){
	mo[[i]]<-abc(target=obs,param=sims[[i]][,1],sumstat=sims[[i]][,16:18],method="rejection",tol=4e-05)
}

pp<-matrix(NA,nrow=3,ncol=4)
for(i in 1:3){
	pp[i,]<-table(mo[[i]]$unadj.values)/length(mo[[i]]$unadj.values)
}
#0 = both directional, 1 = bs on RW, 2 = bs on C, 3 = both bs
#           [,1]       [,2]      [,3]      [,4]
#[1,] 0.08891109 0.12587413 0.3406593 0.4445554
#[2,] 0.02597403 0.06293706 0.2797203 0.6313686
#[3,] 0.02297702 0.07992008 0.2937063 0.6033966

pp[,2]+pp[,4]
#[1] 0.5704296 0.6943057 0.6833167
pp[,3]+pp[,4]
#[1] 0.7852148 0.9110889 0.8971029

for(i in 1:3){
	bs<-which(sims[[i]][,1]==3)## model 3, best model
	mo[[i]]<-abc(target=obs,param=sims[[i]][bs,12:15],sumstat=sims[[i]][bs,16:18],transf=rep("log",4),method="ridge",tol=1.6e-04)
}

## estiamtes of s1 and s2 under bs
## est Ne values
# 1001 post. samples
#                          V12    V13    V14    V15
#Min.:                  0.0498 0.0004 0.0007 0.0474
#Weighted 2.5 % Perc.:  0.1362 0.0011 0.0009 0.1898
#Weighted Median:       0.3970 0.0248 0.0086 0.5841
#Weighted Mean:         0.4227 0.0662 0.0265 0.5592
#Weighted Mode:         0.2763 0.0090 0.0045 0.6397
#Weighted 97.5 % Perc.: 0.8387 0.3692 0.1545 0.8794
#Max.:                  0.9935 0.6804 0.3825 0.9569

## hi Ne values
#                          V12    V13    V14    V15
#Min.:                  0.0082 0.0008 0.0007 0.0150
#Weighted 2.5 % Perc.:  0.0147 0.0010 0.0010 0.0326
#Weighted Median:       0.0661 0.0163 0.0086 0.2336
#Weighted Mean:         0.1419 0.0526 0.0266 0.2848
#Weighted Mode:         0.0417 0.0069 0.0046 0.1423
#Weighted 97.5 % Perc.: 0.6859 0.3229 0.1471 0.8179
#Max.:                  1.6055 0.4597 0.4068 0.9412

## wide Ne values
#                          V12    V13    V14    V15
#Min.:                  0.0115 0.0007 0.0009 0.0153
#Weighted 2.5 % Perc.:  0.0248 0.0012 0.0010 0.0700
#Weighted Median:       0.1001 0.0232 0.0086 0.2985
#Weighted Mean:         0.1831 0.0625 0.0244 0.3563
#Weighted Mode:         0.0632 0.0098 0.0037 0.2222
#Weighted 97.5 % Perc.: 0.7590 0.3295 0.1445 0.8164
#Max.:                  1.6824 0.7328 0.3354 0.9229


pdf("tknulliSelEst.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
plot(1-pmin(mo[[3]]$adj.values[,1],1),1-pmin(mo[[3]]$adj.values[,2],1),pch=19,col=alpha("blue",.4),xlim=c(0,1),ylim=c(0,1),xlab="Fitness for RW homozygote",,ylab="Fitness for C homozygote",cex.lab=1.5,cex.axis=1.1)
o<-kde2d(1-pmin(mo[[3]]$adj.values[,1],1),1-pmin(mo[[3]]$adj.values[,2],1))
contour(o,add=TRUE,lwd=2)

points(1-pmin(mo[[3]]$adj.values[,3],1),1-pmin(mo[[3]]$adj.values[,4],1),pch=19,col=alpha("indianred",.4))
o<-kde2d(1-pmin(mo[[3]]$adj.values[,3],1),1-pmin(mo[[3]]$adj.values[,4],1))
contour(o,add=TRUE,lwd=2)
legend(0,.22,c("C","RW"),col=c("blue","indianred"),bty='n')
     
dev.off()

