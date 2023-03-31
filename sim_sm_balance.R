## simulation selection gene flow balance for a pair of populations with different directional selection or different het advantage

## libraries
library(RColorBrewer)

## dp from gene flow, a is recipient, b is donor; m is gene flow prop.
dp_gf<-function(m=0,pa=NA,pb=NA){
    pf<-pa+m*(pb-pa)
    return(pf)
    }
    
## dp from directional selection
## w11 = 1+s, w12 = 1+sh, w22 = 1
dp_ds<-function(p0=NA,h=0.5,s=0){
    pf<-p0+s*p0*(1-p0)*(p0+h*(1-2*p0))
    return(pf)
    }
    
## dp from balancing selectino
## w11 = 1-s1, w12 = 1, w22 = 1-s2
dp_bs<-function(p0,s1=0,s2=0){
    pf<-p0+p0*(1-p0)*(s2 - p0*(s1+s2))
    }
    

## direction selection version    
s<-seq(0,.5,0.01)
m<-seq(0,.1,.01)
Ns<-length(s)
Nm<-length(m)
P_ds_gf_1<-matrix(NA,nrow=Ns,ncol=Nm)
P_ds_gf_2<-matrix(NA,nrow=Ns,ncol=Nm)
for(i in 1:Ns){for(j in 1:Nm){
    pp<-matrix(.5,nrow=1000,ncol=2)
    for(k in 2:1000){
        pt1<-dp_ds(pp[k-1,1],h=.5,s=s[i])
        pt2<-1-dp_ds(1-pp[k-1,2],h=.5,s=s[i])
        pp[k,1]<-dp_gf(m[j],pt1,pt2)
        pp[k,2]<-dp_gf(m[j],pt2,pt1)
    }
    P_ds_gf_1[i,j]<-pp[1000,1]
    P_ds_gf_2[i,j]<-pp[1000,2]
}  }  

cs<-brewer.pal(9,"RdBu")
brks<-seq(-.51,.51,length.out=10)
image(abs(P_ds_gf_1-P_ds_gf_2)-.5,axes=FALSE,col=cs,breaks=brks)
axis(1,at=(1:Ns-1)/(Ns-1),s,las=2)
axis(2,at=(1:Nm-1)/(Nm-1),m,las=2)
box()

##
