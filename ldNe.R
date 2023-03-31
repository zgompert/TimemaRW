## burrow's delta for a two SNPs
compDelta<-function(g1=NA,g2=NA){ ## genotype vectors for twl loci
    nAABB<-sum(g1==2 & g2==2)
    nAABb<-sum(g1==2 & g2==1)
    nAaBB<-sum(g1==1 & g2==2)
    nAaBb<-sum(g1==1 & g2==1)
    
    pA<-mean(g1)/2
    pB<-mean(g2)/2
    
    nN<-length(g1)
    ## Notation from Zaykin 2004 following Weir 1979
    DAB<-1/nN * (2*nAABB + nAABb + nAaBB + 0.5 * nAaBb) - 2 * pA * pB
    
    ## unbiased
    DAB = DAB * (nN/(nN-1))

    ## Eqn 1 Waples and Do 2008
    DA<-(mean(g1==2) - pA^2)
    DB<-(mean(g2==2) - pB^2)
    rAB<-DAB/sqrt((pA * (1-pA) + DA) * (pB*(1-pB) + DB))
    return(rAB)
}

## use with fewer than 30 individuals
NeSmallS<-function(S=NA,r2=NA){

    er2<-0.0018 + 0.907/S + 4.44/S^2
    r2prime<-r2 - er2
    Ne<-(0.308 + sqrt(0.308^2-2.08*r2prime))/(2*r2prime)
    return(Ne)
}

## use with more than 30 individuals
NeLargeS<-function(S=NA,r2=NA){

    er2<-1/S + 3.19/S^2
    r2prime<-r2 - er2
    Ne<-(1/3+ sqrt(1/9-2.76*r2prime))/(2*r2prime)
    return(Ne)
}
