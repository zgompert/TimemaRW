library(data.table)

## read synteny dat
dat<-fread("/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/comp_aligns/out_synteny_knulli.psl",header=FALSE)
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


invs<-read.table("~/../gompert-group3/data/Tknulli_nanopore/alignment_BCEC-22-4/LG11_inversions.txt",header=FALSE)
invs_st<-invs[,2]
invs_len<-as.numeric(gsub(pattern="SVLEN=",x=invs[,10],replacement=""))
keep<-which(invs_len>1000000)

big_st<-invs_st[keep]
big_end<-invs_st[keep]+invs_len[keep]
i<-11
tcr<-grep(x=dfdat[,14],pattern=chtab[i,1]) 
tkn<-grep(x=dfdat[,10],pattern=chtab[i,2])
cc<-tcr[tcr %in% tkn]
subd<-dfdat[cc,]
xub<-max(subd[,13]);yub<-max(subd[,17]) 
bnds_inv<-c(13650000,bnds[2])


library(ape);library(phytools);library(phangorn)
## read in the trees
trs1<-read.nexus(file="o_knulli_1/sub_perform_comb_og.trees")

## midpoint root
rtct<-midpoint(trs1)

## consensus tree, 50 pp or >
#ct<-consensus(rtct,.5) ## does this work? or just from trs1
ct<-consensus(trs1,p=.5) ## does this work? or just from trs1


## get labels by species/type
nms<-ct$tip.label
nms<-gsub(pattern="timemaN",replacement="",nms)

taxa<-rep(1,length(nms))

bcecc<-read.table("bce_cc.filelist")
taxa[which(nms %in% bcecc[,1])]<-2
bcerw<-read.table("bce_rw.filelist")
taxa[which(nms %in% bcerw[,1])]<-3
taxa[1:8]<-4
taxa[72:80]<-5

ct$node.label<-round(ct$node.label,2)


## summarize posteriors
lf<-list.files(pattern="og.log",recursive=TRUE,include.dirs=TRUE)
lf<-lf[c(7:12,1:6)]

pp<-vector("list",12)
for(i in 1:12){
	pp[[i]]<-read.table(lf[i],header=TRUE,comment.char="#")
}

pp_knulli<-c(pp[[1]]$mrca.age.Tknulli.[-c(1:1000)],pp[[2]]$mrca.age.Tknulli.[-c(1:1000)],pp[[3]]$mrca.age.Tknulli.[-c(1:1000)])
summary(pp_knulli)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6712  5.5234  7.4602  7.8307  9.7199 23.7754 
quantile(pp_knulli,probs=c(0.5,0.05,0.95))
#      50%        5%       95% 
# 7.460225  3.426531 13.452333 

pp_rw<-c(pp[[4]]$mrca.age.TRW.[-c(1:1000)],pp[[5]]$mrca.age.TRW.[-c(1:1000)],pp[[6]]$mrca.age.TRW.[-c(1:1000)])
plot(density(pp_rw))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6433  3.3833  4.7312  5.1862  6.5187 21.8673
quantile(pp_rw,probs=c(0.5,0.05,0.95))
#     50%       5%      95% 
#4.731152 2.050723 9.835396

## TMRCA for knulli and poppensis from Lindtke et al. 2017 Table S10
# shape = 10.0184
# scale = 0.4256

pp_tktp<-rgamma(length(pp_rw),shape=10.0184,scale=0.4256)
mean(pp_rw > pp_tktp)
#[1] 0.6011184


## TMRCA for each inversion haplotype within knulli
pp_knrw<-c(pp[[10]]$mrca.age.TknRW.[-c(1:1000)],pp[[11]]$mrca.age.TknRW.[-c(1:1000)],pp[[12]]$mrca.age.TknRW.[-c(1:1000)])
quantile(pp_knrw,probs=c(0.5,0.05,0.95))
#     50%       5%      95% 
#3.369245 1.477328 7.485583 

pp_knc<-c(pp[[7]]$mrca.age.TknC.[-c(1:1000)],pp[[8]]$mrca.age.TknC.[-c(1:1000)],pp[[9]]$mrca.age.TknC.[-c(1:1000)])
quantile(pp_knc,probs=c(0.5,0.05,0.95))
#     50%       5%      95% 
#2.521746 1.149609 5.353409 

pps<-cbind(pp_knc,pp_knrw,pp_knulli,pp_rw,pp_tktp)
est<-apply(pps,2,quantile,probs=c(.5,.05,.95))


## plots
save(list=ls(),file="beast.rdat")
cs<-c("darkgray","cadetblue","firebrick","gray40","forestgreen")

pdf("fig_treePlus.pdf",width=12,height=5)
par(mfrow=c(1,3))
par(mar=c(4.5,2,2.5,1))
cl<-1.4;cm<-1.4;ca<-1.1

plot(bnds,c(1,1),xlim=c(0,xub),ylim=c(.5,3.5),type='l',lwd=4,axes=FALSE,xlab="Position (bp)",ylab="",cex.lab=1.2,cex.lab=cl)
lines(bnds_inv,c(2,2),lwd=4)
lines(c(big_st[1],big_end[1]),c(3,3),lwd=4)
text(rep(3e7,3),.3+1:3,c("Perform locus from PCA of GBS","Inversion between species from genome alignment","Inversion within T. knulli from nanopore reads"),cex=1.1)
axis(1,cex.axis=.9,cex.axis=ca)
title(main="(D) Evidence for T. knulli inversion",cex.main=cm)

#par(mar=c(1,1,1,1))
plot(ct,show.tip.label=FALSE)
tiplabels(pch=19,col=cs[taxa],offset=1)
legend(-1,15,c("T. knulli RW","T. knulli C","T. poppensis","T. petita","T. californicum"),pch=19,col=cs[c(3,2,5,4,1)],bty='n')
nodelabels(pie=ct$node.label,piecol=c("black","white"),cex=.25)
title(main="(E) Tree for Perform",cex.main=cm)


#par(mar=c(5,5,1,1))
dotchart(x=est[1,],pch=19,xlim=c(0,15),labels=c("T. knulli C","T. knulli RW","T. knulli RW + C","T. knulli RW + T. poppensis","T. knulli + T. poppensis"),xlab="Time (MYA)",cex.lab=cl,cex.axis=ca)
for(j in 1:5){
	lines(c(est[2,j],est[3,j]),rep(j,2))
}
title(main="(F) Divergence times",cex.main=cm)
dev.off()
