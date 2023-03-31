## summarize posteriors
lf<-list.files(pattern="og.log",recursive=TRUE,include.dirs=TRUE)

pp<-vector("list",6)
for(i in 1:6){
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


d_knulli<-density(pp_knulli,adj=1.3)
d_rw<-density(pp_rw,adj=1.3)

library(scales)
pdf("DivTimePost.pdf",width=6,height=5)
par(mar=c(5,5,1,1))
plot(d_knulli$x,d_knulli$y,type='n',ylim=c(0,0.2),xlab="Time (MYA)",ylab="Posterior density",cex.lab=1.5,cex.axis=1.2)
polygon(c(d_knulli$x,rev(d_knulli$x)),c(d_knulli$y,rep(0,length(d_knulli$y))),col="gray30",border=NA)
polygon(c(d_rw$x,rev(d_rw$x)),c(d_rw$y,rep(0,length(d_rw$y))),col=alpha("firebrick3",.6),border=NA)
legend(11,0.2,c("T. knulli RW, T. knulli C","T. knulli RW, T. poppensis"),fill=c("gray30",alpha("firebrick3",.6)))
dev.off()
