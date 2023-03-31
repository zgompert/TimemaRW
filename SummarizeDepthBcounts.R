library(data.table)

dat<-fread("~/../gompert-group3/data/timema_clines_rw_SV/align_rw_plus/t_rw_knulli_perform_reg_depth.txt",header=FALSE)
dm<-as.matrix(dat[,-c(1:2)])
mn<-apply(dm,1,mean)
prop<-apply(dm > 0,1,mean)

keep<-mn > 2 & prop > 0.8

pos<-as.numeric(unlist(dat[keep,2]))
length(pos)
#[1] 160703 = number of sites


save(list=ls(),file="positions.rdat")

## scan all initial called SNPs on scaffold 500
snps_all<-scan("prefilt_var_pos.txt")
## scan snps remaining
snps_kept<-scan("filt_var_pos.txt")

## remove true snps from initial set
drop<-which(snps_all %in% snps_kept)

## get initially called but then dropped, also includes outside of Perform
bad_snps <-snps_all[-drop]

pos_keep<-pos[-which(pos %in% bad_snps)]
## 157212 = number of sites that wouldn't have been filtered out

