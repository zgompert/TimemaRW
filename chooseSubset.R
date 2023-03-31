# chosing subset with maximum data

dat<-read.table("miss",header=FALSE)
## cali, has lowest coverage, go with 350 missing, gives us 8
cali<-dat[ which(1:29 %in% grep(x=dat[,1],pattern="cali") & dat[1:29,2] < 350),1]
write.table(cali,"sub_cali.filelist",quote=FALSE,row.names=FALSE,col.names=FALSE)

## popp, < 88 gives us 8
popp<-dat[c(237:266)[which(237:266 %in% grep(x=dat[,1],pattern="popp") & dat[237:266,2] < 88)],1]
write.table(popp,"sub_popp.filelist",quote=FALSE,row.names=FALSE,col.names=FALSE)

## petita, 10 < 79
xx<-168:236
petita<-dat[xx[sample(which(dat[xx,2] < 79),10,replace=FALSE)],1]
write.table(petita,"sub_petita.filelist",quote=FALSE,row.names=FALSE,col.names=FALSE)
