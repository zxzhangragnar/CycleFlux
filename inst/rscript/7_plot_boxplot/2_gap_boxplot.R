## 此文件中的内容并非 主流程中的内容 可以不执行


###############################################
#1.看有gap的环,除了gap以外的rct是否为up 
#2.重新做boxplot等图分析(用cycle_flux_list)
#3.使用genecards.org更新缺失的基因value
###############################################





#####################################################################################
#####################################################################################
##1.foldchange

#cycle_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_final/cycle_flux_list.RData")

#预处理
#ifgap
for(i in 1:length(cycle_flux_list)){
  cycle_flux_list[[i]][which(cycle_flux_list[[i]]$ifgap=="gap"),"ifgap"] = as.integer(1)
  cycle_flux_list[[i]][which(cycle_flux_list[[i]]$ifgap=="normal"),"ifgap"] = as.integer(0)
}
#ifup
for(i in 1:length(cycle_flux_list)){
  cycle_flux_list[[i]][which(cycle_flux_list[[i]]$ifup=="up"),"ifup"] = as.integer(1)
  cycle_flux_list[[i]][which(cycle_flux_list[[i]]$ifup=="normal"),"ifup"] = as.integer(0)
}


###################################################################
##################### for 248 cycle ###############################
pdf("foldchangesign_stc_gap1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  boxplot(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nongap","WL-gap","HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()



##################### abs for 248 cycle ###############################
#取绝对值做boxplot
pdf("foldchangesign_stc_gap2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  boxplot(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nongap","WL-gap","HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()




#######################################################################################
#同时具有up现象的环怎么样

pdf("foldchangesign_stc_up1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonup","WL-up","HF-noup","HF-up","isolated-noup","isolated-up","others-noup","others-up"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()

#abs
pdf("foldchangesign_stc_up2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonup","WL-up","HF-noup","HF-up","isolated-noup","isolated-up","others-noup","others-up"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()


#######################################################################################
#同时具有gap up现象的环怎么样
#ug

#预处理
for(i in 1:length(cycle_flux_list)){
  cycle_flux_list[[i]][,"ifug"] = as.integer(0)
}
for(i in 1:length(cycle_flux_list)){
  cycle_flux_list[[i]][which((cycle_flux_list[[i]]$ifup==1) & (cycle_flux_list[[i]]$ifgap==1)),"ifug"] = as.integer(1)
}



pdf("foldchangesign_stc_ug1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonug","WL-ug","HF-noug","HF-ug","isolated-noug","isolated-ug","others-noug","others-ug"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()

#abs
pdf("foldchangesign_stc_ug2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonug","WL-ug","HF-noug","HF-ug","isolated-noug","isolated-ug","others-noug","others-ug"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()







#######################################################################################
#分析gap和up的关系 不考虑struct_id

pdf("foldchangesign_gap_up1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),las=2,names=c("nogap_down","nogap_up","gap_down","gap_up"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_flux_list[[i]][,"updownsign"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()

pdf("foldchangesign_gap_up2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),las=2,names=c("nogap_down","nogap_up","gap_down","gap_up"),xlab="",ylab="updownsign",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(abs(cycle_flux_list[[i]][,"updownsign"])~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()

#####################################################################################
#####################################################################################
##2.mutation

#############################
##
#gap

# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")

pdf("mutation_stc_gap1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  boxplot(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nongap","WL-gap","HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="mutation level",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()


#############################
##
#up

# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")


pdf("mutation_stc_up1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  boxplot(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonup","WL-up","HF-noup","HF-up","isolated-noup","isolated-up","others-noup","others-up"),xlab="",ylab="mutation level",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifup"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()

#############################
##
#ug

# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")


pdf("mutation_stc_ug1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  boxplot(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifug"])+cycle_flux_list[[i]][,"struct_id"]*10),las=2,names=c("WL-nonug","WL-ug","HF-noug","HF-ug","isolated-noug","isolated-ug","others-noug","others-ug"),xlab="",ylab="mutation level",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifug"])+cycle_flux_list[[i]][,"struct_id"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()



#######################################################################################
#分析gap和up的关系 不考虑struct_id

pdf("mutation_gap_up1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(cycle_flux_list)){
  print(i)
  boxplot(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),las=2,names=c("nogap_down","nogap_up","gap_down","gap_up"),xlab="",ylab="mx_mutation_rate",col="lightblue",main=names(cycle_flux_list)[i])
  stripchart(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(as.integer(cycle_flux_list[[i]][,"ifgap"])*10+as.integer(cycle_flux_list[[i]][,"ifup"])),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()









