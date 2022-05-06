
##################################### part4 analysis plot #####################################
###############################################################################################


#####################################################################################
#####################################################################################

load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")

###################################################################
##################### for 248 cycle ###############################
pdf("foldchangesign_stc_gap_248_b1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  boxplot(tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),las=2,names=c("WL-nongap","WL-gap",
                                                                                                                                                                     "HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="fc_updown_sign",col="lightblue",main=names(tumor_cyc_foldchangesign_merged)[i])
  stripchart(tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()




##################### abs for 248 cycle ###############################
pdf("foldchangesign_stc_gap_248_b2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  boxplot(abs(tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"])~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),las=2,names=c("WL-nongap","WL-gap",
                                                                                                                                                                                           "HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="fc_updown_sign",col="lightblue",main=names(tumor_cyc_foldchangesign_merged)[i])
  stripchart(abs(tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"])~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()




##################### negative positive ###############################
positive_tumor_cyc_foldchangesign_merged = list()
negative_tumor_cyc_foldchangesign_merged = list()

for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  positive_tumor_cyc_foldchangesign_merged[[i]] = tumor_cyc_foldchangesign_merged[[i]][which(tumor_cyc_foldchangesign_merged[[i]]$fc_updown_sign>0),]
}

for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  negative_tumor_cyc_foldchangesign_merged[[i]] = tumor_cyc_foldchangesign_merged[[i]][which(tumor_cyc_foldchangesign_merged[[i]]$fc_updown_sign<0),]
}

##################### negative positive ###############################




##################### positive value for 248 cycle ###############################
pdf("foldchangesign_stc_gap_248_b3.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(positive_tumor_cyc_foldchangesign_merged)){
  boxplot(positive_tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(positive_tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+positive_tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),las=2,names=c("WL-nongap","WL-gap",
                                                                                                                                                                                                "HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="fc_updown_sign",col="lightblue",main=names(positive_tumor_cyc_foldchangesign_merged)[i])
  stripchart(positive_tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(positive_tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+positive_tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()



##################### negative value for 248 cycle ###############################
pdf("foldchangesign_stc_gap_248_b4.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(negative_tumor_cyc_foldchangesign_merged)){
  boxplot(negative_tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(negative_tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+negative_tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),las=2,names=c("WL-nongap","WL-gap",
                                                                                                                                                                                                "HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="fc_updown_sign",col="lightblue",main=names(negative_tumor_cyc_foldchangesign_merged)[i])
  stripchart(negative_tumor_cyc_foldchangesign_merged[[i]][,"fc_updown_sign"]~c(negative_tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+negative_tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()
































######   debug  ###########  
thca_table = as.data.frame(table(tumor_cyc_foldchangesign_merged[["THCA"]][,c("stcid", "isgap")]))









