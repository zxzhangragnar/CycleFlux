
##################################### part5 mutation analysis plot #####################################
###############################################################################################



###################################################################
##################### isgap ###################################
# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")


pdf("mutation_stc_gap_248_a1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  boxplot(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),las=2,names=c("WL-nongap","WL-gap","HF-nogap","HF-gap","isolated-nogap","isolated-gap","others-nogap","others-gap"),xlab="",ylab="mutation level",col="lightblue",main=names(tumor_cyc_foldchangesign_merged)[i])
  stripchart(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]~c(tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]+tumor_cyc_foldchangesign_merged[[i]][,"stcid"]*10),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
  
}
dev.off()









####################################################################
##################### foldchangesign ###################################

# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")

pdf("mutation_foldchangesign_a1.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  boxplot(as.numeric(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")])~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*100),las=2,xlab="mx_mutation_rate",ylab="fc_updown_sign",col="lightblue",main=names(tumor_cyc_foldchangesign_merged)[i])
  stripchart(as.numeric(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")])~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*100),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
}
dev.off()



####################################################################
##################### line foldchangesign ###################################

# # cycle_mx_mutation_rate
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")
# # tumor_cyc_foldchangesign_merged.RData
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")
# 
# pdf("mutation_foldchangesign_b1.pdf")
# par(mar=c(10,5,5,5))
# for(i in 1:length(tumor_cyc_foldchangesign_merged)){
#   plot(as.numeric(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")])~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*20),las=2,xlab="mx_mutation_rate",ylab="fc_updown_sign",main=names(tumor_cyc_foldchangesign_merged)[i])
#   abline(h=0,col=2,lty=2)
# }
# dev.off()


####################################################################
##################### abs(foldchangesign) ###################################

# cycle_mx_mutation_rate
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")

pdf("mutation_foldchangesign_a2.pdf")
par(mar=c(10,5,5,5))
for(i in 1:length(tumor_cyc_foldchangesign_merged)){
  boxplot(as.numeric(abs(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")]))~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*100),las=2,xlab="mx_mutation_rate",ylab="fc_updown_sign",col="lightblue",main=names(tumor_cyc_foldchangesign_merged)[i])
  stripchart(as.numeric(abs(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")]))~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*100),
             method = "jitter",
             pch = 19,cex=0.5,
             vertical = TRUE,
             add = TRUE)
  abline(h=0,col=2,lty=2)
}
dev.off()



####################################################################
##################### line abs(foldchangesign) ###################################

# 
# # cycle_mx_mutation_rate
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_Analysis/Analysis_mutation_D/cycle_mx_mutation_rate.RData")
# # tumor_cyc_foldchangesign_merged.RData
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")
# 
# pdf("mutation_foldchangesign_b2.pdf")
# par(mar=c(10,5,5,5))
# for(i in 1:length(tumor_cyc_foldchangesign_merged)){
#   plot(as.numeric(abs(tumor_cyc_foldchangesign_merged[[i]][,c("fc_updown_sign")]))~c(cycle_mx_mutation_rate[[i]][,"mx_mutation_rate"]+tumor_cyc_foldchangesign_merged[[i]][,"fc_gap_judge"]*20),las=2,xlab="mx_mutation_rate",ylab="fc_updown_sign",main=names(tumor_cyc_foldchangesign_merged)[i])
#   abline(h=0,col=2,lty=2)
# }
# dev.off()























