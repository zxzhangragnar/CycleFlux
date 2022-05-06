
##################################### part6 heatmap analysis #####################################
###############################################################################################

# 用heatmap来表示得到的结果
# 如:用红色表示断了,蓝色表示没断
# 
# 直观的比较各种模型
# 看哪种模型断的环的数量适中 

####################################################################
##################### cyc gap updownsign ###########################
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")


# 测试数据
COAD_df = tumor_cyc_foldchangesign_merged[["COAD"]]


COAD_df_gap = COAD_df[,c("fc_gap_judge","fc_updown_sign")]

COAD_ma_gap = as.matrix(COAD_df_gap)

pdf("heatmap_cyc_gap_a1.pdf")
par(mar=c(10,5,5,5))
#代码1
require(graphics); require(grDevices)
x  <- COAD_ma_gap
colnames(x) = c("gap","fc")
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              main = "gap fc")

dev.off()













####################################################################
##################### cyc stcid updownsign ###########################
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")

COAD_df = tumor_cyc_foldchangesign_merged[["COAD"]]

pdf("heatmap_cyc_gap_stcid_a1.pdf")
par(mar=c(10,5,5,5))

for(i in 1:4){
  COAD_df_gap = COAD_df[which(COAD_df$stcid == i),c("fc_gap_judge","fc_updown_sign")]
  COAD_ma_gap = as.matrix(COAD_df_gap)
  require(graphics); require(grDevices)
  x  <- COAD_ma_gap
  colnames(x) = c("gap","fc")
  rc <- rainbow(nrow(x), start = 0, end = .3)
  cc <- rainbow(ncol(x), start = 0, end = .3)
  heatmap(x, col = cm.colors(256), scale = "column",
          RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
          main = "gap fc")
}

dev.off()




####################################################################
##################### cyc stcid updownsign ###########################
# tumor_cyc_foldchangesign_merged.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/7_plot_boxplot/plot_analysis/plot_codes/tumor_cyc_foldchangesign_merged.RData")

COAD_df = tumor_cyc_foldchangesign_merged[["COAD"]]

COAD_df_gap = COAD_df[,c("fc_gap_judge","fc_updown_sign","stcid")]

COAD_ma_gap = as.matrix(COAD_df_gap)

pdf("heatmap_cyc_gap_stcid_relation1.pdf")
#par(mar=c(10,5,5,5))
#代码1
require(graphics); require(grDevices)
x  <- COAD_ma_gap
colnames(x) = c("gap","fc","stcid")
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
heatmap(x, col = cm.colors(256), scale = "column",
        RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
        ylab =  "cycid",
        main = "stcid fc_updown_sign isgap")

dev.off()












