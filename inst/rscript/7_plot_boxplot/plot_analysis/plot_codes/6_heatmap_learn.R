
##################################### part6 heatmap analysis #####################################
###############################################################################################

# 用heatmap来表示得到的结果
# 如:用红色表示断了,蓝色表示没断
# 
# 直观的比较各种模型
# 看哪种模型断的环的数量适中 




# 测试数据
a=c(12,14,17,11,16)
b=c(4,20,15,11,9)
c=c(5,7,19,8,18)
d=c(15,13,11,17,16)
e=c(12,19,16,7,9)

A=cbind(a,b,c,d,e)
B=rbind(a,b,c,d,e)

pdf("heatmap_a1.pdf")
par(mar=c(10,5,5,5))
#代码1
require(graphics); require(grDevices)
x  <- as.matrix(A)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "heatmap(<BIOTREE>, ..., scale = \"column\")")
utils::str(hv)

dev.off()


























