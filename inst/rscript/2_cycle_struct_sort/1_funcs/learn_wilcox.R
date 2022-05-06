#
# reference
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test
# https://zhuanlan.zhihu.com/p/150760480
#

##
# x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
# y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
# wilcox.test(x, y, paired = TRUE, alternative = "greater")





women_weight <- c(88.9, 81.2, 73.3, 21.8, 63.4, 84.6, 28.4, 28.8, 28.5)
men_weight <- c(37.8, 80, 33.4, 36, 89.4, 83.3, 97.3, 81.3, 92.4)
# 建立一个数据框
my_data <- data.frame( 
  group = rep(c("Woman", "Man"), each = 9),
  weight = c(women_weight,  men_weight)
)

#install.packages("dplyr")
library(dplyr)


group_by(my_data, group) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE)
  )


# 1.初步检验两独立样本是否满足正态分布

# 假设1：两个样本是否独立？
# 是的，因为来自男性和女性的样本无关。
# 假设2：两组中每组的数据是否服从正态分布？
# 我们将使用with()和shapiro.test()的函数来为每组样本计算Shapiro-Wilk测试。

# Shapiro-Wilk normality test for Men's weights
with(my_data, shapiro.test(weight[group == "Man"]))# p = 0.017
# Shapiro-Wilk normality test for Women's weights
with(my_data, shapiro.test(weight[group == "Woman"])) # p = 0.045

# 输出结果中，两个p值小于显着性水平0.05，说明两组数据的分布与正态分布有显着差异。数据分布不符合正态分布的假设检验成立。
# 请注意，如果两组数据中只有一组不是正态分布，也要使用非参数两样本Wilcoxon秩检验。
# 假设3：这两个总体是否符合方差齐性？
# 我们将使用F检验来检验方差齐性。可以使用var.test()函数执行以下操作：

res.ftest <- var.test(weight ~ group, data = my_data)
res.ftest


# data:  weight by group
# F = 0.88062, num df = 8, denom df = 8, p-value = 0.8617


#
# F检验为p = 0.8617。它大于显着性水平alpha = 0.05。因此，两组数据的方差之间没有显著差异。因此我们认为男女两组方差相等（方差齐性）。
# 由于以上3个假设综合，由于数据不符合正态分布，因此，我们不可以使用student-t检验。需要使用两独立样本Wilcoxon检验。
# 



# 2.计算两独立样本Wilcoxon检验

# 问题：男女体重之间有显着差异吗？

res <- wilcox.test(weight ~ group, data = my_data, var.equal = TRUE)
res <- wilcox.test(weight ~ group, data = my_data)
res
# var.equal：一个逻辑变量，指示是否将两个方差视为相等。
# 如果为 TRUE，则使用合并方差来估计方差，否则使用 Welch 检验。

# W = 59, p-value = 0.1135
# 在上面的结果中：
# p值是wilcoxon检验的显着性水平（p值= 0.1135）。
# 如果要检验男性的体重是否小于女性的体重，请输入以下内容：
wilcoxon.test(weight ~ group, data = my_data,
              var.equal = TRUE, alternative = "less")
# 或者，如果您想测试男性的体重是否大于女性的体重，请输入
wilcoxon.test(weight ~ group, data = my_data,
              var.equal = TRUE, alternative = "greater")


# 3.结果解释
# 检验的p值为 0.1135，大于显着性水平alpha = 0.05
# 我们可以得出结论，男性的体重与女性的体重没有显著不同。


