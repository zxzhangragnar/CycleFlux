# ---【其它方法找离群点】
# 离群点检测：
# https://towardsdatascience.com/a-brief-overview-of-outlier-detection-techniques-1e0b2c19e561?gi=3fff9c375e17
# 
# 在R语言中的离群点检测：
# https://www.r-bloggers.com/2016/12/outlier-detection-and-treatment-with-r/
#   https://statsandr.com/blog/outliers-detection-in-r/
#
# 在R语言中的“outlier”包：(1维数据)
# https://www.rdocumentation.org/packages/outliers/versions/0.14/topics/outlier
# https://www.rdocumentation.org/packages/outliers/versions/0.14/topics/rm.outlier
# 
# 
# 在R语言中的“HDoutliers”包：Leland Wilkinson 检测多维异常值的算法 (多维数据)
# https://rdrr.io/cran/HDoutliers/
# https://www.rdocumentation.org/packages/HDoutliers/versions/1.0.3/topics/HDoutliers
#



#test
distance_matrix = get_distance_matrix(652)
coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
plot(coor_by_distance_matrix)

#suspicious_idv_struct
outlier_coor_by_distance_matrix <- HDoutliers(coor_by_distance_matrix, alpha=.62, transform = TRUE)

outlier_coor_by_distance_matrix
length(outlier_coor_by_distance_matrix)
plotHDoutliers(coor_by_distance_matrix, outlier_coor_by_distance_matrix)




##################################################################################
#################################### learn #######################################
################### HDoutliers Test
#在R语言中的“HDoutliers”包：Leland Wilkinson 检测多维异常值的算法 (多维数据)
#https://www.rdocumentation.org/packages/HDoutliers/versions/1.0.3/topics/HDoutliers
install.packages("HDoutliers")
library(HDoutliers)

data(ex2D)
out.ex2D <- HDoutliers(ex2D)
plotHDoutliers(ex2D,out.ex2D)






