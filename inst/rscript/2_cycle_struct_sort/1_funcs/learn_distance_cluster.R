
# from input_stcid_topology
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
# 按距离聚类
# Clustering with a distance matrix 根据距离矩阵进行聚类
# 详见：https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
# 距离矩阵


##################################################################################
#################################### learn #######################################
#https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
### build_coor_by_distance
tmp <- matrix(c(0,20,20,20,40,60,60,60,100,120,120,120,
                20,0,20,20,60,80,80,80,120,140,140,140,
                20,20,0,20,60,80,80,80,120,140,140,140,
                20,20,20,0,60,80,80,80,120,140,140,140,
                40,60,60,60,0,20,20,20,60,80,80,80,
                60,80,80,80,20,0,20,20,40,60,60,60,
                60,80,80,80,20,20,0,20,60,80,80,80,
                60,80,80,80,20,20,20,0,60,80,80,80,
                100,120,120,120,60,40,60,60,0,20,20,20,
                120,140,140,140,80,60,80,80,20,0,20,20,
                120,140,140,140,80,60,80,80,20,20,0,20,
                120,140,140,140,80,60,80,80,20,20,20,0),
              nr=12, dimnames=list(LETTERS[1:12], LETTERS[1:12]))
d <- as.dist(tmp)
mds.coor <- cmdscale(d)
plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
     rownames(mds.coor), cex=0.8)
abline(h=0,v=0,col="gray75")

##################################################################################
## 1.《hierarchical clustering》(hclust)
clust_res = hclust(dist(tmp), method="single")
# 层次聚类 hierarchical clustering
plot(hclust(dist(tmp), method="single"))


##################################################################################
## 2.《k-medoids clustering algorithm》(clara)
### clara
#https://www.rdocumentation.org/packages/cluster/versions/2.1.2/topics/clara
#install.packages("cluster")
library("cluster")
x <- rbind(cbind(rnorm(200,0,8), rnorm(200,0,8)),
           cbind(rnorm(300,50,8), rnorm(300,50,8)))
clarax <- clara(x, 2, samples=50)
clarax
clarax$clusinfo

### plot
mds.coor = x
rownames(mds.coor) = c(0:499)
plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
     rownames(mds.coor), cex=0.8)
abline(h=0,v=0,col="gray75")


#test
clarax <- clara(coor_by_distance_matrix, 5)
clus_cluster = clarax$clustering
plot_cluster_with_color(coor_by_distance_matrix, clus_cluster)

##################################################################################
## 3.《partitioning around medoids》(pam)
#https://www.rdocumentation.org/packages/cluster/versions/2.1.2/topics/pam
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)),
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5)))
pamx <- pam(x, 2)
pamx # Medoids: '7' and '25' ...
summary(pamx)
plot(pamx)

### plot
mds.coor = x
rownames(mds.coor) = c(0:24)
plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
     rownames(mds.coor), cex=0.8)
abline(h=0,v=0,col="gray75")



##################################################################################
# https://datascience.stackexchange.com/questions/761/clustering-geo-location-coordinates-lat-long-pairs
## 4. 《DBSCAN》
#https://www.rdocumentation.org/packages/dbscan/versions/1.1-8/topics/dbscan


##################################################################################
## 5. 《OPTICS》








