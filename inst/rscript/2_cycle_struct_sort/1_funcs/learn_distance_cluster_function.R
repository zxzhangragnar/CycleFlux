
# from input_stcid_topology
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
# 按距离聚类
# Clustering with a distance matrix 根据距离矩阵进行聚类
# 详见：https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
# 距离矩阵

#########################################################################
# get coor_by_distance_matrix
get_distance_matrix <- function(cycle_num) {
  #cycle_num = 248
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
  cyc_all_distance_df = cyc_all_distance_df
  distance_scale <- cyc_all_distance_df$distance
  #距离矩阵
  #distance_matrix =matrix(distance_scale, nr=248, byrow=TRUE, dimnames=list(c(0:651), c(0:651)))
  distance_matrix =matrix(distance_scale, nr=cycle_num, byrow=TRUE, dimnames=list(c(0:(cycle_num-1)), c(0:(cycle_num-1))))

  return(distance_matrix)
}




#########################################################################
### build_coor_by_distance
#https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
build_coor_by_distance <- function(distance_matrix) {
  tmp = distance_matrix
  d <- as.dist(tmp)
  mds.coor <- cmdscale(d)
  plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
  text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
       rownames(mds.coor), cex=0.8)
  abline(h=0,v=0,col="gray75")
  
  return(mds.coor)
}

# test
distance_matrix = get_distance_matrix(520)
coor_by_distance_matrix = build_coor_by_distance(distance_matrix)

## 1.《hierarchical clustering》(hclust)
clust_res = hclust(dist(distance_matrix), method="single")
# 层次聚类 hierarchical clustering
plot(hclust(dist(distance_matrix), method="single"))




#########################################################################
### plot color of clusters 
#https://stackoverflow.com/questions/36067510/specify-color-for-cluster-elements-in-r
plot_cluster_with_color <- function(input_data, clus_cluster) {
  library("cluster")
  #solve of the bug:'from' must be a finite number."
  input_data = 10-input_data
  clusplot(input_data, clus_cluster,color=TRUE,shade=TRUE,col.p = clus_cluster, lines = 0)
  # clusplot(data1, kmres$cluster,color=TRUE,shade=TRUE,col.p = kmres$cluster)
}

#########################################################################
### dbscan cluster function
cluster_by_dbscan <- function(input_data) {
  #install.packages("dbscan")
  library("dbscan")
  clus <- dbscan(input_data, eps = .18, minPts = 3)
  clus_cluster = clus$cluster
  #names(clus_cluster) = c(0:651)
  names(clus_cluster) = c(rownames(input_data))

  table(clus_cluster)
  return(clus_cluster)
}


#test
distance_matrix = get_distance_matrix(520)
coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
clus_cluster = cluster_by_dbscan(coor_by_distance_matrix) 
#plot_cluster_with_color(coor_by_distance_matrix, clus_cluster)
plot(coor_by_distance_matrix, col = clus_cluster+1L)
table(clus_cluster)

#clus <- kmeans(coor_by_distance_matrix, centers = 10)



#########################################################################
### optics cluster function
cluster_by_optics <- function(input_data) {
  library("dbscan")
  ### run OPTICS
  res <- optics(input_data, eps = 10,  minPts = 10)
  
  # optics_cut = extractDBSCAN
  res <- extractDBSCAN(res, eps_cl = .15)
  
  plot(res)
  plot(input_data, col = res$cluster+1L)
  
  clus_cluster = res$cluster
  table(clus_cluster)
  return(clus_cluster)  
} 
#test
distance_matrix = get_distance_matrix(520)
coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
clus_cluster = cluster_by_optics(coor_by_distance_matrix) 

