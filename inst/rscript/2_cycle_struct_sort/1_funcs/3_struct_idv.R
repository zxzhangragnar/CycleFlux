
#########################################################################
# get coor_by_distance_matrix
get_distance_matrix <- function(cycle_num, cyc_all_distance_df) {
  #cycle_num = 248

  distance_scale <- cyc_all_distance_df$distance
  #距离矩阵
  #distance_matrix =matrix(distance_scale, nr=248, byrow=TRUE, dimnames=list(c(0:651), c(0:651)))
  distance_matrix =matrix(distance_scale, nr=cycle_num, byrow=TRUE, dimnames=list(c(1:cycle_num), c(1:cycle_num)))
  
  return(distance_matrix)
}




#########################################################################
### build_coor_by_distance
#https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
build_coor_by_distance <- function(distance_matrix) {
  tmp = distance_matrix
  d <- as.dist(tmp)
  mds.coor <- cmdscale(d, k = length(d-1))
  # plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
  # text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
  #      rownames(mds.coor), cex=0.8)
  # abline(h=0,v=0,col="gray75")
  
  return(mds.coor)
}






#########################################################################
### dbscan cluster function
cluster_by_dbscan_for_idv <- function(input_data, dbscan_eps, dbscan_mpts) {
  #install.packages("dbscan")
  library("dbscan")
  #clus <- dbscan(input_data, eps = .16, minPts = 1)
  clus <- dbscan(input_data, eps = dbscan_eps, minPts = dbscan_mpts)
  clus_cluster = clus$cluster
  #names(clus_cluster) = c(0:651)
  names(clus_cluster) = c(rownames(input_data))
  
  table(clus_cluster)
  return(clus_cluster)
}

##################################################################################
# 函数位于 2_distance_cluster_function.R
##################################################################################
# 筛选出可疑的idv结构的聚簇   (fre<3)
filter_suspicious_idv_struct <- function(clus_cluster, filter_cluster_size, if_plot) {
  #new sort class/// filter small cluster
  sort_table = table(clus_cluster)
  sort_names = names(sort_table)
  index_tobe_remove = c()
  for (i in 1:length(clus_cluster)) {
    print(i)
    sort_class = as.character(clus_cluster[[i]]) #0:50
    srttab_index = which(sort_names == sort_class)
    #filter_cluster_size = 20,     remove which < 20
    if(sort_table[[srttab_index]] > filter_cluster_size) {
      clus_cluster[[i]] = 800    #other classes all = 800
      index_tobe_remove = append(index_tobe_remove, names(clus_cluster[i]))
    }
  }
  
  if(!if_plot) {
    clus_cluster = clus_cluster[-which(names(clus_cluster) %in% index_tobe_remove)]
  }
  return(clus_cluster)
}


##################################################################################
### tool func
##################################################################################
# 【去掉hub环】
# 【去掉web环】
remove_other_struct <- function(suspicious_idv_cycle_ids, web_struct_cycle_id, hub_struct_cycle_id) {
  suspicious_idv_cycle_ids = setdiff(suspicious_idv_cycle_ids, web_struct_cycle_id)
  suspicious_idv_cycle_ids = setdiff(suspicious_idv_cycle_ids, hub_struct_cycle_id)
  idv_struct_cycle_id = suspicious_idv_cycle_ids
  return(idv_struct_cycle_id)
}


##################################################################################
# 检验
test_by_cycle_id_neighcyc_shared_cycle <- function(idv_struct_cycle_id, cycle_directed) {
  cycle_id_neighcyc_shared_cycle = cycle_directed[,c("cycle_id", "neighcyc", "shared_cycle")]
  cycle_id_neighcyc_shared_cycle = cycle_id_neighcyc_shared_cycle[which(cycle_id_neighcyc_shared_cycle$cycle_id %in% idv_struct_cycle_id),]
}
##################################################################################




##################################################################################
#test
struct_idv_main <- function(output_path, res_path, prm_1=.16, prm_2=1, prm_3=2, prm_4=20, prm_5=.56, prm_6=20) {
  
  load(file.path(output_path, "cycle_directed.RData"))
  load(file.path(res_path, "1_cycle_topology/result_topo/cyc_all_distance_df.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/hub_struct_cycle_id.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/web_struct_cycle_id.RData"))
  
  cycle_num = length(cycle_directed[,1])
  distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
  coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
  clus_cluster = cluster_by_dbscan_for_idv(coor_by_distance_matrix, prm_1, prm_2) 
  # plot(coor_by_distance_matrix, col = clus_cluster+1L)
  # table(clus_cluster)
  
  #suspicious_idv_struct
  suspicious_idv_struct_cluster = filter_suspicious_idv_struct(clus_cluster, prm_3, FALSE)
  # table(suspicious_idv_struct_cluster)
  # plot(coor_by_distance_matrix, col = suspicious_idv_struct_cluster+1L)
  #plot_cluster_with_color(coor_by_distance_matrix, suspicious_idv_struct_cluster)
  
  suspicious_idv_cycle_ids = names(suspicious_idv_struct_cluster)
  #remove web hub
  suspicious_idv_cycle_ids = remove_other_struct(suspicious_idv_cycle_ids, web_struct_cycle_id, hub_struct_cycle_id)
  # 检验
  cluster_idv_ngh_df = test_by_cycle_id_neighcyc_shared_cycle(suspicious_idv_cycle_ids, cycle_directed)
  cluster_idv_ngh_df = cluster_idv_ngh_df[which(cluster_idv_ngh_df$shared_cycle < prm_4),]
  avg_shc1 = mean(cluster_idv_ngh_df$shared_cycle)
  
  idv_struct_cycle_id_1 = cluster_idv_ngh_df$cycle_id
  
  
  
  ##################################################################################
  # 在R语言中的“HDoutliers”包：Leland Wilkinson 检测多维异常值的算法 (多维数据)
  #test
  cycle_num = length(cycle_directed[,1])
  distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
  coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
  # plot(coor_by_distance_matrix)
  
  #suspicious_idv_struct
  library(HDoutliers)
  outlier_coor_by_distance_matrix <- HDoutliers(coor_by_distance_matrix, alpha=prm_5, transform = TRUE)
  # plotHDoutliers(coor_by_distance_matrix, outlier_coor_by_distance_matrix)
  outlier_idv_cycle_ids = outlier_coor_by_distance_matrix
  
  #remove web hub
  outlier_idv_cycle_ids = remove_other_struct(outlier_idv_cycle_ids, web_struct_cycle_id, hub_struct_cycle_id)
  # 检验
  outlier_idv_ngh_df = test_by_cycle_id_neighcyc_shared_cycle(outlier_idv_cycle_ids, cycle_directed)
  outlier_idv_ngh_df = outlier_idv_ngh_df[which(outlier_idv_ngh_df$shared_cycle < prm_6),]
  avg_shc2 = mean(outlier_idv_ngh_df$shared_cycle)
  
  idv_struct_cycle_id_2 = outlier_idv_ngh_df$cycle_id
  
  
  ##################################################################################
  ### final result (merged)
  idv_struct_cycle_id = append(idv_struct_cycle_id_1, idv_struct_cycle_id_2)
  idv_struct_cycle_id = unique(idv_struct_cycle_id)
  
  # 检验
  idv_merge_compare_df = cycle_directed[which(cycle_directed$cycle_id %in% idv_struct_cycle_id), c("cycle_id","neighcyc","shared_cycle")]
  
  
  #save
  res_sub_path = "2_cycle_struct_sort/result_struct"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "idv_struct_cycle_id.RData")
  
  #save(idv_struct_cycle_id, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/idv_struct_cycle_id.RData")
  save(idv_struct_cycle_id, file=res_file_path)
  
}











#######################################################################################################################                 
# # test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# 
# struct_idv_main(output_path, res_path, .16, 1, 2, 20, .56, 20)







