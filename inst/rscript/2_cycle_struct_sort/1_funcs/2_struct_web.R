# 判断标准：
# hub 在去掉1个(泛化成2个3个) 枢纽节点后，很多环就不相连了(泛化成连通性大幅度减少了)
# web 去掉1个枢纽节点后 仍然成环，有很多overlap相连的点
#
#
# 先找到 Hub  ---> 去掉(不考虑)这些环 ---> 再找 web
#
#
# （根据距离distance分布进行聚类）
#
# 【可疑聚簇】
# 根据距离聚类筛选出来的 结果clus_cluster
# 筛选出可疑的web结构的聚簇   （fre>=20）
# 筛选出可疑的idv结构的聚簇   (fre<=2)
#
#
# 【判断标准】
# 根据 distance 等指标 找到可疑 Web 结构聚簇
# 【标准】对网内每一个环判断标准：
# （网内每一个环的neighcyc都大于20）
# （网内所有环之间的距离distance都小于2）
#
# 【去掉hub环】
# 去掉之前找到的属于hub结构的环
#



#########################################################################
# get coor_by_distance_matrix
get_distance_matrix <- function(cycle_num, cyc_all_distance_df) {
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
  # plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
  # text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
  #      rownames(mds.coor), cex=0.8)
  # abline(h=0,v=0,col="gray75")

  return(mds.coor)
}



########################################################################
# 得到 Hub结构及其所属的环


#########################################################################
### dbscan cluster function
cluster_by_dbscan_for_web <- function(input_data, p_eps, p_minpts) {
  #install.packages("dbscan")
  library("dbscan")
  #clus <- dbscan(input_data, eps = .18, minPts = 3)
  clus <- dbscan(input_data, eps = p_eps, minPts = p_minpts)
  clus_cluster = clus$cluster
  #names(clus_cluster) = c(0:651)
  names(clus_cluster) = c(rownames(input_data))

  table(clus_cluster)
  return(clus_cluster)
}

##################################################################################
##################################################################################
#【可疑聚簇】
# 函数位于 2_distance_cluster_function.R
# 筛选出可疑的web结构的聚簇   （fre>=20）
filter_suspicious_web_struct <- function(clus_cluster, filter_cluster_size, if_plot) {
  #new sort class/// filter small cluster
  sort_table = table(clus_cluster)
  sort_names = names(sort_table)
  index_tobe_remove = c()
  for (i in 1:length(clus_cluster)) {
    print(i)
    sort_class = as.character(clus_cluster[[i]]) #0:50
    srttab_index = which(sort_names == sort_class)
    #filter_cluster_size = 20,     remove which < 20
    if(sort_table[[srttab_index]] < filter_cluster_size) {
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
##################################################################################
# 【判断标准】
# 【标准】对网内每一个环判断标准：
get_suspicious_web_table_list <- function(suspicious_web_struct_cluster) {
  suspicious_web_table_list = list()
  if (length(suspicious_web_struct_cluster) > 0) {
    suspicious_web_table_list = names(table(suspicious_web_struct_cluster))
    for (i in 1:length(suspicious_web_table_list)) {
      temp_clus_name = suspicious_web_table_list[[i]]
      temp_df = data.frame()
      temp_df_index = 1
      for (j in 1:length(suspicious_web_struct_cluster)) {
        if(as.integer(suspicious_web_struct_cluster[j]) == temp_clus_name) {
          cycid = names(suspicious_web_struct_cluster[j])
          temp_df[temp_df_index, "cycid"] = cycid
          temp_df_index = temp_df_index + 1
        }
      }
      suspicious_web_table_list[[i]] = temp_df
    }
  }

  return(suspicious_web_table_list)
}



##################################################################################
# （网内所有环之间的距离distance都小于4）
get_far_cycid_pairs_df_list <- function(suspicious_web_table_list, threshold, cycle_directed, cyc_all_distance_df) {
  cycle_num = length(cycle_directed[,1])
  distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
  far_cycid_pairs_df_list = list()
  if (length(suspicious_web_table_list) > 0) {
    for (k in 1:length(suspicious_web_table_list)) {
      suspicious_web_table = suspicious_web_table_list[[k]]
      suspicious_web_cycids = suspicious_web_table$cycid

      far_cycid_pairs_df = data.frame()
      df_row_index = 1
      for (i in 1:length(suspicious_web_cycids)) {
        print(i)
        cyc_a = as.character(suspicious_web_cycids[i])
        for (j in 1:length(suspicious_web_cycids)) {
          cyc_b = as.character(suspicious_web_cycids[j])
          temp_distance = distance_matrix[cyc_a, cyc_b]
          if(temp_distance > threshold) {
            far_cycid_pairs_df[df_row_index, "from"] = cyc_a
            far_cycid_pairs_df[df_row_index, "to"] = cyc_b
            far_cycid_pairs_df[df_row_index, "dist"] = temp_distance
            df_row_index = df_row_index + 1
          }
        }
      }

      far_cycid_pairs_df_list[[k]] = far_cycid_pairs_df
    }
  }

  return(far_cycid_pairs_df_list)
}

#get sum dist
get_sumdist_df_list <- function(suspicious_web_table_list, cycle_directed, cyc_all_distance_df) {
  cycle_num = length(cycle_directed[,1])
  distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
  sumdist_df_list = list()
  if (length(suspicious_web_table_list) > 0) {
    for (k in 1:length(suspicious_web_table_list)) {
      suspicious_web_table = suspicious_web_table_list[[k]]
      suspicious_web_cycids = suspicious_web_table$cycid

      sumdist_df = as.data.frame(suspicious_web_cycids)
      sumdist_df[,"sumdist"] = 0
      for (i in 1:length(suspicious_web_cycids)) {
        print(i)
        cyc_a = as.character(suspicious_web_cycids[i])
        for (j in 1:length(suspicious_web_cycids)) {
          cyc_b = as.character(suspicious_web_cycids[j])
          temp_distance = as.integer(distance_matrix[cyc_a, cyc_b])
          sumdist_df[which(sumdist_df$suspicious_web_cycids == cyc_a),"sumdist"] = sumdist_df[which(sumdist_df$suspicious_web_cycids == cyc_a),"sumdist"] + temp_distance

        }
      }

      sumdist_df_list[[k]] = sumdist_df
    }
  }

  return(sumdist_df_list)
}

#离群点
get_outlier_list <- function(far_cycid_pairs_df_list, sumdist_df_list) {
  outlier_list = c()
  if ((length(far_cycid_pairs_df_list) > 0) & (length(sumdist_df_list) > 0)) {
    # for every pair, remove the dist bigger one
    for (i in 1:length(far_cycid_pairs_df_list)) {
      temp_far_cycid_pairs_df = far_cycid_pairs_df_list[[i]]
      if (length(far_cycid_pairs_df_list[[i]]) == 0) {
        far_cycid_pairs_df_list[[i]] = temp_far_cycid_pairs_df
      }else {
        sumdist_df = sumdist_df_list[[i]]
        for(j in 1:length(rownames(temp_far_cycid_pairs_df))) {
          from_cycid = temp_far_cycid_pairs_df[j, "from"]
          to_cycid = temp_far_cycid_pairs_df[j, "to"]
          from_cycid_dist = sumdist_df[which(sumdist_df$suspicious_web_cycids == from_cycid), "sumdist"]
          to_cycid_dist = sumdist_df[which(sumdist_df$suspicious_web_cycids == to_cycid), "sumdist"]
          if(from_cycid_dist >= to_cycid_dist) {
            outlier_list = append(outlier_list, from_cycid)
          }else {
            outlier_list = append(outlier_list, to_cycid)
          }
        }
      }
    }
  }
  #remove离群点
  outlier_list = unique(outlier_list)
  return(outlier_list)

}




##################################################################################
# （网内每一个环的neighcyc都大于20）
filter_neighcyc_bigger <- function(suspicious_web_cycids, threshold, cycle_directed) {
  cycid_neighcyc_sharedcyc = cycle_directed[,c("cycid", "neighcyc", "sharedcyc")]
  cycid_neighcyc_sharedcyc = cycid_neighcyc_sharedcyc[which(cycid_neighcyc_sharedcyc$cycid %in% suspicious_web_cycids),]
  cycid_neighcyc_sharedcyc = cycid_neighcyc_sharedcyc[which(cycid_neighcyc_sharedcyc$neighcyc > threshold),]
  cycid_neighcyc_sharedcyc = cycid_neighcyc_sharedcyc[which(cycid_neighcyc_sharedcyc$sharedcyc > threshold),]
  suspicious_web_cycids = cycid_neighcyc_sharedcyc[,c("cycid")]
  return(suspicious_web_cycids)
}



##################################################################################






#############################################################################################
# test
struct_web_main <- function(output_path, res_path, par1=.18, par2=3, par3=7, par4=2, par5=10) {

  load(file.path(output_path, "cycle_directed.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/hub_struct_cycid.RData"))
  load(file.path(res_path, "1_cycle_topology/result_topo/cyc_all_distance_df.RData"))

  cycle_num = length(cycle_directed[,1])
  distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
  coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
  clus_cluster = cluster_by_dbscan_for_web(coor_by_distance_matrix, par1, par2)
  #plot_cluster_with_color(coor_by_distance_matrix, clus_cluster)
  # plot(coor_by_distance_matrix, col = clus_cluster+1L)
  # table(clus_cluster)

  #suspicious_web_struct
  #suspicious_web_struct_cluster = filter_suspicious_web_struct(clus_cluster, 15, FALSE)
  suspicious_web_struct_cluster = filter_suspicious_web_struct(clus_cluster, par3, FALSE)

  table(suspicious_web_struct_cluster)
  # plot(coor_by_distance_matrix, col = suspicious_web_struct_cluster+1L)
  #plot_cluster_with_color(coor_by_distance_matrix, suspicious_web_struct_cluster)


  #test
  suspicious_web_table_list = get_suspicious_web_table_list(suspicious_web_struct_cluster)
  far_cycid_pairs_df_list = get_far_cycid_pairs_df_list(suspicious_web_table_list, par4, cycle_directed, cyc_all_distance_df)


  #test
  suspicious_web_table_list = get_suspicious_web_table_list(suspicious_web_struct_cluster)
  sumdist_df_list = get_sumdist_df_list(suspicious_web_table_list, cycle_directed, cyc_all_distance_df)


  outlier_list = get_outlier_list(far_cycid_pairs_df_list, sumdist_df_list)
  ###
  suspicious_web_cycids = names(suspicious_web_struct_cluster)
  ### remove 离群点
  suspicious_web_cycids = setdiff(suspicious_web_cycids, outlier_list)


  #test
  #suspicious_web_cycids = filter_neighcyc_bigger(suspicious_web_cycids, 25)
  suspicious_web_cycids = filter_neighcyc_bigger(suspicious_web_cycids, par5, cycle_directed)


  # 【去掉hub环】
  # 得到 Hub结构及其所属的环
  ### remove Hub struct
  suspicious_web_cycids = setdiff(suspicious_web_cycids, hub_struct_cycid)

  #Final result:
  web_struct_cycid = suspicious_web_cycids


  #save
  res_sub_path = "2_cycle_struct_sort/result_struct"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "web_struct_cycid.RData")

  #save(web_struct_cycid, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/web_struct_cycid.RData")
  save(web_struct_cycid, file=res_file_path)

}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# struct_web_main(output_path, res_path, .18, 3, 7, 2, 10)






##################################################################################
### plot the result
# plot_the_web_cycids <- function(clus_cluster, suspicious_web_cycids) {
#   plot_web_cycids = clus_cluster
#   for (i in 1:length(plot_web_cycids)) {
#     if(names(plot_web_cycids[i]) %in% suspicious_web_cycids) {
#       plot_web_cycids[[i]] = 80
#     }else {
#       plot_web_cycids[[i]] = 90
#     }
#   }
#   plot(coor_by_distance_matrix, col = plot_web_cycids+1L)
# }
#
# #plot
# plot_the_web_cycids(clus_cluster, suspicious_web_cycids)










