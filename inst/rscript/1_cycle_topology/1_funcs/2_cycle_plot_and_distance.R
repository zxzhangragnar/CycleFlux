#   
############################### 在以下计算中，将"环"看成网络中的节点 ####################################
# 
# 
# 以"环"为点  思想 的延伸
# 
# 如何能plot出 以环作为点 的图像?
#   
# 求 环 和 环 之间的距离 distance 
#
#
# 环的距离数据 from python Functions_new/cycle_distance.py
#
#
#



##################################### part1 build graph ################################################
#############################################################################################

get_cyc_all_distance_df <- function(cycle_directed) {
  cycle_distance = as.data.frame(cycle_directed[,c("cycid", "distance")])
  cyc_all_distance_df = data.frame()
  for (i in 1:length(rownames(cycle_distance))) {
    print(paste0("cyc+",i))
    temp_cycid = cycle_directed[i, "cycid"]
    temp_distance_str = as.character(cycle_distance[i,"distance"])
    temp_distance_str = substring(temp_distance_str, 2,nchar(temp_distance_str)-1)
    temp_distance_str <- unlist(strsplit(temp_distance_str,split = ", "))
    for (j in 1:length(temp_distance_str)) {
      temp_distance_str_cyc_val = unlist(strsplit(temp_distance_str[j],split = ": "))
      temp_distance_cycid = temp_distance_str_cyc_val[1]
      temp_distance_val = temp_distance_str_cyc_val[2]
      temp_df = data.frame(
        from = c(temp_cycid), 
        to = c(temp_distance_cycid),
        distance = c(temp_distance_val)
      )
      cyc_all_distance_df = rbind(cyc_all_distance_df, temp_df)
    }
  }
  return(cyc_all_distance_df)
}


#############################################################################################
# test
cycle_plot_and_distance_main <- function(output_path, res_path) {
  load(file.path(output_path, "res_allpathway_cycle_union_directed.RData"))
  cyc_all_distance_df = get_cyc_all_distance_df(cycle_directed)
  
  # save
  res_sub_path = "1_cycle_topology/result_topo"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cyc_all_distance_df.RData")
  
  #save(cyc_all_distance_df, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
  save(cyc_all_distance_df, file=res_file_path)
  
}



# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# cycle_plot_and_distance_main(output_path, res_path) 





