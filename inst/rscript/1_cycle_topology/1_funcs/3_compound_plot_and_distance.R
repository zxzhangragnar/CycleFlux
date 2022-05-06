#此函数内容不必执行，与结果无关

############################### 在以下计算中，将"compound"看成网络中的节点 ####################################

##################################### part1 build graph ################################################
#############################################################################################

#build graph


# cycle_directed
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_compound_union_directed.RData")
compound_distance = as.data.frame(compound_directed[,c("cpdname", "distance")])
cpd_all_distance_df = data.frame()
for (i in 1:length(rownames(compound_distance))) {
  print(paste0("cpd+",i))
  temp_cdpname = compound_directed[i, "cpdname"]
  temp_distance_str = as.character(compound_distance[i,"distance"])
  temp_distance_str = substring(temp_distance_str, 2,nchar(temp_distance_str)-1)
  temp_distance_str <- unlist(strsplit(temp_distance_str,split = ", "))
  for (j in 1:length(temp_distance_str)) {
    temp_distance_str_cdp_val = unlist(strsplit(temp_distance_str[j],split = ": "))
    temp_distance_cdpname = substring(temp_distance_str_cdp_val[1], 2,nchar(temp_distance_str_cdp_val[1])-1)
    temp_distance_val = temp_distance_str_cdp_val[2]
    temp_df = data.frame(
      from = c(temp_cdpname), 
      to = c(temp_distance_cdpname),
      distance = c(temp_distance_val)
    )
    cpd_all_distance_df = rbind(cpd_all_distance_df, temp_df)
  }
}


# save 
# save(cpd_all_distance_df, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cpd_all_distance_df.RData")







