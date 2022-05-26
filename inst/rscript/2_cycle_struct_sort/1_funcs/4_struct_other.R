# 
# web       1
# hub       2
# idv       3  
# nrl       4
#

#########################################################################
# get coor_by_distance_matrix
get_distance_matrix <- function(cycle_num, cyc_all_distance_df) {
  distance_scale <- cyc_all_distance_df$distance
  #距离矩阵
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
struct_other_main <- function(output_path, res_path) {
  
  load(file.path(output_path, "cycle_directed.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/web_struct_cycle_id.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/hub_struct_cycle_id.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/idv_struct_cycle_id.RData"))
  
  nrl_struct_cycle_id = cycle_directed$cycle_id
  nrl_struct_cycle_id = setdiff(nrl_struct_cycle_id, web_struct_cycle_id)
  nrl_struct_cycle_id = setdiff(nrl_struct_cycle_id, hub_struct_cycle_id)
  nrl_struct_cycle_id = setdiff(nrl_struct_cycle_id, idv_struct_cycle_id)
  
  #存放4类结构的环序号的list
  struct_cycle_sort_list = list()
  
  struct_cycle_sort_list[["net"]] = web_struct_cycle_id
  struct_cycle_sort_list[["hub"]] = hub_struct_cycle_id
  struct_cycle_sort_list[["idv"]] = idv_struct_cycle_id
  struct_cycle_sort_list[["nrl"]] = nrl_struct_cycle_id
  
  
  #save
  res_sub_path = "2_cycle_struct_sort/result_struct"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "struct_cycle_sort_list.RData")
  
  save(struct_cycle_sort_list, file=res_file_path)
}

#########################################################################################################################                 
## test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# struct_other_main(output_path, res_path)


