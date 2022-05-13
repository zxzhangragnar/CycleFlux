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

#########################################################################
struct_other_main <- function(output_path, res_path) {

  load(file.path(output_path, "cycle_directed.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/web_struct_cycid.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/hub_struct_cycid.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/idv_struct_cycid.RData"))

  nrl_struct_cycid = cycle_directed$cycid
  nrl_struct_cycid = setdiff(nrl_struct_cycid, web_struct_cycid)
  nrl_struct_cycid = setdiff(nrl_struct_cycid, hub_struct_cycid)
  nrl_struct_cycid = setdiff(nrl_struct_cycid, idv_struct_cycid)

  #存放4类结构的环序号的list
  struct_cycle_sort_list = list()

  struct_cycle_sort_list[["net"]] = web_struct_cycid
  struct_cycle_sort_list[["hub"]] = hub_struct_cycid
  struct_cycle_sort_list[["idv"]] = idv_struct_cycid
  struct_cycle_sort_list[["nrl"]] = nrl_struct_cycid


  #save
  res_sub_path = "2_cycle_struct_sort/result_struct"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "struct_cycle_sort_list.RData")

  #save(struct_cycle_sort_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/struct_cycle_sort_list.RData")
  save(struct_cycle_sort_list, file=res_file_path)

}






# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# struct_other_main(output_path, res_path)













##############
## plot cycle
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/cycle_directed.RData")
#
# ## cyc_all_distance_df
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
#
# cycle_num = length(cycle_directed[,1])
# distance_matrix = get_distance_matrix(cycle_num, cyc_all_distance_df)
# coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
# struct_vector = c(0:(cycle_num-1))
# names(struct_vector) = c(0:(cycle_num-1))
# for (i in 1:length(struct_vector)) {
#   if(names(struct_vector[i]) %in% web_struct_cycid) {
#     struct_vector[[i]] = 1
#   }else if(names(struct_vector[i]) %in% hub_struct_cycid) {
#     struct_vector[[i]] = 2
#   }else if(names(struct_vector[i]) %in% idv_struct_cycid) {
#     struct_vector[[i]] = 3
#   }else {
#     struct_vector[[i]] = 4
#   }
# }
#
# plot(coor_by_distance_matrix, col = struct_vector)




