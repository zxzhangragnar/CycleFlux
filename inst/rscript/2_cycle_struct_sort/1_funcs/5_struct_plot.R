
#此函数内容不必执行，与结果无关

#########################################################################
# get coor_by_distance_matrix
get_distance_matrix_compound <- function(cpd_num) {
  #compound_num = 2615
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/all_cpd_distance_frame.RData")
  all_cpd_distance_frame = all_cpd_distance_frame
  distance_scale <- all_cpd_distance_frame$distance
  #距离矩阵
  #distance_matrix =matrix(distance_scale, nr=652, byrow=TRUE, dimnames=list(c(0:651), c(0:651)))
  distance_matrix =matrix(distance_scale, nr=cpd_num, byrow=TRUE, dimnames=list(c(0:(cpd_num-1)), c(0:(cpd_num-1))))
  
  return(distance_matrix)
}


#########################################################################
### build_coor_by_distance
build_coor_by_distance <- function(distance_matrix) {
  tmp = distance_matrix
  d <- as.dist(tmp)
  for (i in 1:length(d)) {
    if(is.na(d[i])) {
      d[i] = 20 #inf
    }
  }
  mds.coor <- cmdscale(d)
  plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="")
  text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
       rownames(mds.coor), cex=0.8)
  abline(h=0,v=0,col="gray75")
  
  return(mds.coor)
}









############################## all_cpd_distance_frame ###############################

get_temp_frame <- function(i, cpdnames) {
  temp_frame = as.data.frame(cpdnames)
  temp_frame = cbind(temp_frame, c(NA))
  colnames(temp_frame) = c("to", "distance")
  
  part_df = cpd_all_distance_df[which(cpd_all_distance_df$from == cpdnames[i]), c("to","distance")]
  part_df = unique(part_df)
  temp_frame = rbind(part_df, temp_frame)
  
  temp_df = as.data.frame(cpdnames)
  temp_df[,"temp"] = 1 #any value not useful
  colnames(temp_df) = c("to", "temp")
  
  require(data.table)
  setDT(temp_frame)
  setDT(temp_df); 
  temp_frame = temp_frame[temp_df, mult = "first", on = "to", nomatch=0L]
  temp_frame = temp_frame[, c("to","distance")]
  temp_frame = cbind(cpdnames[i], temp_frame)
  colnames(temp_frame) = c("from", "to", "distance")
  
  return(temp_frame)
}


load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_compound_union_directed.RData")
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cpd_all_distance_df.RData")
cpd_all_distance_df = unique(cpd_all_distance_df)
cpdnames = compound_directed$cpdname
cpdnames = unique(cpdnames)  # 331
all_cpd_distance_frame = data.frame()
for (i in 1:length(cpdnames)) {
  temp_frame = get_temp_frame(i, cpdnames)
  all_cpd_distance_frame = rbind(all_cpd_distance_frame, temp_frame)
}
# save
save(all_cpd_distance_frame, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/all_cpd_distance_frame.RData")




##############
## plot compound
cpd_num = length(cpdnames)
distance_matrix = get_distance_matrix_compound(cpd_num)
coor_by_distance_matrix = build_coor_by_distance(distance_matrix)
plot(distance_matrix)
#plot(coor_by_distance_matrix, col = struct_vector)






















######################################## learn ###################################
# https://www.cnblogs.com/jessepeng/p/11072533.html
# data = data.frame(id = c(1,2,3,4,5),
#                   state = c("KS","MN","AL","FL","CA"))
# scores = data.frame(id = c(1,1,1,2,2,3,3,3),
#                     score = c(66,75,78,86,85,76,75,90))
# 
# require(data.table)
# setDT(scores);
# setDT(data) # convert to data.tables by reference
# scores[data, mult = "first", on = "id", nomatch=0L]








