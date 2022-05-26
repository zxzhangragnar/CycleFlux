
#1.21 joint meeting 

#cycle_edge_flux_list
#load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")


# gap_c<-unique(cycle_edge_flux_list[[1]][which((ddd<(-2))|(ddd2<1.71)),1])
# up_c<-unique(cycle_edge_flux_list[[1]][which((ddd>(1))&(ddd2>10)),1])#p.value of diff<0.001, log(fc)>0.5, and the gene show significant up regulation are larger than 10
# gap_up_c = intersect(up_c,gap_c)
####################

########################################################
##########
#analysis
# R中主要是两个rdata文件:
# 1. 每一列为环中的1个反应  (cycle_edge_flux_list)
# 2.每一列为1个环 (cycle_flux_list)
# 
# (始终保持这两个文件为最新，并用最新的这两个文件做分析 即可)

###############################################
#1.看有gap的环,除了gap以外的rct是否为up 
#2.重新做boxplot等图分析(用cycle_flux_list)
#3.使用genecards.org更新缺失的基因value
###############################################
## 此代码生成3个结果: (供其它代码使用)
# 1.never_considered_compounds
# 2.gapup_cycle_chain_list
# 均保存在 result_analysis/ 目录下

###############################################################################
########################## 1.never_considered_compounds #############################

##1.在统计化合物时, 删去never_considered_comp.csv中的化合物
get_never_considered_comp<-function(never_considered_comp, compounds_dict) {
  never_considered_comp_arr = never_considered_comp$x
  compounds_dict_cids = compounds_dict$ENTRY
  never_considered_comp_arr_cids = intersect(never_considered_comp_arr, compounds_dict_cids)
  never_considered_comp_names = never_considered_comp_arr_cids
  
  return(never_considered_comp_names)
}






###############################################################################
################################ 2.gapup_cycle_chain_list #######################
get_ug_chain_list<-function(tumors_array, cycle_edge_flux_list, all_chain_list_cid) {
  ug_chain_list_dict = list()
  for (j in 1:length(tumors_array)) {
    tumor_name = tumors_array[j]
    temp_tumor_df = cycle_edge_flux_list[[tumor_name]]
    #有gap现象的环的cycle_id
    gap_c = unique(temp_tumor_df[which(temp_tumor_df$ifgap == "gap"),'cycle_id'])
    up_c = unique(temp_tumor_df[which(temp_tumor_df$ifup == "up"),'cycle_id'])
    ug_c = intersect(gap_c, up_c)
    
    temp_chain_list = all_chain_list_cid[[tumor_name]]
    cycle_len = length(temp_chain_list)
    #names(temp_chain_list) = c(0:247)
    names(temp_chain_list) = c(1:cycle_len)
    ug_c = as.integer(ug_c)
    ug_chain_list = temp_chain_list[ug_c]
    
    ug_chain_list_dict[[tumor_name]] = ug_chain_list
  }
  return(ug_chain_list_dict)
}



##############################################################################
single_cycle_gap_analysis_main <- function(output_path, res_path, package_path, input_tumor_name) {
  
  # init
  library(readr)
  never_considered_comp = read_csv(file.path(package_path, "tool_data/never_considered_comp.csv"))
  
  load(file.path(res_path, "4_flux_edge/result_final/compounds_dict.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/TCGA_gap_all_chain_list_cid.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))

  ##
  never_considered_comp_names = get_never_considered_comp(never_considered_comp, compounds_dict)
  tumors_array = c(input_tumor_name)
  gapup_cycle_chain_list = get_ug_chain_list(tumors_array, cycle_edge_flux_list, all_chain_list_cid)
  
  
  #save
  res_sub_path = "6_graph_single_cycle/result_analysis"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path_nc = file.path(res_path, res_sub_path, "never_considered_compounds.RData")
  res_file_path_chain = file.path(res_path, res_sub_path, "gapup_cycle_chain_list.RData")
  
  # save(never_considered_comp_arr,never_considered_comp_names, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/never_considered_compounds.RData")
  # save(gapup_cycle_chain_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData")
  save(never_considered_comp_names, file=res_file_path_nc)
  save(gapup_cycle_chain_list, file=res_file_path_chain)
  
}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
# 
# single_cycle_gap_analysis_main(output_path, res_path, package_path, input_tumor_name)

