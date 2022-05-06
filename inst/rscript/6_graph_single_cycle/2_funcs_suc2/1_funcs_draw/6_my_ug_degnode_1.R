
#情况1:  e2<--up--e1<--up--c2<--up--c1
####################################################################





####################################################################
suc2_my_ug_degnode_1_main <- function(output_path, res_path, input_tumor_name) {
  
  # init
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_data/cycle_edge_flux_list_2_1.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_data/cycle_edge_flux_list_2_2.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(output_path, "res_cycle_successors_2.RData"))
  
  ##
  tumor_name = input_tumor_name
  st2_1_gapinfo = cycle_edge_flux_list_2_1[[tumor_name]]
  st2_2_gapinfo = cycle_edge_flux_list_2_2[[tumor_name]]
  st2_1_gapinfo = cbind(cycle_degnode_situation_2$cycid, st2_1_gapinfo)
  st2_2_gapinfo = cbind(cycle_degnode_situation_2$cycid, st2_2_gapinfo)
  colnames(st2_1_gapinfo)[1] = "cycid"
  colnames(st2_2_gapinfo)[1] = "cycid"
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  up_cout_node = ug_cyc_pct_upnode$c_out
  
  #  e1<----c2<--up--c1
  select_st2_1_gapinfo = data.frame()
  select_st2_2_gapinfo = data.frame()
  
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    
    tmp_row1 = st2_1_gapinfo[which((st2_1_gapinfo$cycid == tmp_cycid) & (st2_1_gapinfo$node == tmp_cout)),]
    tmp_row2 = st2_2_gapinfo[which((st2_1_gapinfo$cycid == tmp_cycid) & (st2_1_gapinfo$node == tmp_cout)),]
    
    select_st2_1_gapinfo = rbind(select_st2_1_gapinfo, tmp_row1)
    select_st2_2_gapinfo = rbind(select_st2_2_gapinfo, tmp_row2)
    
  }
  
  
  #  e1<--up--c2<--up--c1
  select_st2_1_gapinfo_e1 = select_st2_1_gapinfo[which(select_st2_1_gapinfo$ifup == "up"),]
  #取相同的索引
  select_st2_2_gapinfo_e1 = select_st2_2_gapinfo[which(select_st2_1_gapinfo$ifup == "up"),]
  
  #  e2<--up--e1<--up--c2<--up--c1
  select_st2_2_gapinfo_e1e2 = select_st2_2_gapinfo_e1[which(select_st2_2_gapinfo_e1$ifup == "up"),]
  #取相同的索引
  select_st2_1_gapinfo_e1e2 = select_st2_1_gapinfo_e1[which(select_st2_2_gapinfo_e1$ifup == "up"),]
  
  
  # save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc2/cyc_successors_result"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path_st2_1_gapinfo_e1 = file.path(res_path, res_sub_path, "select_st2_1_gapinfo_e1.RData")
  res_file_path_st2_2_gapinfo_e1 = file.path(res_path, res_sub_path, "select_st2_2_gapinfo_e1.RData")
  res_file_path_st2_1_gapinfo_e1e2 = file.path(res_path, res_sub_path, "select_st2_1_gapinfo_e1e2.RData")
  res_file_path_st2_2_gapinfo_e1e2 = file.path(res_path, res_sub_path, "select_st2_2_gapinfo_e1e2.RData")
  
  save(select_st2_1_gapinfo_e1, file=res_file_path_st2_1_gapinfo_e1)
  save(select_st2_2_gapinfo_e1, file=res_file_path_st2_2_gapinfo_e1)
  save(select_st2_1_gapinfo_e1e2, file=res_file_path_st2_1_gapinfo_e1e2)
  save(select_st2_2_gapinfo_e1e2, file=res_file_path_st2_2_gapinfo_e1e2)
  
}





# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# 
# suc2_my_ug_degnode_1_main(output_path, res_path, input_tumor_name)




































