
#情况1:  e2<--up--c2<--up--c1
#主要负责出度
####################################################################

get_select_out_upinfo <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"
  
  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  
  #  e1<----c2<--up--c1
  select_out_upinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    
    tmp_row1 = out_gapinfo[which((out_gapinfo$cycid == tmp_cycid) & (out_gapinfo$node == tmp_cout)),]
    select_out_upinfo = rbind(select_out_upinfo, tmp_row1)
  }
  
  #  e1<--up--c2<--up--c1
  select_out_upinfo = select_out_upinfo[which(select_out_upinfo$ifup == "up"),]
  return(select_out_upinfo)
}


get_select_out_upinfo_both <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"
  
  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  
  #  e1<----c2<--up--c1
  select_out_upinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]
    
    tmp_row1 = out_gapinfo[which((out_gapinfo$cycid == tmp_cycid) & ((out_gapinfo$node == tmp_cout) | (out_gapinfo$node == tmp_cin))),]
    select_out_upinfo = rbind(select_out_upinfo, tmp_row1)
  }
  
  #  e1<--up--c2<--up--c1
  select_out_upinfo = select_out_upinfo[which(select_out_upinfo$ifup == "up"),]
  return(select_out_upinfo)
}



get_select_out_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"
  
  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  
  #不仅取up那条边周围的 degree, 而是取所有的 不管是环中的哪个反应
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct
  
  #  e1<----c2<--up--c1
  select_out_upinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    
    tmp_row1 = out_gapinfo[which((out_gapinfo$cycid == tmp_cycid) & (out_gapinfo$node == tmp_cout)),]
    select_out_upinfo = rbind(select_out_upinfo, tmp_row1)
  }
  
  #  e1<--up--c2<--up--c1
  select_out_upinfo = select_out_upinfo[which(select_out_upinfo$ifup == "up"),]
  return(select_out_upinfo)
}

##################################################################################
### 以下为 gap 的内容


get_select_out_gapinfo <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"
  
  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifgap == "gap"),]
  #up_cout_node %in% st2_1_gapinfo$node
  
  #  e1<----c2<--up--c1
  select_out_gapinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    
    tmp_row1 = out_gapinfo[which((out_gapinfo$cycid == tmp_cycid) & (out_gapinfo$node == tmp_cout)),]
    select_out_gapinfo = rbind(select_out_gapinfo, tmp_row1)
  }
  
  #  e1<--up--c2<--up--c1
  select_out_gapinfo = select_out_gapinfo[which(select_out_gapinfo$ifgap == "gap"),]
  return(select_out_gapinfo)
}



get_select_out_gapinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"
  
  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]
  
  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifgap == "gap"),]
  ug_cyc_pct_upnode = ug_cyc_pct
  
  #  e1<----c2<--up--c1
  select_out_gapinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]
    
    tmp_row1 = out_gapinfo[which((out_gapinfo$cycid == tmp_cycid) & (out_gapinfo$node == tmp_cout)),]
    select_out_gapinfo = rbind(select_out_gapinfo, tmp_row1)
  }
  
  #  e1<--up--c2<--up--c1
  select_out_gapinfo = select_out_gapinfo[which(select_out_gapinfo$ifgap == "gap"),]
  return(select_out_gapinfo)
}




#####################################################################################
my_ug_degnode_1_main <- function(output_path, res_path, input_tumor_name) {
  
  # init
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_out.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))

  ##
  tumors_array = c(input_tumor_name)
  select_out_upinfo = list()
  select_out_upinfo_both = list()
  select_out_upinfo_eachedge = list()
  select_out_gapinfo = list()
  select_out_gapinfo_eachedge = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    
    ##up
    temp_select_out_upinfo = get_select_out_upinfo(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_out_upinfo_both = get_select_out_upinfo_both(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_out_upinfo_eachedge = get_select_out_upinfo_eachedge(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    ##gap
    temp_select_out_gapinfo = get_select_out_gapinfo(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_out_gapinfo_eachedge = get_select_out_gapinfo_eachedge(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    
    select_out_upinfo[[tumor_name]] = temp_select_out_upinfo
    select_out_upinfo_both[[tumor_name]] = temp_select_out_upinfo_both
    select_out_upinfo_eachedge[[tumor_name]] = temp_select_out_upinfo_eachedge
    select_out_gapinfo[[tumor_name]] = temp_select_out_gapinfo
    select_out_gapinfo_eachedge[[tumor_name]] = temp_select_out_gapinfo_eachedge
    
  }
  
  # save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc1/cyc_successors_result"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  ## up
  res_file_path_upinfo = file.path(res_path, res_sub_path, "select_out_upinfo.RData")
  res_file_path_upinfo_both = file.path(res_path, res_sub_path, "select_out_upinfo_both.RData")
  res_file_path_upinfo_eachedge = file.path(res_path, res_sub_path, "select_out_upinfo_eachedge.RData")
  ## gap
  res_file_path_gapinfo = file.path(res_path, res_sub_path, "select_out_gapinfo.RData")
  res_file_path_gapinfo_eachedge = file.path(res_path, res_sub_path, "select_out_gapinfo_eachedge.RData")
  
  ## up
  save(select_out_upinfo, file=res_file_path_upinfo)
  save(select_out_upinfo_both, file=res_file_path_upinfo_both)
  save(select_out_upinfo_eachedge, file=res_file_path_upinfo_eachedge)
  ## gap
  save(select_out_gapinfo, file=res_file_path_gapinfo)
  save(select_out_gapinfo_eachedge, file=res_file_path_gapinfo_eachedge)
  
}












# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# 
# my_ug_degnode_1_main(output_path, res_path, input_tumor_name)

























