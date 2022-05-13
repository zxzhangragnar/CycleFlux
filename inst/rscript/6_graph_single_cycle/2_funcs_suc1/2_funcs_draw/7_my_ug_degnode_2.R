
#情况2:  c2<--up--c1<--up--e1
#主要负责入度
####################################################################

get_select_in_upinfo <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]

  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_upinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycid == tmp_cycid) & (in_gapinfo$node == tmp_cin)),]
    select_in_upinfo = rbind(select_in_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_upinfo = select_in_upinfo[which(select_in_upinfo$ifup == "up"),]

  return(select_in_upinfo)
}

get_select_in_upinfo_both <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]

  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_upinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycid == tmp_cycid) & ((in_gapinfo$node == tmp_cin) | (in_gapinfo$node == tmp_cout))),]
    select_in_upinfo = rbind(select_in_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_upinfo = select_in_upinfo[which(select_in_upinfo$ifup == "up"),]

  return(select_in_upinfo)
}




get_select_in_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_upinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycid == tmp_cycid) & (in_gapinfo$node == tmp_cin)),]
    select_in_upinfo = rbind(select_in_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_upinfo = select_in_upinfo[which(select_in_upinfo$ifup == "up"),]

  return(select_in_upinfo)
}


##################################################################################
### 以下为 gap 的内容

get_select_in_gapinfo <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]

  #  c2<--up--c1
  ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifgap == "gap"),]
  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_gapinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycid == tmp_cycid) & (in_gapinfo$node == tmp_cin)),]
    select_in_gapinfo = rbind(select_in_gapinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_gapinfo = select_in_gapinfo[which(select_in_gapinfo$ifgap == "gap"),]

  return(select_in_gapinfo)
}


get_select_in_gapinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycid %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifgap == "gap"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  #  e1<----c2<--up--c1
  select_in_gapinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycid = ug_cyc_pct_upnode[i, "cycid"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycid == tmp_cycid) & (in_gapinfo$node == tmp_cin)),]
    select_in_gapinfo = rbind(select_in_gapinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_gapinfo = select_in_gapinfo[which(select_in_gapinfo$ifgap == "gap"),]

  return(select_in_gapinfo)
}



#####################################################################################
my_ug_degnode_2_main <- function(output_path, res_path, input_tumor_name) {

  # init
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_in.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))


  ##
  tumors_array = c(input_tumor_name)
  select_in_upinfo = list()
  select_in_upinfo_both = list()
  select_in_upinfo_eachedge = list()
  select_in_gapinfo = list()
  select_in_gapinfo_eachedge = list()

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ## up
    temp_select_in_upinfo = get_select_in_upinfo(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_in_upinfo_both = get_select_in_upinfo_both(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_in_upinfo_eachedge = get_select_in_upinfo_eachedge(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)
    ## gap
    temp_select_in_gapinfo = get_select_in_gapinfo(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)
    temp_select_in_gapinfo_eachedge = get_select_in_gapinfo_eachedge(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)

    select_in_upinfo[[tumor_name]] = temp_select_in_upinfo
    select_in_upinfo_both[[tumor_name]] = temp_select_in_upinfo_both
    select_in_upinfo_eachedge[[tumor_name]] = temp_select_in_upinfo_eachedge
    select_in_gapinfo[[tumor_name]] = temp_select_in_gapinfo
    select_in_gapinfo_eachedge[[tumor_name]] = temp_select_in_gapinfo_eachedge
  }



  # save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc1/cyc_successors_result"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)

  ## up
  res_file_path_upinfo = file.path(res_path, res_sub_path, "select_in_upinfo.RData")
  res_file_path_upinfo_both = file.path(res_path, res_sub_path, "select_in_upinfo_both.RData")
  res_file_path_upinfo_eachedge = file.path(res_path, res_sub_path, "select_in_upinfo_eachedge.RData")
  ## gap
  res_file_path_gapinfo = file.path(res_path, res_sub_path, "select_in_gapinfo.RData")
  res_file_path_gapinfo_eachedge = file.path(res_path, res_sub_path, "select_in_gapinfo_eachedge.RData")


  ## up
  save(select_in_upinfo, file=res_file_path_upinfo)
  save(select_in_upinfo_both, file=res_file_path_upinfo_both)
  save(select_in_upinfo_eachedge, file=res_file_path_upinfo_eachedge)
  ## gap
  save(select_in_gapinfo, file=res_file_path_gapinfo)
  save(select_in_gapinfo_eachedge, file=res_file_path_gapinfo_eachedge)

}






# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
#
# my_ug_degnode_2_main(output_path, res_path, input_tumor_name)





