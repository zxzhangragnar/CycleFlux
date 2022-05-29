
#e2<--up--c2<--up--c1
#outdegree
get_select_out_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"

  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycle_id %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  #  e1<----c2<--up--c1
  select_out_upinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycle_id = ug_cyc_pct_upnode[i, "cycle_id"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]

    tmp_row1 = out_gapinfo[which((out_gapinfo$cycle_id == tmp_cycle_id) & (out_gapinfo$node == tmp_cout)),]
    select_out_upinfo = rbind(select_out_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_out_upinfo = select_out_upinfo[which(select_out_upinfo$ifup == "up"),]
  return(select_out_upinfo)
}



get_select_in_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycle_id %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$ifup == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_upinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycle_id = ug_cyc_pct_upnode[i, "cycle_id"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycle_id == tmp_cycle_id) & (in_gapinfo$node == tmp_cin)),]
    select_in_upinfo = rbind(select_in_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_upinfo = select_in_upinfo[which(select_in_upinfo$ifup == "up"),]

  return(select_in_upinfo)
}



#####################################################################################
my_ug_degnode_1_main <- function(output_path, res_path, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_out_upinfo_eachedge = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ##up
    temp_select_out_upinfo_eachedge = get_select_out_upinfo_eachedge(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)

    select_out_upinfo_eachedge[[tumor_name]] = temp_select_out_upinfo_eachedge
  }

  return(select_out_upinfo_eachedge)
}



my_ug_degnode_2_main <- function(output_path, res_path, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_in_upinfo_eachedge = list()

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ## up
    temp_select_in_upinfo_eachedge = get_select_in_upinfo_eachedge(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)

    select_in_upinfo_eachedge[[tumor_name]] = temp_select_in_upinfo_eachedge
  }

  return(select_in_upinfo_eachedge)
}








# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
#
# my_ug_degnode_1_main(output_path, res_path, input_tumor_name)

























