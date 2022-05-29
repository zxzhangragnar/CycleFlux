
get_never_considered_comp<-function(never_considered_comp, compounds_dict) {
  never_considered_comp_arr = never_considered_comp$x
  compounds_dict_cids = compounds_dict$ENTRY
  never_considered_comp_arr_cids = intersect(never_considered_comp_arr, compounds_dict_cids)
  never_considered_comp_names = never_considered_comp_arr_cids

  return(never_considered_comp_names)
}



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


single_cycle_gap_analysis_main <- function(output_path, res_path, package_path, input_tumor_name) {

  ##
  tumors_array = c(input_tumor_name)
  gapup_cycle_chain_list = get_ug_chain_list(tumors_array, cycle_edge_flux_list, all_chain_list_cid)
  return(gapup_cycle_chain_list)
}

never_considered_comp_names_main <- function(output_path, res_path, package_path, input_tumor_name) {
  # init
  library(readr)
  never_considered_comp = read_csv(file.path("tool_data/never_considered_comp.csv"))

  ##
  never_considered_comp_names = get_never_considered_comp(never_considered_comp, compounds_dict)
  return(never_considered_comp_names)
}


# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
#
# single_cycle_gap_analysis_main(output_path, res_path, package_path, input_tumor_name)

