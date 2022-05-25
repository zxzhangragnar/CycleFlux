

##################################################################
refesh_cycle_flux_list<-function(tumor_name, cycle_flux_list, cycle_edge_flux_list){
  temp_df = cycle_edge_flux_list[[tumor_name]]
  
  gap_c = unique(temp_df[which(temp_df$ifgap=="gap"), "cycle_id"])
  up_c = unique(temp_df[which(temp_df$ifup=="up"), "cycle_id"])

  temp_stat_df = temp_df[,c("cycle_id","mean_fc")]
  temp_mean_fc = aggregate(temp_stat_df$mean_fc, list(temp_stat_df$cycle_id), mean)
  cycle_flux_list[[tumor_name]][, "mean_fc"] = temp_mean_fc$mean_fc
  
  cycle_flux_list[[tumor_name]][, "ifup"] = "normal"
  cycle_flux_list[[tumor_name]][, "ifgap"] = "normal"
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% up_c), "ifup"] = "up"
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% gap_c), "ifgap"] = "gap"
  
  return(cycle_flux_list)
}


cyc_info_gap_main <- function(res_path, input_tumor_name) {
  
  load(file.path(res_path, "/5_flux_cycle/result_final/cycle_flux_list.RData"))
  load(file.path(res_path, "/4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_flux_list = refesh_cycle_flux_list(tumor_name, cycle_flux_list, cycle_edge_flux_list)
  }
  
  #save
  res_sub_path = "5_flux_cycle/result_final"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_flux_list.RData")
  save(cycle_flux_list, file=res_file_path)
  
}



# 
# 
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
# 
# cyc_info_gap_main(res_path, input_tumor_name)







