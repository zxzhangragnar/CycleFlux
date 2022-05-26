
get_gap_cycleid_freq_list <- function(input_tumor_name, cycle_edge_flux_list) {
  gap_cycleid_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    gap_cycleid = cycle_edge_flux_list[[tumor_name]][which(cycle_edge_flux_list[[tumor_name]]$ifgap == "gap"),]
    gap_cycleid = unique(gap_cycleid)
    gap_cycleid_freq = as.data.frame(table(gap_cycleid$cycle_id))
    
    gap_cycleid_freq = gap_cycleid_freq[order(gap_cycleid_freq$Freq, decreasing = T),]
    
    colnames(gap_cycleid_freq) = c("cycle_id", "gap_freq")
    gap_cycleid_freq_list[[tumor_name]] = gap_cycleid_freq
  }
  return(gap_cycleid_freq_list)
}






get_shift_node_freq_list <- function(input_tumor_name, cycle_shift_path_df_list) {
  shift_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    shift_node_freq = as.data.frame(table(cycle_shift_path_df_list[[tumor_name]]$shift_node))
    shift_node_freq = shift_node_freq[order(shift_node_freq$Freq, decreasing = T),]
    
    shift_node_df = cycle_shift_path_df_list[[tumor_name]]
    if (length(rownames(shift_node_df)) != 0) {
      cycle_ids =  aggregate(shift_node_df$cycle_id, list(shift_node_df$shift_node), paste, collapse = ";")
      rownames(cycle_ids) = cycle_ids[,1]
      
      colnames(shift_node_freq) = c("node", "shift_freq")
      
      for (j in 1:length(shift_node_freq[,1])) {
        shift_node_freq[j,"cycle_id"] = cycle_ids[shift_node_freq[j,"node"],2]
      }
    }
    
    shift_node_freq_list[[tumor_name]] = shift_node_freq
  }
  return(shift_node_freq_list)
}



get_shift_edge_node_freq_list <- function(input_tumor_name, cycle_shift_path_df_list) {
  shift_edge_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    temp_cycle_shift_path_df = cycle_shift_path_df_list[[tumor_name]]
    shift_edgenode_freq = data.frame()
    if (length(rownames(temp_cycle_shift_path_df)) != 0) {
      shift_edgenode_freq = as.data.frame(table(temp_cycle_shift_path_df[,c("endpoint_node","shift_node")]))
      shift_edgenode_freq = shift_edgenode_freq[order(shift_edgenode_freq$Freq, decreasing = T),]
      
      colnames(shift_edgenode_freq) = c("edge_node", "node", "shift_freq")
      
      shift_edgenode_freq = shift_edgenode_freq[which(shift_edgenode_freq$shift_freq != 0),]
    }
    
    shift_edge_node_freq_list[[tumor_name]] = shift_edgenode_freq
  }
  return(shift_edge_node_freq_list)
}




#####################################################################################
# cycle_edge_expression
freq_stat_main <- function(res_path, input_tumor_name) {
  
  #init
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_shift_path_df_list.RData"))
  
  tumors_array = c(input_tumor_name)
  
  ##
  gap_cycleid_freq_list = get_gap_cycleid_freq_list(input_tumor_name, cycle_edge_flux_list)
  shift_node_freq_list = get_shift_node_freq_list(input_tumor_name, cycle_shift_path_df_list)
  shift_edge_node_freq_list = get_shift_edge_node_freq_list(input_tumor_name, cycle_shift_path_df_list)
  
  #save
  res_sub_path = "8_analysis_metric/metric_result"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "freq_statistics.RData")
  
  
  save(gap_cycleid_freq_list, shift_node_freq_list, shift_edge_node_freq_list, file=res_file_path)
  
}



# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# freq_stat_main(res_path, input_tumor_name)





