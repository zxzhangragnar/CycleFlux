


##################################### part1 build graph ################################################
#############################################################################################

get_cyc_all_distance_df <- function(g, cycle_directed) {
  cyc_all_distance_df = data.frame()
  
  for (i in 1:length(cycle_directed[,1])) {
    cycle_1_id = cycle_directed[i, "cycle_id"]
    compound_1_chain = cycle_directed[i, "compound_chain"]
    
    compound_1_chain = unlist(strsplit(compound_1_chain, split = ";"))
    cycle_1_node = unique(unlist(strsplit(compound_1_chain[1], split = "->")))
    
    for (j in 1:length(cycle_directed[,1])) {
      cycle_2_id = cycle_directed[j, "cycle_id"]
      compound_2_chain = cycle_directed[j, "compound_chain"]
      
      compound_2_chain = unlist(strsplit(compound_2_chain, split = ";"))
      cycle_2_node = unique(unlist(strsplit(compound_2_chain[1], split = "->")))
      
      distance_cycle = min(distances(g, cycle_1_node, cycle_2_node))
      
      temp_row = data.frame()
      temp_row[1, "from"] = cycle_1_id
      temp_row[1, "to"] = cycle_2_id
      temp_row[1, "distance"] = distance_cycle
      
      cyc_all_distance_df = rbind(cyc_all_distance_df, temp_row)
    }
    
  }
  
  
  return(cyc_all_distance_df)
}


#############################################################################################
# test
cycle_plot_and_distance_main <- function(output_path, res_path) {
  load(file.path(output_path, "g.RData"))
  load(file.path(output_path, "cycle_directed.RData"))
  cyc_all_distance_df = get_cyc_all_distance_df(g, cycle_directed)
  
  # save
  res_sub_path = "1_cycle_topology/result_topo"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cyc_all_distance_df.RData")
  
  #save(cyc_all_distance_df, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
  save(cyc_all_distance_df, file=res_file_path)
  
}



# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# cycle_plot_and_distance_main(output_path, res_path)





