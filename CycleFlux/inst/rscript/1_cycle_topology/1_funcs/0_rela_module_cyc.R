
##################################### part1 module_cycid ################################################
#############################################################################################


#############################################################################################
get_cyc_elements_df<-function(cycle_directed) {
  cyc_elements_rct_enzy = cycle_directed[,c("cycid","ord_cpds_str")]
  for (i in 1:length(cyc_elements_rct_enzy[,1])) {
    # rid string getter
    temp_cycid = cyc_elements_rct_enzy[i,1]
    temp_reaction_array = cyc_elements_rct_enzy[i,2]
    temp_reaction_array = substring(temp_reaction_array, 2,nchar(temp_reaction_array)-1)
    temp_reaction_array <- unlist(strsplit(temp_reaction_array,split = ","))[1]
    temp_reaction_array = substring(temp_reaction_array, 2,nchar(temp_reaction_array)-1)
    temp_reaction_array <- unlist(strsplit(temp_reaction_array,split = "->"))

    temp_cycrct_array = list()
    for (j in 1:length(temp_reaction_array)) {
      if(j%%2 == 0) {
        temp_cycrct_array = append(temp_cycrct_array,temp_reaction_array[j])
      }
    }
    # all reactions of the cycle
    temp_cycrct_array = unlist(temp_cycrct_array)

    temp_cycrct = paste("'",temp_cycrct_array,"'", sep = "")
    temp_cycrct = paste(temp_cycrct, collapse  = ", ")
    temp_cycrct = paste0("{",temp_cycrct,"}")

    cyc_elements_rct_enzy[i, "rids"] = temp_cycrct

  }

  cyc_elements_df = cyc_elements_rct_enzy[,c(1,3)]
  return(cyc_elements_df)
}

#############################################################################################
# deal module string
deal_module_string_not_unique <- function(cycle_rela_module, colname) {
  for (i in 1:length(rownames(cycle_rela_module))) {
    module_str = cycle_rela_module[i,colname]
    module_str = substring(module_str, 2,)
    module_str <- unlist(strsplit(module_str,split = ";"))
    #module_str = unique(module_str)
    module_str = paste(module_str, collapse  = ";")
    cycle_rela_module[i,colname] = module_str
  }
  return(cycle_rela_module)
}

#############################################################################################
# get module by reaction
get_cycle_rela_module_by_rct<-function(cycle_rela_module) {

  cycle_rela_module[, "reaction"] = ""
  for (i in 1:length(rownames(cycle_rela_module))) {
    temp_cyc_rids = cycle_rela_module[i,"rids"]
    temp_cyc_rids = substring(temp_cyc_rids, 2,nchar(temp_cyc_rids)-1)
    temp_cyc_rids <- unlist(strsplit(temp_cyc_rids,split = ", "))
    temp_rid_cell = ""
    for (j in 1:length(temp_cyc_rids)) {
      temp_rid = substring(temp_cyc_rids[j], 2,nchar(temp_cyc_rids[j])-1)

      temp_rid_cell = paste0(temp_rid_cell, ";", temp_rid)
    }
    cycle_rela_module[i, "reaction"] = temp_rid_cell
  }
  return(cycle_rela_module)
}









#############################################################################################
# test
rela_module_cyc_main <- function(output_path, res_path) {
  load(file.path(output_path, "cycle_directed.RData"))

  cycle_elements_df = get_cyc_elements_df(cycle_directed)
  cycle_rela_module = get_cycle_rela_module_by_rct(cycle_elements_df)
  cycle_rela_module = deal_module_string_not_unique(cycle_rela_module,"reaction")
  cycle_rela_module = cycle_rela_module[,c("cycid", "reaction")]


  # save
  res_sub_path = "1_cycle_topology/result_topo"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_rela_module.RData")

  #save(cycle_rela_module, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cycle_rela_module.RData")
  save(cycle_rela_module, file=res_file_path)

}





# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# rela_module_cyc_main(output_path, res_path)



