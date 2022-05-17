########################################### 1 ########################################################
# "mean_fc"
########################################### 2 ########################################################
# "gene_fold_change"
# "gene_pvalue"
# "t_gene_mean_val"
# "n_gene_mean_val"
# "t_gene_0_pct"
# "n_gene_0_pct"
# "t_gene_log_lss1_pct"
# "n_gene_log_lss1_pct"
########################################### 3 ########################################################
# "ifgap"
# "ifup"
#######################################################################################################


#################################### 2 #############################################
####################################################################################
get_mean_fc<-function(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  
  temp_cycle_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]
  temp_cycle_edge_flux[,"mean_fc"] = 0
  
  for (i in 1:length(temp_cycle_edge_flux[,1])) {
    temp_edge_sum_fc = 0
    temp_gene_symbol = temp_cycle_edge_flux[i,"gene_symbol"]
    temp_gene_symbol_arr = unlist(strsplit(temp_gene_symbol,split = ";"))
    
    for (j in 1:length(temp_gene_symbol_arr)) {
      temp_gene = temp_gene_symbol_arr[j]
      if (!(temp_gene %in% missing_genes)) {
        genesign = as.double(temp_gene_stat[temp_gene,"FC"])
        temp_edge_sum_fc = temp_edge_sum_fc + genesign
      }
    }
    temp_cycle_edge_flux[i,"mean_fc"] = temp_edge_sum_fc/length(temp_gene_symbol_arr) #more equal
  }
  
  return(temp_cycle_edge_flux)
}


#################################### 2 #############################################
####################################################################################
add_gene_stat<-function(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  
  temp_cycle_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]
  
  temp_cycle_edge_flux[, "gene_fc"] = NA
  temp_cycle_edge_flux[, "gene_pvalue"] = NA
  temp_cycle_edge_flux[, "gene_mt"] = NA
  temp_cycle_edge_flux[, "gene_mn"] = NA
  
  for (i in 1:length(temp_cycle_edge_flux[,1])) {
    temp_gene_symbol = temp_cycle_edge_flux[i,"gene_symbol"]
    temp_gene_symbol_arr = unlist(strsplit(temp_gene_symbol,split = ";"))
    
    gene_fc_arr = c()
    gene_pvalue_arr = c()
    gene_mt_arr = c()
    gene_mn_arr = c()
    for (j in 1:length(temp_gene_symbol_arr)) {
      temp_gene = temp_gene_symbol_arr[j]
      if (!(temp_gene %in% missing_genes)) {
        tmp_fc = temp_gene_stat[temp_gene,"FC"]
        tmp_pvalue = temp_gene_stat[temp_gene,"p.value"]
        tmp_mt = temp_gene_stat[temp_gene,"mt"]
        tmp_mn = temp_gene_stat[temp_gene,"mn"]
        gene_fc_arr = append(gene_fc_arr, paste0(temp_gene,":",tmp_fc))
        gene_pvalue_arr = append(gene_pvalue_arr, paste0(temp_gene,":",tmp_pvalue))
        gene_mt_arr = append(gene_mt_arr, paste0(temp_gene,":",tmp_mt))
        gene_mn_arr = append(gene_mn_arr, paste0(temp_gene,":",tmp_mn))
      }
      
    }
    
    gene_fc_str = paste(gene_fc_arr, collapse  = ";")
    gene_pvalue_str = paste(gene_pvalue_arr, collapse  = ";")
    gene_mt_str = paste(gene_mt_arr, collapse  = ";")
    gene_mn_str = paste(gene_mn_arr, collapse  = ";")
    
    temp_cycle_edge_flux[i, "gene_fc"] = gene_fc_str
    temp_cycle_edge_flux[i, "gene_pvalue"] = gene_pvalue_str
    temp_cycle_edge_flux[i, "gene_mt"] = gene_mt_str
    temp_cycle_edge_flux[i, "gene_mn"] = gene_mn_str
  }
  
  return(temp_cycle_edge_flux)
}

#################################### 3 #############################################
####################################################################################
get_up_edges <- function(tumor_name, cycle_edge_flux_list, gene_gapup_info) {
  
  temp_gene_gapup_info = gene_gapup_info[[paste0("TCGA-",tumor_name)]]
  temp_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_edge_flux[, "ifup"] = "normal"
  for (i in 1:length(temp_edge_flux[,1])) {
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr = unlist(strsplit(gene_arr, split = ";"))
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if(temp_gene %in% temp_gene_gapup_info$up_genes) {
        temp_edge_flux[i, "ifup"] = "up"
      }
    }
  }
  return(temp_edge_flux)
}


get_gap_edges <- function(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  
  temp_gene_gapup_info = gene_gapup_info[[paste0("TCGA-",tumor_name)]]
  temp_all_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]
  temp_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_edge_flux[, "ifgap"] = "normal"
  for (i in 1:length(temp_edge_flux[,1])) {
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr <- unlist(strsplit(gene_arr, split = ";"))
    
    gene_mt_df = data.frame()
    ## find mt(mean tumor value) gene
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if (temp_gene %in% missing_genes) {
        temp_mt = 0
      }else {
        temp_mt = temp_all_gene_stat[temp_gene, "mt"]
      }
      gene_mt_df[j, "gene"] = temp_gene
      gene_mt_df[j, "mt"] = temp_mt
    }
    
    gene_mt_df = gene_mt_df[order(gene_mt_df$mt,decreasing = TRUE),]
    max_mt_gene = gene_mt_df[1, "gene"]
    if(max_mt_gene %in% temp_gene_gapup_info$gap_genes) {
      temp_edge_flux[i, "ifgap"] = "gap"
    }
  }
  return(temp_edge_flux)
}



## main function
get_cycle_edge_flux_list <- function(edge_expression, tumors_array, all_gene_stat, gene_gapup_info, gene_missing_list) {
  cycle_edge_flux_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_edge_flux_list[[tumor_name]] = edge_expression
    
    cycle_edge_flux_list[[tumor_name]] = get_up_edges(tumor_name, cycle_edge_flux_list, gene_gapup_info)
    cycle_edge_flux_list[[tumor_name]] = get_gap_edges(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list)
    cycle_edge_flux_list[[tumor_name]] = get_mean_fc(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list)
    cycle_edge_flux_list[[tumor_name]] = add_gene_stat(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list)
  }
  return(cycle_edge_flux_list)
}

#####################################################################################
# cycle_successors_2_expression_1 cycle_successors_3_expression_1
cycle_successors_expression_main <- function(output_path, res_path, package_path, input_tumor_name, par1=1, par2=10, par3=-2, par4=1.71) {
  #init
  load(file.path(output_path, "cycle_successors_2_expression_1.RData"))
  load(file.path(output_path, "cycle_successors_2_expression_2.RData"))
  load(file.path(output_path, "cycle_successors_3_expression_1.RData"))
  load(file.path(output_path, "cycle_successors_3_expression_2.RData"))

  load(file.path(res_path, "/tool_data/TCGA_upgap_genes.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))
  
  tumors_array = c(input_tumor_name)
  
  ##
  cycle_edge_flux_list_2_1 = get_cycle_edge_flux_list(cycle_successors_2_expression_1, tumors_array, stat_all, DEG_selected, gene_missing_list)
  cycle_edge_flux_list_2_2 = get_cycle_edge_flux_list(cycle_successors_2_expression_2, tumors_array, stat_all, DEG_selected, gene_missing_list)
  cycle_edge_flux_list_3_1 = get_cycle_edge_flux_list(cycle_successors_3_expression_1, tumors_array, stat_all, DEG_selected, gene_missing_list)
  cycle_edge_flux_list_3_2 = get_cycle_edge_flux_list(cycle_successors_3_expression_2, tumors_array, stat_all, DEG_selected, gene_missing_list)
  
  
  #save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc2/cyc_successors_data"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  
  res_file_path_2_1 = file.path(res_path, res_sub_path, "cycle_edge_flux_list_2_1.RData")
  res_file_path_2_2 = file.path(res_path, res_sub_path, "cycle_edge_flux_list_2_2.RData")
  res_file_path_3_1 = file.path(res_path, res_sub_path, "cycle_edge_flux_list_3_1.RData")
  res_file_path_3_2 = file.path(res_path, res_sub_path, "cycle_edge_flux_list_3_2.RData")
  
  save(cycle_edge_flux_list_2_1, file=res_file_path_2_1)
  save(cycle_edge_flux_list_2_2, file=res_file_path_2_2)
  save(cycle_edge_flux_list_3_1, file=res_file_path_3_1)
  save(cycle_edge_flux_list_3_2, file=res_file_path_3_2)
  
}






# test
output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
input_tumor_name = "COAD"

cycle_successors_expression_main(output_path, res_path, package_path, input_tumor_name)




