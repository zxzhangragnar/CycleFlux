##############################################################
# cycle_sum_stat
#
##############################################################
get_cycle_sum_stat <- function(tumor_name, all_gene_stat, cycle_expression, gene_missing_list)
{
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  temp_all_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]

  cycle_stat = as.data.frame(cbind(cycle_expression[,"cycid"],0))
  for (i in 1:length(cycle_expression[,"cycid"])) {
    gene_num_arr = cycle_expression[i,"gene_express"]
    gene_num_arr = substring(gene_num_arr, 2,nchar(gene_num_arr)-1)
    gene_num_arr <- unlist(strsplit(gene_num_arr,split = ", "))
    temp_cyc_allgenenum_sum = 0
    temp_cyc_genesign_sum = 0
    for (j in 1:length(gene_num_arr)) {
      temp_str = unlist(strsplit(gene_num_arr[j],split = ": "))
      temp_str_gene = substring(temp_str[1], 2,nchar(temp_str[1])-1)
      temp_str_num = as.integer(temp_str[2])
      if (!(temp_str_gene %in% missing_genes)) {
        temp_genesign = as.double(temp_all_gene_stat[temp_str_gene,"sign"])*as.double(temp_str_num)
        temp_cyc_allgenenum_sum = temp_cyc_allgenenum_sum + temp_str_num
        temp_cyc_genesign_sum = as.double(temp_cyc_genesign_sum) + temp_genesign #temp_genesign
      }
      
    }
    cycle_stat[i,2] = round(temp_cyc_genesign_sum/temp_cyc_allgenenum_sum, 2) #各个基因的sign总和/基因出现的总次数
  }
  colnames(cycle_stat) = c("cycid", "updownsign")
  
  return(cycle_stat)
}


get_cycle_obvs <- function(tumor_name, cycle_edge_flux_list) {
  temp_cycle_edge_obvs = cycle_edge_flux_list[[tumor_name]]
  
  for (i in 1:length(temp_cycle_edge_obvs[,"cycid"])) {
    if (temp_cycle_edge_obvs[i,"DE_cof"] > 0) {
      temp_cycle_edge_obvs[i,"ifDE"] = 1
    }else {
      temp_cycle_edge_obvs[i,"ifDE"] = 0
    }
  }
  
  cycle_obvs = aggregate(temp_cycle_edge_obvs$ifDE, list(temp_cycle_edge_obvs$cycid), mean)
  colnames(cycle_obvs) = c("cycid", "DE_cof")
  cycle_obvs[,"DE_cof"] = round(cycle_obvs[,"DE_cof"], 2)
  
  #obvs_cycid = cycle_obvs[which(cycle_obvs$DE_cof>0.5), "cycid"]
  
  return(cycle_obvs)
}


merge_obvs<-function(tumor_name, cycle_merge_stat, cycle_edge_flux_list) {
  temp_tumor_df = cycle_merge_stat[[tumor_name]]
  cycle_obvs = get_cycle_obvs(tumor_name, cycle_edge_flux_list)
  temp_tumor_df[,"DE_cof"] = cycle_obvs$DE_cof
  return(temp_tumor_df)
}

merge_struct<-function(tumor_name, cycle_merge_stat, struct_cycle_sort_list) {
  #struct
  temp_tumor_df = cycle_merge_stat[[tumor_name]]
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["net"]]),"stcid"] = 1
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["hub"]]),"stcid"] = 2
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["idv"]]),"stcid"] = 3
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["nrl"]]),"stcid"] = 4
  return(temp_tumor_df)
}





cyc_info_1_main <- function(output_path, res_path, package_path, input_tumor_name) {
  
  load(file.path(output_path, "cycle_expression.RData"))
  load(file.path(package_path, "/tool_data/TCGA_upgap_genes.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/struct_cycle_sort_list.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  
  #merge
  cycle_merge_stat = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    temp_cycle_stat = get_cycle_sum_stat(tumor_name, stat_all, cycle_expression, gene_missing_list)
    cycle_merge_stat[[tumors_array[i]]] = temp_cycle_stat
    
    temp_tumor_df = merge_obvs(tumor_name, cycle_merge_stat, cycle_edge_flux_list)
    cycle_merge_stat[[tumors_array[i]]] = temp_tumor_df
    
    temp_tumor_df = merge_struct(tumors_array[i], cycle_merge_stat, struct_cycle_sort_list)
    cycle_merge_stat[[tumors_array[i]]] = temp_tumor_df
  }
  
  #save
  res_sub_path = "5_flux_cycle/result_enzyme"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_merge_stat.RData")
  
  save(cycle_merge_stat, file=res_file_path)
}




# 
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
# 
# cyc_info_1_main(output_path, res_path, package_path, input_tumor_name)





