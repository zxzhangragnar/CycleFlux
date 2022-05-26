# gene表达量高的cycle take care obvs
# gene不表达的cycle not obvs
#########################################################################


get_cycle_edge_obvs <- function(tumor_name, all_gene_stat, cycle_edge_expression, gene_missing_list, prm_1=0.05, prm_2=1) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  temp_all_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]
  
  cycle_edge_obvs = as.data.frame(cycle_edge_expression[,"cycle_id"])
  colnames(cycle_edge_obvs) = c("cycle_id")
  for (i in 1:length(cycle_edge_expression[,"cycle_id"])) {
    gene_str = cycle_edge_expression[i,"gene_symbol"]
    gene_arr = unlist(strsplit(gene_str,split = ";"))
    
    gene_num = 0
    obvs_num = 0
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      temp_pvalue = 0
      temp_fc = 0
      if (!(temp_gene %in% missing_genes)) {
        temp_pvalue = temp_all_gene_stat[temp_gene,"p.value"]
        temp_fc = temp_all_gene_stat[temp_gene,"FC"]
        gene_num = gene_num + 1
      }
      
      if ((temp_pvalue < prm_1) & (abs(temp_fc) > prm_2)) {
        obvs_num = obvs_num + 1
      }
    }
    
    obvs_rat = obvs_num/gene_num
    # cycle_edge_obvs[i,"DE_num"] = obvs_num
    # cycle_edge_obvs[i,"gene_num"] = gene_num
    
    cycle_edge_obvs[i,"DE_cof"] = round(obvs_rat, 2)
    
  }
  return(cycle_edge_obvs)
}


get_cycle_obvs <- function(tumor_name, cycle_edge_obvs_list, prm_1=0.5) {
  temp_cycle_edge_obvs = cycle_edge_obvs_list[[tumor_name]]
  
  for (i in 1:length(temp_cycle_edge_obvs[,"cycle_id"])) {
    if (temp_cycle_edge_obvs[i,"DE_cof"] > 0) {
      temp_cycle_edge_obvs[i,"ifDE"] = 1
    }else {
      temp_cycle_edge_obvs[i,"ifDE"] = 0
    } 
  }
  
  cycle_obvs = aggregate(temp_cycle_edge_obvs$ifDE, list(temp_cycle_edge_obvs$cycle_id), mean)
  colnames(cycle_obvs) = c("cycle_id", "DE_cof")
  cycle_obvs[,"DE_cof"] = round(cycle_obvs[,"DE_cof"], 2)
  
  #obvs_cycle_id = cycle_obvs[which(cycle_obvs$DE_cof>0.5), "cycle_id"]
  cycle_obvs[,"DE"] = "normal"
  cycle_obvs[which(cycle_obvs$DE_cof>prm_1),"DE"] = "obvious"
  
  return(cycle_obvs)
}


#############

cycle_edge_obvs_main <- function(output_path, res_path, package_path, input_tumor_name) {
  
  load(file.path(output_path, "cycle_edge_expression.RData"))
  load(file.path(package_path, "/tool_data/TCGA_upgap_genes.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))
  
  tumors_array = c(input_tumor_name)
  
  cycle_edge_obvs_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_edge_obvs = get_cycle_edge_obvs(tumor_name, stat_all, cycle_edge_expression, gene_missing_list)
    cycle_edge_obvs_list[[tumor_name]] = cycle_edge_obvs
  }
  cycle_obvs_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_obvs = get_cycle_obvs(tumor_name, cycle_edge_obvs_list)
    cycle_obvs_list[[tumor_name]] = cycle_obvs
  }

  
  #save
  
}










# test
output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
input_tumor_name = "COAD"

obvs_cycles_main(output_path, res_path, package_path, input_tumor_name)






