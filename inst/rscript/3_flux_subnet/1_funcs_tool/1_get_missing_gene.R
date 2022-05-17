####################################################################################################


####################################################################################################
# 找到环中所包含的所有基因
get_gene_expressinfo<-function(subnet_edge_expression){
  cyc_gene_express_list = c()
  
  for (i in 1:length(subnet_edge_expression[,1])) {
    gene_num_arr = as.character(subnet_edge_expression[i,"gene_symbol"])
    gene_num_arr = unlist(strsplit(gene_num_arr,split = ";"))
    for (j in 1:length(gene_num_arr)) {
      cyc_gene_express_list = append(cyc_gene_express_list,gene_num_arr[j])
    }
  }
  cyc_gene_express_list = unique(cyc_gene_express_list)
  return(cyc_gene_express_list)
}

####################################################################################################


####################################################################################################
get_cyc_gene_not_in_tgca<-function(cyc_gene_expressinfo, gene_foldchange_list, to_find_correct_symbol_v1_1){
  cyc_gene_kegg = c()
  cyc_gene_tcga = c()
  for (i in 1:length(rownames(to_find_correct_symbol_v1_1))) {
    kegg_gene = to_find_correct_symbol_v1_1[i,1]
    tcga_gene = to_find_correct_symbol_v1_1[i,2]
    if(kegg_gene %in% cyc_gene_expressinfo){
      cyc_gene_kegg = append(cyc_gene_kegg,kegg_gene)
    }
    if(tcga_gene %in% cyc_gene_expressinfo){
      cyc_gene_tcga = append(cyc_gene_tcga,tcga_gene)
    }
  }
  
  cyc_gene_not_in_tgca = c()
  cyc_gene_in_tgca = c()
  tcga_gene_list = rownames(gene_foldchange_list[[1]])
  cyc_gene_expressinfo = unique(cyc_gene_expressinfo)
  for (i in 1:length(cyc_gene_expressinfo)) {
    gene = cyc_gene_expressinfo[i]
    if(gene %in% tcga_gene_list){
      cyc_gene_in_tgca = append(cyc_gene_in_tgca,gene)
    }else{
      cyc_gene_not_in_tgca = append(cyc_gene_not_in_tgca,gene)
    }
  }
  
  return(cyc_gene_not_in_tgca)
}


get_gene_missing_list<-function(cyc_gene_expressinfo, gene_foldchange_list, to_find_correct_symbol_v1_1){
  #test
  cyc_gene_not_in_tgca = get_cyc_gene_not_in_tgca(cyc_gene_expressinfo, gene_foldchange_list, to_find_correct_symbol_v1_1)
  cyc_gene_not_in_tgca_but_in_tofind = intersect(cyc_gene_not_in_tgca,to_find_correct_symbol_v1_1$X1)
  cyc_gene_not_in_tgca_not_in_tofind = setdiff(cyc_gene_not_in_tgca, cyc_gene_not_in_tgca_but_in_tofind)
  
  #result
  gene_missing_list = list()
  gene_missing_list[["cyc_gene_expressinfo"]] = cyc_gene_expressinfo
  gene_missing_list[["cyc_gene_not_in_tgca"]] = cyc_gene_not_in_tgca
  gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]] = cyc_gene_not_in_tgca_but_in_tofind
  gene_missing_list[["cyc_gene_not_in_tgca_not_in_tofind"]] = cyc_gene_not_in_tgca_not_in_tofind
  
  return(gene_missing_list)
}
####################################################################################################



#########################################################################
get_gene_and_tofind_list_main <- function(output_path, res_path, package_path) {

  load(file.path(output_path, "subnet_edge_expression.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_foldchange_list.RData"))
  
  cyc_gene_expressinfo = get_gene_expressinfo(subnet_edge_expression)
  
  library(readr)
  to_find_correct_symbol_v1_1 = read_csv(file.path(package_path, "tool_data/to_find_correct_symbol_v1.1.csv"), col_names = FALSE)
  
  gene_missing_list = get_gene_missing_list(cyc_gene_expressinfo, gene_foldchange_list, to_find_correct_symbol_v1_1)
  
  # save
  res_sub_path = "3_flux_subnet/result_tool"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "gene_missing_list.RData")

  save(gene_missing_list, file=res_file_path)
}





# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# get_gene_and_tofind_list_main(output_path, res_path, package_path)



