##################################################################################
##################################################################################
# enzy找：
# 1.哪些酶中全部基因 up down
# 2.哪些酶表达量最高的 gene up down
get_enzyme_stat <- function(tumor_name, all_gene_stat, hsa_net, cycle_edge_expression, gene_missing_list) {
  temp_gene_stat = all_gene_stat[[paste0("TCGA-",tumor_name)]]
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  
  cyc_enzyme_genecount = data.frame()
  for (i in 1:length(rownames(cycle_edge_expression))) {
    #这个环的这条边(相同cin,cout)所对应的所有酶
    temp_cin = cycle_edge_expression[i,"c_in"]
    temp_cout = cycle_edge_expression[i,"c_out"]
    
    temp_enzyme_array = cycle_edge_expression[i,"enzyme"]
    temp_enzyme_array = unlist(strsplit(temp_enzyme_array,split = ";"))
    
    for (j in 1:length(temp_enzyme_array)) {
      temp_enzyme = temp_enzyme_array[j]
      #这个酶所对应的所有基因
      temp_gene_array = as.character(unlist(hsa_net[which((hsa_net$EC == temp_enzyme) & ((hsa_net$C_in == temp_cin)|(hsa_net$C_in == temp_cout)) & ((hsa_net$C_out == temp_cin)|(hsa_net$C_out == temp_cout))), "Gene_symbol"]))

      temp_gene_array = unique(temp_gene_array)
      temp_gene_array = unlist(strsplit(temp_gene_array,split = ";"))
      for (k in 1:length(temp_gene_array)) {
        temp_gene = temp_gene_array[k]
        
        temp_row = data.frame()
        if (!(temp_gene %in% missing_genes)) {
          temp_row[1,"cycid"] = cycle_edge_expression[i, "cycid"]
          temp_row[1,"rid"] = cycle_edge_expression[i, "rid"]
          temp_row[1,"enzyme"] = temp_enzyme
          temp_row[1,"gene_name"] = temp_gene
          
          temp_row[1,"tumor_gene_count"] = temp_gene_stat[c(temp_gene), "mt"]
          temp_row[1,"normal_gene_count"] = temp_gene_stat[c(temp_gene), "mn"]
          temp_row[1,"foldchange"] = temp_gene_stat[c(temp_gene), "FC"]
          
          if (temp_gene_stat[c(temp_gene), "sign"] > 0) {
            temp_row[1,"updown"] = "up"
          }else {
            temp_row[1,"updown"] = "down"
          }
          
        }
        
        cyc_enzyme_genecount = rbind(cyc_enzyme_genecount, temp_row)  
      }
      
    }
  }
  return(cyc_enzyme_genecount)
}


#############################################################                 
# 1.酶中所有基因gene全上调，全下调
# 2.酶中表达量最高的基因gene上调，下调
get_enzyme_gene_condition <- function(tumor_name, enzyme_stat) {
  temp_enzyme_genecount = enzyme_stat[[tumor_name]]
  unique_enzyme_genecount = unique(temp_enzyme_genecount[,c("enzyme","gene_name","tumor_gene_count","normal_gene_count","foldchange")])
  temp_enzyme_genecount = temp_enzyme_genecount[rownames(unique_enzyme_genecount),]
  
  highstcount_gene_up = data.frame()
  highstcount_gene_down = data.frame()
  enzyme_table = unique(temp_enzyme_genecount$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i] 
    temp_enzyme_df = temp_enzyme_genecount[which(temp_enzyme_genecount$enzyme == temp_enzyme),]
    if (!is.na(temp_enzyme_df[1,"updown"])) { #BCMO1 NA
      if (temp_enzyme_df[1,"updown"] == "up") { # row 1 is highstcount
        highstcount_gene_up = rbind(highstcount_gene_up, temp_enzyme_df)      
      }else if(temp_enzyme_df[1,"updown"] == "down") {
        highstcount_gene_down = rbind(highstcount_gene_down, temp_enzyme_df)      
      }    
    }
  }
  
  ## all_gene
  all_gene_up = data.frame()
  all_gene_down = data.frame()
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i] 
    temp_enzyme_df = temp_enzyme_genecount[which(temp_enzyme_genecount$enzyme == temp_enzyme),]
    is_all_up = TRUE
    is_all_down = TRUE
    for (j in 1:length(temp_enzyme_df)) {
      if (!is.na(temp_enzyme_df[j,"updown"])) { #BCMO1 NA
        if (temp_enzyme_df[j,"updown"] != "up") { # row 1 is highstcount
          is_all_up = FALSE
        }else if(temp_enzyme_df[j,"updown"] != "down") {
          is_all_down = FALSE
        }    
      }    
    }
    if (is_all_up) {
      all_gene_up = rbind(all_gene_up, temp_enzyme_df)
    }
    if (is_all_down) {
      all_gene_down = rbind(all_gene_down, temp_enzyme_df)
    }
  }
  
  type_result = list()
  type_result[["highstcount_gene_up"]] = unique(highstcount_gene_up$enzyme)
  type_result[["highstcount_gene_down"]] = unique(highstcount_gene_down$enzyme)
  type_result[["all_gene_up"]] = unique(all_gene_up$enzyme)
  type_result[["all_gene_down"]] = unique(all_gene_down$enzyme)
  
  return(type_result)
  
}







################################### merge part2 #########################################                 
#########################################################################################

#在ec中  在aup adown中 增加上  每个 enzyme 多少gene up 多少gene down
merge_enzyme_situation <- function(tumor_name, enzyme_stat, enzyme_gene_condition) {
  temp_ec_stat = enzyme_stat[[tumor_name]]
  temp_ec_condition = enzyme_gene_condition[[tumor_name]]
  
  temp_ec_stat[,"ec_situation"] = 0
  temp_ec_stat[which(temp_ec_stat$enzyme %in% temp_ec_condition[["highstcount_gene_up"]]),"ec_situation"] = "hup"
  temp_ec_stat[which(temp_ec_stat$enzyme %in% temp_ec_condition[["highstcount_gene_down"]]),"ec_situation"] = "hdown"
  temp_ec_stat[which(temp_ec_stat$enzyme %in% temp_ec_condition[["all_gene_up"]]),"ec_situation"] = "aup"
  temp_ec_stat[which(temp_ec_stat$enzyme %in% temp_ec_condition[["all_gene_down"]]),"ec_situation"] = "adown"
  
  for (i in 1:length(rownames(temp_ec_stat))) {
    temp_cyc = temp_ec_stat[i,"cycid"]
    temp_rid = temp_ec_stat[i,"rid"]
    temp_ezy = temp_ec_stat[i,"enzyme"]
    
    temp_df = temp_ec_stat[which(temp_ec_stat$cycid == temp_cyc & temp_ec_stat$rid == temp_rid & temp_ec_stat$enzyme == temp_ezy),]
    ec_gene_num = length(rownames(temp_df))
    ec_gene_up_num = 0
    ec_gene_down_num = 0
    
    ec_gene_up_df = temp_df[which(temp_df$updown == "up"),]
    ec_gene_up_num = length(rownames(ec_gene_up_df))
    ec_gene_down_df = temp_df[which(temp_df$updown == "down"),]
    ec_gene_down_num = length(rownames(ec_gene_down_df))
    
    
    temp_ec_stat[i, "ec_gene_num"] = ec_gene_num
    temp_ec_stat[i, "ec_gene_up_num"] = ec_gene_up_num
    temp_ec_stat[i, "ec_gene_down_num"] = ec_gene_down_num
  }
  return(temp_ec_stat)  
}






get_enzyme_info_1_main <- function(input_net_file, output_path, res_path, package_path, input_tumor_name) {
  
  load(file.path(output_path, "cycle_edge_expression.RData"))
  load(file.path(input_net_file))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))
  load(file.path(package_path, "/tool_data/TCGA_upgap_genes.RData"))
  
  # test
  tumors_array = c(input_tumor_name)

  enzyme_stat = list()
  enzyme_gene_condition = list()
  cycle_enzyme_stat = list()
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    
    enzyme_stat[[tumor_name]] = as.data.frame(get_enzyme_stat(tumor_name, stat_all, hsa_net, cycle_edge_expression, gene_missing_list))

    enzyme_gene_condition[[tumor_name]] = get_enzyme_gene_condition(tumor_name, enzyme_stat)
    
    cycle_enzyme_stat[[tumor_name]] = as.data.frame(merge_enzyme_situation(tumor_name, enzyme_stat, enzyme_gene_condition))
  }

  #save
  res_sub_path = "5_flux_cycle/result_enzyme"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_enzyme_stat.RData")

  save(cycle_enzyme_stat, file=res_file_path)
}


# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# 
# input_tumor_name = "COAD"
# get_enzyme_info_1_main(input_net_file, output_path, res_path, package_path, input_tumor_name)
# 





