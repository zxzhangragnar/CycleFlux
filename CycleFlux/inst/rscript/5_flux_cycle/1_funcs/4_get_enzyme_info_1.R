

##################################################################################
##################################################################################
# enzy 规则 造规则
# enzy找：1.有多少酶中全部基因 up down
# 2.有多少酶表达量最高的gene up down

# 给cycle_edge_expression增加一条 gene_count
# 将tumor_count\normal_count\tumor_FC添加到 tumor_cyc_enzyme_genecount_updown 中
get_cyc_enzyme_genecount<-function(tumor_name, tumor_gene_count, normal_gene_count, hsa_net, cycle_edge_expression, gene_foldchange_list, cyc_gene_not_in_tgca_but_in_tofind) {

  tumor_gene_FC = gene_foldchange_list[[tumor_name]]
  rownames_tumor_gene_FC = rownames(tumor_gene_FC)


  gene_count_df = tumor_gene_count[[tumor_name]]
  normal_gene_df = normal_gene_count[[tumor_name]]
  gene_FC_df = as.data.frame(tumor_gene_FC[,"foldchange"])
  colnames(gene_FC_df) = c("foldchange")
  rownames(gene_FC_df) = rownames_tumor_gene_FC

  cyc_enzyme_genecount = data.frame()
  for (i in 1:length(cycle_edge_expression[,1])) {
    #这个环的这个反应(相同cin,cout)所对应的所有酶
    temp_cin = cycle_edge_expression[i,"c_in"]
    temp_cout = cycle_edge_expression[i,"c_out"]

    temp_enzyme_array = cycle_edge_expression[i,"enzyme"]
    temp_enzyme_array <- unlist(strsplit(temp_enzyme_array,split = ";"))

    for (j in 1:length(temp_enzyme_array)) {
      temp_enzyme = temp_enzyme_array[j]
      #这个酶所对应的所有基因
      temp_gene_array = as.character(unlist(hsa_net[which((hsa_net$EC == temp_enzyme) & ((hsa_net$C_in == temp_cin)|(hsa_net$C_in == temp_cout)) & ((hsa_net$C_out == temp_cin)|(hsa_net$C_out == temp_cout))), "Gene_symbol"]))

      temp_gene_array = unique(temp_gene_array)
      temp_gene_array <- unlist(strsplit(temp_gene_array,split = ";"))
      for (k in 1:length(temp_gene_array)) {
        temp_gene = temp_gene_array[k]
        temp_gene_old = temp_gene
        if (temp_gene %in% cyc_gene_not_in_tgca_but_in_tofind) {
          temp_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_gene),2][[1]]
        }
        temp_row = data.frame()
        temp_row[1,"cycid"] = cycle_edge_expression[i,1]
        temp_row[1,"rid"] = cycle_edge_expression[i,2]
        temp_row[1,"enzyme"] = temp_enzyme
        temp_row[1,"gene_name"] = temp_gene_old

        temp_row[1,"gene_count"] = gene_count_df[c(temp_gene),]
        temp_row[1,"normal_gene_count"] = normal_gene_df[c(temp_gene),]
        temp_row[1,"FoldChange"] = gene_FC_df[c(temp_gene),]

        cyc_enzyme_genecount = rbind(cyc_enzyme_genecount, temp_row)
      }

    }
  }
  print(tumor_name)
  return(cyc_enzyme_genecount)
}






#############################################################
# 计算环中某个酶ec的某个基因gene的sign 上调下调值
get_tumor_gap_enzyme_genecount_updowntcga<-function(tumor_name, gene_foldchange_list, datatype, tumor_cyc_enzyme_genecount, tumor_cyc_enzyme_genecount_pure, cyc_gene_not_in_tgca_but_in_tofind)
{

  stat_r = as.data.frame(gene_foldchange_list[tumor_name])
  names(stat_r) = c("gene_id","log2FoldChange","up_down","sign_up_down","pvalue")

  #based on tumor_cyc_enzyme_genecount
  if (datatype == "pure") {
    cyc_enzyme_genecount_df = tumor_cyc_enzyme_genecount_pure[[tumor_name]]
  }else {
    cyc_enzyme_genecount_df = tumor_cyc_enzyme_genecount[[tumor_name]]
  }
  cyc_up_down = cyc_enzyme_genecount_df
  cyc_up_down[,"sign"] = 0
  for (i in 1:length(cyc_up_down[,1])) {
    temp_str_gene = cyc_up_down[i, "gene_name"]

    if (temp_str_gene %in% cyc_gene_not_in_tgca_but_in_tofind) {
      temp_str_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_str_gene),2][[1]]
    }
    genesign = as.double(stat_r[temp_str_gene,"sign_up_down"])
    cyc_up_down[i,"sign"] = genesign

  }
  cyc_up_down[which(cyc_up_down$sign > 0),"updown"] = "up"
  cyc_up_down[which(cyc_up_down$sign < 0),"updown"] = "down"

  return(cyc_up_down)
}

#########################################




#############################################################
# 将每个酶的各个gene按count大小排序
order_tumor_gap_enzyme_genecount_updowntcga<-function(tumor_name, tumor_cyc_enzyme_genecount_updown_pure)
{
  #new data frame
  count_pure_temp_ordered = data.frame()
  count_pure_temp = tumor_cyc_enzyme_genecount_updown_pure[[tumor_name]]

  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i]
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
    temp_enzyme_df = temp_enzyme_df[order(temp_enzyme_df$gene_count, decreasing = T),]
    #add to new data frame
    count_pure_temp_ordered = rbind(count_pure_temp_ordered, temp_enzyme_df)
  }

  return(count_pure_temp_ordered)
}







#############################################################
# 分成几类type:
# 1.酶中所有基因gene全上调，全下调
# 2.酶中表达量最高的基因gene上调，下调

#get type function
get_tumor_gap_enzyme_genecount_updowntcga_type<-function(tumor_name, tumor_cyc_count_pure_temp_ordered) {

  ## highstcount_gene
  highstcount_gene_up = data.frame()
  highstcount_gene_down = data.frame()
  count_pure_temp = tumor_cyc_count_pure_temp_ordered[[tumor_name]]
  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i]
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
    if (!is.na(temp_enzyme_df[1,"updown"])) { #BCMO1 NA
      if (temp_enzyme_df[1,"updown"] == "up") { # row 1 is highstcount
        #add to new data frame
        highstcount_gene_up = rbind(highstcount_gene_up, temp_enzyme_df)
      }else if(temp_enzyme_df[1,"updown"] == "down") {
        #add to new data frame
        highstcount_gene_down = rbind(highstcount_gene_down, temp_enzyme_df)
      }
    }
  }

  ## all_gene
  all_gene_up = data.frame()
  all_gene_down = data.frame()
  count_pure_temp = tumor_cyc_count_pure_temp_ordered[[tumor_name]]
  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i]
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
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







################################### merge part2 ############################################
###############################################################################



merge_enzyme_count_ptype<-function(tumor_name, tumor_cyc_enzyme_genecount_updown, tumor_cyc_count_pure_type) {
  temp_ec_count = tumor_cyc_enzyme_genecount_updown[[tumor_name]]
  temp_ec_ptype = tumor_cyc_count_pure_type[[tumor_name]]

  temp_ec_count[,"ec_situation"] = 0
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["highstcount_gene_up"]]),"ec_situation"] = "hup"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["highstcount_gene_down"]]),"ec_situation"] = "hdown"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["all_gene_up"]]),"ec_situation"] = "aup"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["all_gene_down"]]),"ec_situation"] = "adown"

  return(temp_ec_count)
}


#在ec中  在aup adown中 增加上  每个 enzyme 多少gene up 多少gene down
add_enzyme_count_gene_num <- function(tumor_name, tumor_cyc_enzyme_count_type) {
  #tumor_name = "COAD"
  tumor_enzyme_info = tumor_cyc_enzyme_count_type[[tumor_name]]

  for (i in 1:length(rownames(tumor_enzyme_info))) {
    temp_sit = tumor_enzyme_info[i,"ec_situation"]
    temp_cyc = tumor_enzyme_info[i,"cycid"]
    temp_rid = tumor_enzyme_info[i,"rid"]
    temp_ezy = tumor_enzyme_info[i,"enzyme"]

    temp_df = tumor_enzyme_info[which(tumor_enzyme_info$cycid == temp_cyc & tumor_enzyme_info$rid == temp_rid & tumor_enzyme_info$enzyme == temp_ezy),]
    ec_gene_num = length(rownames(temp_df))
    ec_gene_up_num = 0
    ec_gene_down_num = 0

    ec_gene_up_df = temp_df[which(temp_df$updown == "up"),]
    ec_gene_up_num = length(rownames(ec_gene_up_df))
    ec_gene_down_df = temp_df[which(temp_df$updown == "down"),]
    ec_gene_down_num = length(rownames(ec_gene_down_df))


    tumor_enzyme_info[i, "ec_gene_num"] = ec_gene_num
    tumor_enzyme_info[i, "ec_gene_up_num"] = ec_gene_up_num
    tumor_enzyme_info[i, "ec_gene_down_num"] = ec_gene_down_num
  }
  print(tumor_name)
  return(tumor_enzyme_info)
}






get_enzyme_info_1_main <- function(input_net_file, output_path, res_path, package_path, input_tumor_name) {

  #tumor_gene_count normal_gene_count
  load(file.path(output_path, "cycle_edge_expression.RData"))
  load(file.path(input_net_file))
  load(file.path(res_path, "5_flux_cycle/result_enzyme/tumor_normal_gene_count.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_foldchange_list.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))

  library(readr)
  to_find_correct_symbol_v1_1 = read_csv(file.path(package_path, "tool_data/to_find_correct_symbol_v1.1.csv"), col_names = FALSE)

  # test
  cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]
  #tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
  tumors_array = c(input_tumor_name)


  tumor_cyc_enzyme_genecount = list()
  tumor_cyc_enzyme_genecount_pure = list()

  tumor_cyc_enzyme_genecount_updown = list()
  tumor_cyc_enzyme_genecount_updown_pure = list()

  tumor_cyc_count_pure_temp_ordered = list()
  tumor_cyc_count_pure_type = list()

  tumor_cyc_enzyme_count_type = list()

  for (i in 1:length(tumors_array)) {
    #保存到tumor字典
    enzyme_genecount = as.data.frame(get_cyc_enzyme_genecount(tumors_array[i], tumor_gene_count, normal_gene_count, hsa_net, cycle_edge_expression, gene_foldchange_list, cyc_gene_not_in_tgca_but_in_tofind))
    enzyme_genecount_pure = unique(enzyme_genecount[,c("enzyme","gene_name","gene_count","normal_gene_count","FoldChange")])
    tumor_cyc_enzyme_genecount[[tumors_array[i]]] = enzyme_genecount
    tumor_cyc_enzyme_genecount_pure[[tumors_array[i]]] = enzyme_genecount_pure

    ec_count_updown = get_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i], gene_foldchange_list, "normal", tumor_cyc_enzyme_genecount, tumor_cyc_enzyme_genecount_pure, cyc_gene_not_in_tgca_but_in_tofind)
    tumor_cyc_enzyme_genecount_updown[[tumors_array[i]]] = ec_count_updown
    ec_count_updown_pure = get_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i], gene_foldchange_list, "pure", tumor_cyc_enzyme_genecount, tumor_cyc_enzyme_genecount_pure, cyc_gene_not_in_tgca_but_in_tofind)
    tumor_cyc_enzyme_genecount_updown_pure[[tumors_array[i]]] = ec_count_updown_pure

    count_pure_temp_ordered = order_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i], tumor_cyc_enzyme_genecount_updown_pure)
    tumor_cyc_count_pure_temp_ordered[[tumors_array[i]]] = count_pure_temp_ordered

    type_result = get_tumor_gap_enzyme_genecount_updowntcga_type(tumors_array[i], tumor_cyc_count_pure_temp_ordered)
    tumor_cyc_count_pure_type[[tumors_array[i]]] = type_result

    temp_ec_count = as.data.frame(merge_enzyme_count_ptype(tumors_array[i], tumor_cyc_enzyme_genecount_updown, tumor_cyc_count_pure_type))
    tumor_cyc_enzyme_count_type[[tumors_array[i]]] = temp_ec_count

    tumor_enzyme_info = as.data.frame(add_enzyme_count_gene_num(tumors_array[i], tumor_cyc_enzyme_count_type))
    tumor_cyc_enzyme_count_type[[tumors_array[i]]] = tumor_enzyme_info

  }


  #save
  res_sub_path = "5_flux_cycle/result_enzyme"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "tumor_cyc_enzyme_count_type.RData")

  #save(tumor_cyc_enzyme_count_type, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
  save(tumor_cyc_enzyme_count_type, file=res_file_path)

}




# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
#
# input_tumor_name = "COAD"
# get_enzyme_info_1_main(input_net_file, output_path, res_path, package_path, input_tumor_name)






