##############################################################
# 目的: 
# tumor_cyc_updownsign: 每个环的updown和pvalue
#
##############################################################


############################################# 1 ####################################################
####################################################################################################



############################################# 2 ####################################################
####################################################################################################
#使用获取到的gene_foldchange_list数据,获取cyc_up_down信息
get_tumor_updownsign<-function(tumor_name, gene_foldchange_list, cyc_gene_not_in_tgca, if_filter)
{
  #tumor_name = "COAD"
  #注意要加上as.dataframe,否则无法判断e-3比较大小
  stat_r = as.data.frame(gene_foldchange_list[tumor_name])
  names(stat_r) = c("gene_id","log2FoldChange","up_down","sign_up_down","pvalue")
  et = stat_r
  # if_filter
  if(if_filter) {
    #etSig = et[which(et$pvalue < 0.05 & abs(et$log2FoldChange) > 1),]
    #etSig = et[which(et$pvalue < 0.001),]
    etSig = et
  }else {
    etSig = et
  }
  
  #cyc_up_down
  cyc_up_down = as.data.frame(cbind(cycle_expression[,"cycid"],0))
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
      if (temp_str_gene %in% cyc_gene_not_in_tgca) {
        print(paste0("keggname:",temp_str_gene))
        temp_str_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_str_gene),2][[1]]
      }
      temp_genesign = as.double(etSig[temp_str_gene,"sign_up_down"])*as.double(temp_str_num)
      temp_cyc_allgenenum_sum = temp_cyc_allgenenum_sum + temp_str_num
      #cyc_up_down[i,2] = as.double(cyc_up_down[i,2]) + temp_genesign #temp_genesign
      temp_cyc_genesign_sum = as.double(temp_cyc_genesign_sum) + temp_genesign #temp_genesign
      
      if (j == 1) {
        cyc_up_down[i,3] = as.character(etSig[temp_str_gene,"pvalue"])
      }else {
        cyc_up_down[i,3] = paste0(cyc_up_down[i,3],";",etSig[temp_str_gene,"pvalue"])  
      }
    }
    
    cyc_up_down[i,2] = temp_cyc_genesign_sum/temp_cyc_allgenenum_sum #各个基因的sign总和/基因出现的总次数
    
  }
  colnames(cyc_up_down) = c("cycid", "updownsign", "PValues")
  
  
  ## if pvalue>0.01 
  for (k in 1:length(rownames(cyc_up_down))) {
    is_pvalue_sd = 0
    pvalue_arr = unlist(strsplit(cyc_up_down[k,"PValues"],split = ";"))
    if(is.na(cyc_up_down[k,"PValues"])) {
      cyc_up_down[k, "is_pvalue_sd"] = 0
    } else {
      for (h in 1:length(pvalue_arr)) {
        if(pvalue_arr[h] == "NA" | pvalue_arr[h] == "NaN") { is_pvalue_sd = 0 } 
        else if(as.double(pvalue_arr[h]) < 0.001) { is_pvalue_sd = 1 }
      }
      cyc_up_down[k, "is_pvalue_sd"] = is_pvalue_sd
    }
  }
  cyc_up_down = cyc_up_down[,c("cycid", "updownsign", "is_pvalue_sd")]
  ################ 

  print(tumor_name)
  return(cyc_up_down)
}

#test
# gene_foldchange_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_foldchange_list.RData")
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_missing_list.RData")
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/cycle_expression.RData")

library(readr)
to_find_correct_symbol_v1_1 <- read_csv("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/to_find_correct_symbol_v1.1.csv", 
                                        col_names = FALSE)
cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]
#结合数据to_find_correct_symbol,得到新数据
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_updownsign = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  cyc_up_down = get_tumor_updownsign(tumors_array[i], gene_foldchange_list, cyc_gene_not_in_tgca_but_in_tofind, TRUE)
  tumor_cyc_updownsign[[tumors_array[i]]] = cyc_up_down
}





################################### 增加 structid ############################################                 
###############################################################################
# 环结构划分在res_cycle_struct_sort
# 结合环结构
#struct_cycle_sort_list 
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/2_cycle_struct_sort/result_struct/struct_cycle_sort_list.RData")


#继续向tumor_cyc_updownsign添加内容
merge_struct<-function(tumor_name, tumor_cyc_updownsign) {
  #struct
  temp_tumor_df = tumor_cyc_updownsign[[tumor_name]]
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["net"]]),"stcid"] = 1
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["hub"]]),"stcid"] = 2
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["idv"]]),"stcid"] = 3
  temp_tumor_df[which(temp_tumor_df$cycid %in% struct_cycle_sort_list[["nrl"]]),"stcid"] = 4
  #ifgap
  #temp_tumor_df[,"isgap"] = 'default'
  return(temp_tumor_df)
}


#merge
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_tumor_df = as.data.frame(merge_struct(tumors_array[i], tumor_cyc_updownsign))
  tumor_cyc_updownsign[[tumors_array[i]]] = temp_tumor_df
}

save(tumor_cyc_updownsign, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_updownsign.RData")











