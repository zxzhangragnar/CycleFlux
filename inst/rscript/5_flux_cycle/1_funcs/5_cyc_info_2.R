
##############################################################################################
##############################################################################################
# 
# 4列（5列） 
# 
# [1]          [2]          [3]             [4]               [5]             
# isgap        count        count         FC1(maxmin)       FC2(maxmin)
#            (Cancer)	     (Normal)
# 
# 
# ------------------------------------------------------------------------------------
#   [5]
# C的最大表达值所对应的FC取最小值
# 
# 
# (见手机照片： [4]5/2  [5]1/2)  
#
# 环 --->  反应 --->  g1:  10/20     (1/2)
#                     g2:  5/2       (5/2)
#
# [4]选择上面值最大的(10最大，虽然它最后计算出来的值(FC)1/2不是最大的)
# [5]选择最后计算出来的值(FC)最大的即5/2，虽然它上面值 5不是最大的 
# ------------------------------------------------------------------------------------
#   每个Rct中取最大gene值 （作为这个Rct的值）
# 每个环中取Rct值最小的那个
# 
# 
# ------------------------------------------------------------------------------------
#   
#   FC值计算
# 为了防止有0影响计算结果
# (Cancer+0.01)
# ——————————————
# (Normal+0.01)
# 
# 每个gene 的FoldChange 
# 
# Rct中的gene的FC最大的作为FC的值
#
#
#
#








######################################################################################################
# 将tumor_count\normal_count\tumor_FC添加到 tumor_cyc_enzyme_genecount_updown 中
# 见 find_gap_count : get_cyc_enzyme_genecount()





######################################################################################################
# 第一列:maxmin by count
# highest count of rct, smallest rct of cycid 
get_cyc_rid_maxmin_bycount <- function(tumor_name, tumor_cyc_enzyme_count_type) {
  # init
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
  
  temp_ec_count = tumor_cyc_enzyme_count_type[[tumor_name]]
  cycid_arr = unique(temp_ec_count$cycid)
  
  #result list
  temp_cycid_maxmin_bycount_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycid_maxmin_bycount_df) = c("cycid")
  
  for (i in 1:length(cycid_arr)) {
    cycid = cycid_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycid == cycid),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycid=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "gene_count"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "gene_count"] #5:gene_count
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max count
        rid_count_val_highest = rid_count_val_arr[1]
        #new
        rid_count_list[[rid_name]]["val"] = rid_count_val_highest
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$gene_count == rid_count_val_highest), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$gene_count == rid_count_val_highest), "gene_name"]           
      }
    }
    
    #2.in cyc: get min Rct
    ridnames = names(rid_count_list)
    rid_count_list_value = c()
    #new
    rid_count_list_ec = c()
    rid_count_list_gene = c()
    for (j in 1:length(rid_count_list)) {
      #new
      rid_count_list_value = append(rid_count_list_value,rid_count_list[[j]]["val"])
      rid_count_list_ec = append(rid_count_list_ec,rid_count_list[[j]]["ec"])
      rid_count_list_gene = append(rid_count_list_gene,rid_count_list[[j]]["gene"])
    }
    min_rct_idx = which(rid_count_list_value== min(rid_count_list_value), arr.ind = TRUE)[1]  
    
    #highest count of rct, smallest rct of cycid 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_count_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_count_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_count_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycid_maxmin_bycount_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycid_maxmin_bycount_df[i,"count"] = cyc_selected_rct_count_val
    #new
    temp_cycid_maxmin_bycount_df[i,"sel_ec"] = cyc_selected_rct_count_ec
    temp_cycid_maxmin_bycount_df[i,"sel_gene"] = cyc_selected_rct_count_gene
  }
  
  print(tumor_name)
  return(temp_cycid_maxmin_bycount_df)
}

#tumor_cyc_enzyme_count_type
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
cycid_maxmin_bycount_list = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_list = get_cyc_rid_maxmin_bycount(tumors_array[i], tumor_cyc_enzyme_count_type)
  cycid_maxmin_bycount_list[[tumors_array[i]]] = temp_list
}



######################################################################################################
# 第二列:maxmin by count_N
# highest count_N of rct, smallest rct of cycid 
get_cyc_rid_maxmin_bycount_N <- function(tumor_name, tumor_cyc_enzyme_count_type, normal_gene_count) {
  # init
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
  
  temp_ec_count = tumor_cyc_enzyme_count_type[[tumor_name]]
  cycid_arr = unique(temp_ec_count$cycid)
  
  #result list
  temp_cycid_maxmin_bycount_df_N = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycid_maxmin_bycount_df_N) = c("cycid")
  
  for (i in 1:length(cycid_arr)) {
    cycid = cycid_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycid == cycid),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycid=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "normal_gene_count"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "normal_gene_count"] #6:normal_count
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max FC
        rid_count_val_highest = rid_count_val_arr[which(rid_count_val_arr== max(rid_count_val_arr))[1]]
        #new
        rid_count_list[[rid_name]]["val"] = rid_count_val_highest
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$normal_gene_count == rid_count_val_highest), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$normal_gene_count == rid_count_val_highest), "gene_name"]    
      }

    }
    
    #2.in cyc: get min Rct
    ridnames = names(rid_count_list)
    rid_count_list_value = c()
    #new
    rid_count_list_ec = c()
    rid_count_list_gene = c()
    for (j in 1:length(rid_count_list)) {
      #new
      rid_count_list_value = append(rid_count_list_value,rid_count_list[[j]]["val"])
      rid_count_list_ec = append(rid_count_list_ec,rid_count_list[[j]]["ec"])
      rid_count_list_gene = append(rid_count_list_gene,rid_count_list[[j]]["gene"])
    }
    min_rct_idx = which(rid_count_list_value== min(rid_count_list_value), arr.ind = TRUE)[1]  
    
    #highest count of rct, smallest rct of cycid
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_count_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_count_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_count_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycid_maxmin_bycount_df_N[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycid_maxmin_bycount_df_N[i,"count"] = cyc_selected_rct_count_val
    #new
    temp_cycid_maxmin_bycount_df_N[i,"sel_ec"] = cyc_selected_rct_count_ec
    temp_cycid_maxmin_bycount_df_N[i,"sel_gene"] = cyc_selected_rct_count_gene
    
  }
  
  print(tumor_name)
  return(temp_cycid_maxmin_bycount_df_N)
}

#tumor_cyc_enzyme_count_type
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
#normal_gene_count
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_normal_gene_count.RData")
#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
cycid_maxmin_bycount_N_list = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_list = get_cyc_rid_maxmin_bycount_N(tumors_array[i], tumor_cyc_enzyme_count_type, normal_gene_count)
  cycid_maxmin_bycount_N_list[[tumors_array[i]]] = temp_list
}




######################################################################################################
# 第三列:maxmin by FC
# highest FC of rct, smallest rct of cycid 

get_cyc_rid_maxmin_byFC <- function(tumor_name, tumor_cyc_enzyme_count_type) {
  #tumor_name = "COAD"
  # init
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
  
  temp_ec_count = tumor_cyc_enzyme_count_type[[tumor_name]]
  cycid_arr = unique(temp_ec_count$cycid)
  
  #result list
  temp_cycid_maxmin_byFC_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycid_maxmin_byFC_df) = c("cycid")
  
  for (i in 1:length(cycid_arr)) {
    cycid = cycid_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycid == cycid),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycid=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "FoldChange"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "FoldChange"] #7:FoldChange
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max FC
        rid_count_val_highest = rid_count_val_arr[which(rid_count_val_arr== max(rid_count_val_arr))[1]]
        #new
        rid_count_list[[rid_name]]["val"] = rid_count_val_highest
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$FoldChange == rid_count_val_highest), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$FoldChange == rid_count_val_highest), "gene_name"]        
      }

    }
    
    #2.in cyc: get min Rct
    ridnames = names(rid_count_list)
    rid_count_list_value = c()
    #new
    rid_count_list_ec = c()
    rid_count_list_gene = c()
    for (j in 1:length(rid_count_list)) {
      #new
      rid_count_list_value = append(rid_count_list_value,rid_count_list[[j]]["val"])
      rid_count_list_ec = append(rid_count_list_ec,rid_count_list[[j]]["ec"])
      rid_count_list_gene = append(rid_count_list_gene,rid_count_list[[j]]["gene"])
    }
    min_rct_idx = which(rid_count_list_value== min(rid_count_list_value), arr.ind = TRUE)[1]  
    
    #highest count of rct, smallest rct of cycid 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_FC_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_FC_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_FC_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycid_maxmin_byFC_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycid_maxmin_byFC_df[i,"FC"] = cyc_selected_rct_FC_val
    #new
    temp_cycid_maxmin_byFC_df[i,"sel_ec"] = cyc_selected_rct_FC_ec
    temp_cycid_maxmin_byFC_df[i,"sel_gene"] = cyc_selected_rct_FC_gene
    
  }
  print(tumor_name)
  return(temp_cycid_maxmin_byFC_df)
}

#tumor_cyc_enzyme_count_type
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")

#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
cycid_maxmin_byFC_list = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_list = get_cyc_rid_maxmin_byFC(tumors_array[i], tumor_cyc_enzyme_count_type)
  cycid_maxmin_byFC_list[[tumors_array[i]]] = temp_list
}






######################################################################################################
# 第四列:maxmin by count_FC
# highest count gene's FC of rct, smallest rct of cycid 

get_cyc_rid_maxmin_by_countFC <- function(tumor_name, tumor_cyc_enzyme_count_type) {
  #tumor_name = "COAD"
  
  temp_ec_count = tumor_cyc_enzyme_count_type[[tumor_name]]
  cycid_arr = unique(temp_ec_count$cycid)
  
  #result list
  temp_cycid_maxmin_by_countFC_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycid_maxmin_by_countFC_df) = c("cycid")
  
  for (i in 1:length(cycid_arr)) {
    cycid = cycid_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycid == cycid),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycid=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "FoldChange"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        part_ec_count_rid = part_ec_count[which(part_ec_count$rid == rid_name),]
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "gene_count"]
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max FC
        rid_count_val_highest = rid_count_val_arr[which(rid_count_val_arr== max(rid_count_val_arr))[1]]
        rid_highest_count_sFC_val = part_ec_count_rid[which(part_ec_count_rid$gene_count == rid_count_val_highest),"FoldChange"]
        #new
        rid_count_list[[rid_name]]["val"] = rid_highest_count_sFC_val
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$FoldChange == rid_highest_count_sFC_val), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$FoldChange == rid_highest_count_sFC_val), "gene_name"]  
      }
    }
    
    #2.in cyc: get min Rct
    ridnames = names(rid_count_list)
    rid_count_list_value = c()
    #new
    rid_count_list_ec = c()
    rid_count_list_gene = c()
    for (j in 1:length(rid_count_list)) {
      #new
      rid_count_list_value = append(rid_count_list_value,rid_count_list[[j]]["val"])
      rid_count_list_ec = append(rid_count_list_ec,rid_count_list[[j]]["ec"])
      rid_count_list_gene = append(rid_count_list_gene,rid_count_list[[j]]["gene"])
    }
    min_rct_idx = which(rid_count_list_value== min(rid_count_list_value), arr.ind = TRUE)[1]  
    
    #highest count of rct, smallest rct of cycid 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_countFC_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_countFC_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_countFC_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycid_maxmin_by_countFC_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycid_maxmin_by_countFC_df[i,"countFC"] = cyc_selected_rct_countFC_val
    #new
    temp_cycid_maxmin_by_countFC_df[i,"sel_ec"] = cyc_selected_rct_countFC_ec
    temp_cycid_maxmin_by_countFC_df[i,"sel_gene"] = cyc_selected_rct_countFC_gene
    
  }
  print(tumor_name)
  return(temp_cycid_maxmin_by_countFC_df)
}


#tumor_cyc_enzyme_count_type
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")

#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
cycid_maxmin_by_countFC_list = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_list = get_cyc_rid_maxmin_by_countFC(tumors_array[i], tumor_cyc_enzyme_count_type)
  cycid_maxmin_by_countFC_list[[tumors_array[i]]] = temp_list
}



# 得到了
# cycid_maxmin_bycount_list
# cycid_maxmin_bycount_N_list
# cycid_maxmin_byFC_list
# cycid_maxmin_by_countFC_list



######################################################################################################
#基于 tumor_cyc_enzyme_count_type 给 tumor_cyc_updownsign 增加4列
# 第一列:count
# 第二列:count_N
# 第三列:FC
# 第四列:count_FC

# tumor_cyc_updownsign
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_updownsign.RData")



#tool function
merge_temp_col_n<-function(temp_list,tumor_name,type_name) {
  tempcol = temp_list[[tumor_name]]
  tempcol[,"merge"] = paste0(tempcol[,"sel_rct"],"-",tempcol[,"sel_ec"],"-",tempcol[,"sel_gene"],": ",tempcol[,3])
  tempcol = as.data.frame(tempcol[,c("cycid", "merge")])
  colnames(tempcol) = c("cycid", type_name)
  return(tempcol)
}

# new list (merged)
tumor_cyc_updown_info = list()
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_tumor_df = tumor_cyc_updownsign[[tumors_array[i]]]
  col_1 = merge_temp_col_n(cycid_maxmin_bycount_list, tumors_array[i], "count")
  col_2 = merge_temp_col_n(cycid_maxmin_bycount_N_list, tumors_array[i], "count_N")
  col_3 = merge_temp_col_n(cycid_maxmin_byFC_list, tumors_array[i], "FC")
  col_4 = merge_temp_col_n(cycid_maxmin_by_countFC_list, tumors_array[i], "countFC")
  
  have_cycid = temp_tumor_df[,"cycid"]
  new_tumor_df = as.data.frame(col_1[which(col_1$cycid %in% have_cycid),2])
  new_tumor_df = cbind(new_tumor_df, col_2[which(col_2$cycid %in% have_cycid),2])
  new_tumor_df = cbind(new_tumor_df, col_3[which(col_3$cycid %in% have_cycid),2])
  new_tumor_df = cbind(new_tumor_df, col_4[which(col_4$cycid %in% have_cycid),2])
  colnames(new_tumor_df) = c("count", "count_N", "FC", "countFC")
  
  temp_tumor_df = cbind(temp_tumor_df, new_tumor_df)
  
  tumor_cyc_updown_info[[tumors_array[i]]] = temp_tumor_df
}


#得到结果tumor_cyc_updown_info




###############################  9_cycle_flux_list  ##########################################
################################################################################################
get_cycle_flux_list<-function(tumor_name, tumor_cyc_updown_info) {
  part1_df = tumor_cyc_updown_info[[tumor_name]][,c("cycid", "stcid","is_pvalue_sd")]
  colnames(part1_df) = c("cycle_id", "struct_id","if_pvalue")
  
  # part2_df = as.data.frame(tumor_cyc_updown_info[[tumor_name]][,c("updownsign")])
  # colnames(part2_df) = c("updown_num")
  
  part3_df = tumor_cyc_updown_info[[tumor_name]][,c("count","count_N","FC","countFC")]
  colnames(part3_df) = c("count","count_N","foldchange","count_foldchange")
  
  # cycle_updown_df = cbind(part1_df, part2_df, part3_df)
  cycle_updown_df = cbind(part1_df, part3_df)
  return(cycle_updown_df)
}

#修改tumor_cyc_updown_info列的顺序

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
cycle_flux_list = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  cycle_updown_df = get_cycle_flux_list(tumors_array[i], tumor_cyc_updown_info)
  cycle_flux_list[[tumors_array[i]]] = cycle_updown_df
}


save(cycle_flux_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_final/cycle_flux_list.RData")

################################################################################################
################################################################################################


