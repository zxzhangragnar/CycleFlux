
##############################################################################################
##############################################################################################
# [1]          [2]          [3]             [4]               [5]             
# isgap        count        count         FC1(maxmin)       FC2(maxmin)
#            (Cancer)	     (Normal)
##############################################################################################

######################################################################################################
# 第一列:maxmin by count
# highest count of rct, smallest rct of cycle_id 
get_cyc_rid_maxmin_bycount <- function(tumor_name, cycle_enzyme_stat) {

  temp_ec_count = cycle_enzyme_stat[[tumor_name]]
  cycle_id_arr = unique(temp_ec_count$cycle_id)
  
  #result list
  temp_cycle_id_maxmin_bycount_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycle_id_maxmin_bycount_df) = c("cycle_id")
  
  for (i in 1:length(cycle_id_arr)) {
    cycle_id = cycle_id_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycle_id == cycle_id),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycle_id=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "tumor_gene_count"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "tumor_gene_count"] #5:gene_count
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max count
        rid_count_val_highest = rid_count_val_arr[1]
        #new
        rid_count_list[[rid_name]]["val"] = rid_count_val_highest
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$tumor_gene_count == rid_count_val_highest), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$tumor_gene_count == rid_count_val_highest), "gene_name"]           
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
    
    #highest count of rct, smallest rct of cycle_id 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_count_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_count_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_count_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycle_id_maxmin_bycount_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycle_id_maxmin_bycount_df[i,"count"] = cyc_selected_rct_count_val
    #new
    temp_cycle_id_maxmin_bycount_df[i,"sel_ec"] = cyc_selected_rct_count_ec
    temp_cycle_id_maxmin_bycount_df[i,"sel_gene"] = cyc_selected_rct_count_gene
  }
  
  return(temp_cycle_id_maxmin_bycount_df)
}





######################################################################################################
# 第二列:maxmin by count_N
# highest count_N of rct, smallest rct of cycle_id 
get_cyc_rid_maxmin_bycount_N <- function(tumor_name, cycle_enzyme_stat, normal_gene_count) {

  temp_ec_count = cycle_enzyme_stat[[tumor_name]]
  cycle_id_arr = unique(temp_ec_count$cycle_id)
  
  #result list
  temp_cycle_id_maxmin_bycount_df_N = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycle_id_maxmin_bycount_df_N) = c("cycle_id")
  
  for (i in 1:length(cycle_id_arr)) {
    cycle_id = cycle_id_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycle_id == cycle_id),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycle_id=0
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
    
    #highest count of rct, smallest rct of cycle_id
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_count_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_count_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_count_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycle_id_maxmin_bycount_df_N[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycle_id_maxmin_bycount_df_N[i,"count"] = cyc_selected_rct_count_val
    #new
    temp_cycle_id_maxmin_bycount_df_N[i,"sel_ec"] = cyc_selected_rct_count_ec
    temp_cycle_id_maxmin_bycount_df_N[i,"sel_gene"] = cyc_selected_rct_count_gene
    
  }
  
  return(temp_cycle_id_maxmin_bycount_df_N)
}





######################################################################################################
# 第三列:maxmin by FC
# highest FC of rct, smallest rct of cycle_id 

get_cyc_rid_maxmin_byFC <- function(tumor_name, cycle_enzyme_stat) {

  temp_ec_count = cycle_enzyme_stat[[tumor_name]]
  cycle_id_arr = unique(temp_ec_count$cycle_id)
  
  #result list
  temp_cycle_id_maxmin_byFC_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycle_id_maxmin_byFC_df) = c("cycle_id")
  
  for (i in 1:length(cycle_id_arr)) {
    cycle_id = cycle_id_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycle_id == cycle_id),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycle_id=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "foldchange"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "foldchange"] #7:foldchange
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max FC
        rid_count_val_highest = rid_count_val_arr[which(rid_count_val_arr== max(rid_count_val_arr))[1]]
        #new
        rid_count_list[[rid_name]]["val"] = rid_count_val_highest
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$foldchange == rid_count_val_highest), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$foldchange == rid_count_val_highest), "gene_name"]        
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
    
    #highest count of rct, smallest rct of cycle_id 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_FC_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_FC_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_FC_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycle_id_maxmin_byFC_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycle_id_maxmin_byFC_df[i,"FC"] = cyc_selected_rct_FC_val
    #new
    temp_cycle_id_maxmin_byFC_df[i,"sel_ec"] = cyc_selected_rct_FC_ec
    temp_cycle_id_maxmin_byFC_df[i,"sel_gene"] = cyc_selected_rct_FC_gene
  }
  
  return(temp_cycle_id_maxmin_byFC_df)
}







######################################################################################################
# 第四列:maxmin by count_FC
# highest count gene's FC of rct, smallest rct of cycle_id 

get_cyc_rid_maxmin_by_countFC <- function(tumor_name, cycle_enzyme_stat) {
  
  temp_ec_count = cycle_enzyme_stat[[tumor_name]]
  cycle_id_arr = unique(temp_ec_count$cycle_id)
  #result list
  temp_cycle_id_maxmin_by_countFC_df = as.data.frame(unique(temp_ec_count[,1]))
  names(temp_cycle_id_maxmin_by_countFC_df) = c("cycle_id")
  
  for (i in 1:length(cycle_id_arr)) {
    cycle_id = cycle_id_arr[i]
    part_ec_count = temp_ec_count[which(temp_ec_count$cycle_id == cycle_id),]
    rid_types = unique(part_ec_count$rid)
    
    # rid_count for cycle_id=0
    rid_count_list = list()
    for (r in 1:length(rid_types)) {
      rid_name = rid_types[r]
      if(is.na(part_ec_count[which(part_ec_count$rid == rid_name), "foldchange"])){
        rid_count_list[[rid_name]]["val"] = NA
        rid_count_list[[rid_name]]["ec"] = NA
        rid_count_list[[rid_name]]["gene"] = NA
      }else {
        part_ec_count_rid = part_ec_count[which(part_ec_count$rid == rid_name),]
        rid_count_val_arr = part_ec_count[which(part_ec_count$rid == rid_name), "tumor_gene_count"]
        rid_count_val_arr <- rid_count_val_arr[!is.na(rid_count_val_arr)] # remove NA
        # 1.in Rct: get max FC
        rid_count_val_highest = rid_count_val_arr[which(rid_count_val_arr== max(rid_count_val_arr))[1]]
        rid_highest_count_sFC_val = part_ec_count_rid[which(part_ec_count_rid$tumor_gene_count == rid_count_val_highest),"foldchange"]
        #new
        rid_count_list[[rid_name]]["val"] = rid_highest_count_sFC_val
        rid_count_list[[rid_name]]["ec"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$foldchange == rid_highest_count_sFC_val), "enzyme"]
        rid_count_list[[rid_name]]["gene"] = part_ec_count[which(part_ec_count$rid == rid_name & part_ec_count$foldchange == rid_highest_count_sFC_val), "gene_name"]  
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
    
    #highest count of rct, smallest rct of cycle_id 
    cyc_selected_rct_name = ridnames[min_rct_idx]
    cyc_selected_rct_countFC_val = rid_count_list_value[min_rct_idx]
    #new
    cyc_selected_rct_countFC_ec = rid_count_list_ec[min_rct_idx]
    cyc_selected_rct_countFC_gene = rid_count_list_gene[min_rct_idx]
    
    temp_cycle_id_maxmin_by_countFC_df[i,"sel_rct"] = cyc_selected_rct_name
    temp_cycle_id_maxmin_by_countFC_df[i,"countFC"] = cyc_selected_rct_countFC_val
    #new
    temp_cycle_id_maxmin_by_countFC_df[i,"sel_ec"] = cyc_selected_rct_countFC_ec
    temp_cycle_id_maxmin_by_countFC_df[i,"sel_gene"] = cyc_selected_rct_countFC_gene
    
  }
  return(temp_cycle_id_maxmin_by_countFC_df)
}

# cycle_id_maxmin_bycount_list
# cycle_id_maxmin_bycount_N_list
# cycle_id_maxmin_byFC_list
# cycle_id_maxmin_by_countFC_list


######################################################################################################
# cycle_enzyme_stat 给 cycle_enzyme_stat 增加4列
# 第一列:count
# 第二列:count_N
# 第三列:FC
# 第四列:count_FC
merge_temp_col_n<-function(temp_list,tumor_name,type_name) {
  tempcol = temp_list[[tumor_name]]
  tempcol[,"merge"] = paste0(tempcol[,"sel_rct"],"-",tempcol[,"sel_ec"],"-",tempcol[,"sel_gene"],": ",tempcol[,3])
  tempcol = as.data.frame(tempcol[,c("cycle_id", "merge")])
  colnames(tempcol) = c("cycle_id", type_name)
  return(tempcol)
}


################################################################################################

cyc_info_2_main <- function(res_path, input_tumor_name) {
  
  load(file.path(res_path, "5_flux_cycle/result_enzyme/cycle_enzyme_stat.RData"))
  load(file.path(res_path, "5_flux_cycle/result_enzyme/cycle_merge_stat.RData"))
  
  #test
  tumors_array = c(input_tumor_name)
  cycle_id_maxmin_bycount_list = list()
  cycle_id_maxmin_bycount_N_list = list()
  cycle_id_maxmin_byFC_list = list()
  cycle_id_maxmin_by_countFC_list = list()
  cycle_flux_list = list()
  for (i in 1:length(tumors_array)) {
    temp_list = get_cyc_rid_maxmin_bycount(tumors_array[i], cycle_enzyme_stat)
    cycle_id_maxmin_bycount_list[[tumors_array[i]]] = temp_list
    
    temp_list = get_cyc_rid_maxmin_bycount_N(tumors_array[i], cycle_enzyme_stat, normal_gene_count)
    cycle_id_maxmin_bycount_N_list[[tumors_array[i]]] = temp_list
    
    temp_list = get_cyc_rid_maxmin_byFC(tumors_array[i], cycle_enzyme_stat)
    cycle_id_maxmin_byFC_list[[tumors_array[i]]] = temp_list
    
    temp_list = get_cyc_rid_maxmin_by_countFC(tumors_array[i], cycle_enzyme_stat)
    cycle_id_maxmin_by_countFC_list[[tumors_array[i]]] = temp_list
    
    ##
    temp_tumor_df = cycle_merge_stat[[tumors_array[i]]]
    col_1 = merge_temp_col_n(cycle_id_maxmin_bycount_list, tumors_array[i], "count")
    col_2 = merge_temp_col_n(cycle_id_maxmin_bycount_N_list, tumors_array[i], "count_N")
    col_3 = merge_temp_col_n(cycle_id_maxmin_byFC_list, tumors_array[i], "FC")
    col_4 = merge_temp_col_n(cycle_id_maxmin_by_countFC_list, tumors_array[i], "countFC")
    
    have_cycle_id = temp_tumor_df[,"cycle_id"]
    new_tumor_df = as.data.frame(col_1[which(col_1$cycle_id %in% have_cycle_id),2])
    new_tumor_df = cbind(new_tumor_df, col_2[which(col_2$cycle_id %in% have_cycle_id),2])
    new_tumor_df = cbind(new_tumor_df, col_3[which(col_3$cycle_id %in% have_cycle_id),2])
    new_tumor_df = cbind(new_tumor_df, col_4[which(col_4$cycle_id %in% have_cycle_id),2])
    colnames(new_tumor_df) = c("count", "count_N", "FC", "countFC")
    
    cycle_flux_list[[tumors_array[i]]] = cbind(temp_tumor_df, new_tumor_df)
  }
  
  #save
  res_sub_path = "5_flux_cycle/result_final"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_flux_list.RData")
  save(cycle_flux_list, file=res_file_path)
}



# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
# 
# cyc_info_2_main(res_path, input_tumor_name)


