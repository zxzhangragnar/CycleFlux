########################################### 1 #########################################################
# 目的:每个环的每个反应Reaction为独立的一行
# 添加1列foldchangesign
#######################################################################################################


########################################### 2 ########################################################
# 目标：为tumor_cyc_gap_foldchangsign.RData添加8列
# "gene_fold_change"
# "gene_pvalue"
# "t_gene_mean_val"
# "n_gene_mean_val"
# "t_gene_0_pct"
# "n_gene_0_pct"
# "t_gene_log_lss1_pct"
# "n_gene_log_lss1_pct"
#######################################################################################################


########################################### 3 ########################################################
# 目标：为cycle_edge_flux_list
# 更新原来的cycle_edge_flux_list的ifgap项的内容
# 去掉4列：
# "ifgap"
# "gene_num"
# "cyc_len"
# "gaprate"
# 增加2列：
# "ifgap"
# "ifup"
#
#######################################################################################################





#################################### 1 #############################################
####################################################################################

# 计算环中某个反应rct的sign 上调下调值
get_tumor_gap_foldchangesign<-function(gene_foldchange_list, tumor_name, cycle_edge_expression, gene_missing_list, is_filter)
{
  cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]
  
  #print(tumor_name)
  stat_r = as.data.frame(gene_foldchange_list[tumor_name])
  names(stat_r) = c("gene_id","log2FoldChange","up_down","sign_up_down","pvalue")
  et = stat_r
  if(is_filter) {
    #etSig = et[which(et$pvalue < 0.05 & abs(et$log2FoldChange) > 1),]
    etSig = et
  }else {
    etSig = et
  }
  cyc_up_down = cycle_edge_expression
  cyc_up_down[,"foldchangesign"] = 0
  
  for (i in 1:length(cycle_edge_expression[,1])) {
    temp_row_fcsign_sum = 0 #当前行的各个gene的fcsign之和
    gene_arr = cycle_edge_expression[i,"gene_symbol"]
    gene_arr <- unlist(strsplit(gene_arr,split = ";"))
    for (j in 1:length(gene_arr)) {
      temp_gene_str = gene_arr[j]
      if (temp_gene_str %in% cyc_gene_not_in_tgca_but_in_tofind) {
        temp_gene_str = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_gene_str),2][[1]]
      }
      genesign = as.double(etSig[temp_gene_str,"log2FoldChange"])
      temp_row_fcsign_sum = temp_row_fcsign_sum + genesign
    }
    cyc_up_down[i,"foldchangesign"] = temp_row_fcsign_sum/length(gene_arr) #更公平,并不是gene多的边fc的值就更大
  }
  
  return(cyc_up_down)
}









#################################### 2 #############################################
####################################################################################

add_tumor_cyc_gap_foldchangsign<-function(tumor_name, gene_missing_list, to_find_correct_symbol_v1_1, tumor_cyc_gap_foldchangsign, gene_foldchange_list, gene_meanval_list)
{
  cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]
  
  temp_tumor_cyc_gap_foldchangsign = tumor_cyc_gap_foldchangsign[[tumor_name]]
  temp_gene_foldchange_list = gene_foldchange_list[[tumor_name]]
  temp_gene_meanval_list = gene_meanval_list[[tumor_name]]
  
  temp_tumor_cyc_gap_foldchangsign_p1 = temp_tumor_cyc_gap_foldchangsign
  
  
  #set new columns total:8 
  temp_tumor_cyc_gap_foldchangsign_p1[, "gene_fold_change"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "gene_pvalue"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "t_gene_mean_val"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "n_gene_mean_val"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "t_gene_0_pct"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "n_gene_0_pct"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "t_gene_log_lss1_pct"] = NA
  temp_tumor_cyc_gap_foldchangsign_p1[, "n_gene_log_lss1_pct"] = NA
  
  # 1.gene_fold_change 2.gene_pvalue
  
  for (i in 1:length(temp_tumor_cyc_gap_foldchangsign_p1[,1])) {
    tmp_gene_symbol = temp_tumor_cyc_gap_foldchangsign_p1[i,"gene_symbol"]
    tmp_gene_symbol_arr = unlist(strsplit(tmp_gene_symbol,split = ";"))
    
    ##
    gene_fc_arr = c()
    gene_pvalue_arr = c()
    
    t_gene_mean_val_arr = c()
    n_gene_mean_val_arr = c()
    t_gene_0_pct_arr = c()
    n_gene_0_pct_arr = c()
    t_gene_log_lss1_pct_arr = c()
    n_gene_log_lss1_pct_arr = c()
    
    ##
    for (j in 1:length(tmp_gene_symbol_arr)) {
      tmp_gene = tmp_gene_symbol_arr[j]
      
      #找不到的gene从cyc_gene_not_in_tgca_but_in_tofind中找
      if (tmp_gene %in% cyc_gene_not_in_tgca_but_in_tofind) {
        tmp_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == tmp_gene),2][[1]]
      }
      
      ## gene_foldchange_list
      tmp_fc = temp_gene_foldchange_list[tmp_gene,"foldchange"]
      tmp_pvalue = temp_gene_foldchange_list[tmp_gene,"pvalue"]
      ## gene_meanval_list
      temp_t_gene_mean_val = temp_gene_meanval_list[tmp_gene, "t_gene_mean_val"]
      temp_n_gene_mean_val = temp_gene_meanval_list[tmp_gene, "n_gene_mean_val"]
      temp_t_gene_0_pct = temp_gene_meanval_list[tmp_gene, "t_gene_0_pct"]
      temp_n_gene_0_pct = temp_gene_meanval_list[tmp_gene, "n_gene_0_pct"]
      temp_t_gene_log_lss1_pct = temp_gene_meanval_list[tmp_gene, "t_gene_log_lss1_pct"]
      temp_n_gene_log_lss1_pct = temp_gene_meanval_list[tmp_gene, "n_gene_log_lss1_pct"]
      
      #put
      gene_fc_arr = append(gene_fc_arr, paste0(tmp_gene,":",tmp_fc))
      gene_pvalue_arr = append(gene_pvalue_arr, paste0(tmp_gene,":",tmp_pvalue))
      
      t_gene_mean_val_arr = append(t_gene_mean_val_arr, paste0(tmp_gene,":",temp_t_gene_mean_val))
      n_gene_mean_val_arr = append(n_gene_mean_val_arr, paste0(tmp_gene,":",temp_n_gene_mean_val))
      t_gene_0_pct_arr = append(t_gene_0_pct_arr, paste0(tmp_gene,":",temp_t_gene_0_pct))
      n_gene_0_pct_arr = append(n_gene_0_pct_arr, paste0(tmp_gene,":",temp_n_gene_0_pct))
      t_gene_log_lss1_pct_arr = append(t_gene_log_lss1_pct_arr, paste0(tmp_gene,":",temp_t_gene_log_lss1_pct))
      n_gene_log_lss1_pct_arr = append(n_gene_log_lss1_pct_arr, paste0(tmp_gene,":",temp_n_gene_log_lss1_pct))
    }
    
    #put in
    gene_fc_str = paste(gene_fc_arr, collapse  = ";")
    gene_pvalue_str = paste(gene_pvalue_arr, collapse  = ";")
    
    t_gene_mean_val_str = paste(t_gene_mean_val_arr, collapse  = ";")
    n_gene_mean_val_str = paste(n_gene_mean_val_arr, collapse  = ";")
    t_gene_0_pct_str = paste(t_gene_0_pct_arr, collapse  = ";")
    n_gene_0_pct_str = paste(n_gene_0_pct_arr, collapse  = ";")
    t_gene_log_lss1_pct_str = paste(t_gene_log_lss1_pct_arr, collapse  = ";")
    n_gene_log_lss1_pct_str = paste(n_gene_log_lss1_pct_arr, collapse  = ";")
    
    # 1.gene_fold_change
    temp_tumor_cyc_gap_foldchangsign_p1[i, "gene_fold_change"] = gene_fc_str
    # 2.gene_pvalue
    temp_tumor_cyc_gap_foldchangsign_p1[i, "gene_pvalue"] = gene_pvalue_str
    
    # 3.t_gene_mean_val
    temp_tumor_cyc_gap_foldchangsign_p1[i, "t_gene_mean_val"] = t_gene_mean_val_str
    # 4.n_gene_mean_val
    temp_tumor_cyc_gap_foldchangsign_p1[i, "n_gene_mean_val"] = n_gene_mean_val_str
    # 5.t_gene_0_pct
    temp_tumor_cyc_gap_foldchangsign_p1[i, "t_gene_0_pct"] = t_gene_0_pct_str
    # 6.n_gene_0_pct
    temp_tumor_cyc_gap_foldchangsign_p1[i, "n_gene_0_pct"] = n_gene_0_pct_str
    # 7.t_gene_log_lss1_pct
    temp_tumor_cyc_gap_foldchangsign_p1[i, "t_gene_log_lss1_pct"] = t_gene_log_lss1_pct_str
    # 8.n_gene_log_lss1_pct
    temp_tumor_cyc_gap_foldchangsign_p1[i, "n_gene_log_lss1_pct"] = n_gene_log_lss1_pct_str
  }
  temp_tumor_cyc_gap_foldchangsign = temp_tumor_cyc_gap_foldchangsign_p1
  
  return(temp_tumor_cyc_gap_foldchangsign)
}






#################################### 3 #############################################
####################################################################################
get_gap_c_max_gene_foldchange<-function(tumor_name, cycle_edge_flux_list) {
  
  aaa<-cycle_edge_flux_list[[tumor_name]][,"gene_fold_change"] # gene_fold_change
  ddd<-c()
  for(i in 1:length(aaa))
  {
    bbb<-unlist(strsplit(aaa[i],";"))
    ccc<-c()
    ccc2<-c()
    for(j in 1:length(bbb))
    {
      ccc<-c(ccc,unlist(strsplit(bbb[j],":"))[1])
      ccc2<-c(ccc2,unlist(strsplit(bbb[j],":"))[2])
    }
    ccc2<-as.numeric(ccc2)
    names(ccc2)<-ccc
    ddd<-c(ddd,max(ccc2,na.rm=T)) #添加ccc2中的最大者
  }
  
  return(ddd)
}

get_gap_c_max_gene_meanval<-function(tumor_name, cycle_edge_flux_list) {
  
  ddd2<-c()
  aaa<-cycle_edge_flux_list[[tumor_name]][,"t_gene_mean_val"] # t_gene_mean_val
  for(i in 1:length(aaa))
  {
    bbb<-unlist(strsplit(aaa[i],";"))
    ccc<-c()
    ccc2<-c()
    for(j in 1:length(bbb))
    {
      ccc<-c(ccc,unlist(strsplit(bbb[j],":"))[1])
      ccc2<-c(ccc2,unlist(strsplit(bbb[j],":"))[2])
    }
    ccc2<-as.numeric(ccc2)
    names(ccc2)<-ccc
    ddd2<-c(ddd2,max(ccc2,na.rm=T))
  }
  
  return(ddd2)
}


get_gap_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  c_max_gene_foldchange = get_gap_c_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  c_max_gene_meanval = get_gap_c_max_gene_meanval(tumor_name, cycle_edge_flux_list)
  
  ddd = c_max_gene_foldchange
  ddd2 = c_max_gene_meanval
  
  gap_c<-unique(cycle_edge_flux_list[[tumor_name]][which((ddd<(-2))|(ddd2<1.71)),1])
  return(gap_c)
}



get_up_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  c_max_gene_foldchange = get_gap_c_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  c_max_gene_meanval = get_gap_c_max_gene_meanval(tumor_name, cycle_edge_flux_list)
  
  ddd = c_max_gene_foldchange
  ddd2 = c_max_gene_meanval
  
  up_c<-unique(cycle_edge_flux_list[[tumor_name]][which((ddd>(1))&(ddd2>10)),1])#p.value of diff<0.001, log(fc)>0.5, and the gene show significant up regulation are larger than 10
  return(up_c)
}


#################################### 4 #############################################
####################################################################################
## main function
get_cycle_edge_flux_list <- function(cycle_edge_expression, tumors_array, gene_foldchange_list, gene_meanval_list, gene_missing_list, to_find_correct_symbol_v1_1, par1, par2, par3, par4) {
  #tumor_cyc_gap_foldchangsign = get_tumor_cyc_gap_foldchangsign(cycle_edge_expression)
  
  #tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
  tumor_cyc_gap_foldchangsign = list()
  for (i in 1:length(tumors_array)) {
    #保存到tumor字典
    gap_updown = get_tumor_gap_foldchangesign(gene_foldchange_list, tumors_array[i], cycle_edge_expression, gene_missing_list, FALSE)
    tumor_cyc_gap_foldchangsign[[tumors_array[i]]] = gap_updown
  }
  
  cycle_edge_flux_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    print(tumor_name)
    statr = as.data.frame(add_tumor_cyc_gap_foldchangsign(tumor_name, gene_missing_list, to_find_correct_symbol_v1_1, tumor_cyc_gap_foldchangsign, gene_foldchange_list, gene_meanval_list))
    cycle_edge_flux_list[[tumor_name]] = statr
    
    #######################################################################
    ## gap up
    c_max_gene_foldchange = get_gap_c_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
    c_max_gene_meanval = get_gap_c_max_gene_meanval(tumor_name, cycle_edge_flux_list)
    
    ddd = c_max_gene_foldchange
    ddd2 = c_max_gene_meanval
    
    #ifgap
    cycle_edge_flux_list[[tumor_name]][, "ifgap"] = "normal"
    cycle_edge_flux_list[[tumor_name]][, "ifup"] = "normal"
    
    # cycle_edge_flux_list[[tumor_name]][which((ddd>(1))&(ddd2>10)), "ifup"] = "up"
    # cycle_edge_flux_list[[tumor_name]][which((ddd<(-2))|(ddd2<1.71)), "ifgap"] = "gap"
    cycle_edge_flux_list[[tumor_name]][which((ddd>(par1))&(ddd2>par2)), "ifup"] = "up"
    cycle_edge_flux_list[[tumor_name]][which((ddd<(par3))|(ddd2<par4)), "ifgap"] = "gap"
    #######################################################################
    
  }
  return(cycle_edge_flux_list)
}







#####################################################################################
# cycle_edgesucs_expression_in cycle_edgesucs_expression_out
cycle_edge_flux_list_inout_main <- function(output_path, res_path, package_path, input_tumor_name, par1=1, par2=10, par3=-2, par4=1.71) {
  #init
  load(file.path(output_path, "cycle_edgesucs_expression_in.RData"))
  load(file.path(output_path, "cycle_edgesucs_expression_out.RData"))
  
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_foldchange_list.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_meanval_list.RData"))
  load(file.path(res_path, "3_flux_subnet/result_tool/gene_missing_list.RData"))
  
  library(readr)
  to_find_correct_symbol_v1_1 = read_csv(file.path(package_path, "tool_data/to_find_correct_symbol_v1.1.csv"), col_names = FALSE)
  tumors_array = c(input_tumor_name)
  
  ##
  cycle_edge_flux_list_in = get_cycle_edge_flux_list(cycle_edgesucs_expression_in, tumors_array, gene_foldchange_list, gene_meanval_list, gene_missing_list, to_find_correct_symbol_v1_1, par1, par2, par3, par4)
  cycle_edge_flux_list_out = get_cycle_edge_flux_list(cycle_edgesucs_expression_out, tumors_array, gene_foldchange_list, gene_meanval_list, gene_missing_list, to_find_correct_symbol_v1_1, par1, par2, par3, par4)
  
  #save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  
  res_file_path_in = file.path(res_path, res_sub_path, "cycle_edge_flux_list_in.RData")
  res_file_path_out = file.path(res_path, res_sub_path, "cycle_edge_flux_list_out.RData")
  
  # save(cycle_edge_flux_list_in, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_in.RData")
  # save(cycle_edge_flux_list_out, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_out.RData")
  
  save(cycle_edge_flux_list_in, file=res_file_path_in)
  save(cycle_edge_flux_list_out, file=res_file_path_out)
}





# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
# 
# cycle_edge_flux_list_inout_main(output_path, res_path, package_path, input_tumor_name)






