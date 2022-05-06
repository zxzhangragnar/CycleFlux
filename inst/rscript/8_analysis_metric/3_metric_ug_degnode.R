##################################################################################
################################################################################## 
# Task:
#             1.环内 （环上某node --> 环外某node --> 环上某node）
# 			      2.出环 （环上某node --> 环外某node ）
# 			      3.进环 （环外某node --> 环上某node ）
#
#
#     往后走一步的Reaction     #这些Reaction是否也在别的环里 1在 0不在
#
##################################################################################






#########################################################################
# 1.统计
# 全部环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# 全部环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#########################################################################

#########################################################################
# 2.统计
# ug环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# ug环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#########################################################################

#########################################################################
# 3.统计(按流方向统计)
# 全部环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# 全部环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#
# od<---c_out<---c_in<---ind
# od<---c_out<---c_in<---ind
#########################################################################

#########################################################################
# 4.统计(按流方向统计)
# ug环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# ug环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#
# od<---c_out<---c_in<---ind
# od<---c_out<---c_in<---ind
#########################################################################

# 5.在suc1_code中做 metric 统计 od 和 ind 的化合物频率(gap化合物频率 up化合物频率)
# 在suc1_code中建立新的文件夹









##################################################################################################################################################
## tool function
compare_with_all_Freq <- function(table_temp) {
  #对比总的出现次数
  ## 化合物在全部度中的频率
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/table_all_deg_cpds.RData")
  table_all_deg_cpds = table_all_deg_cpds
  for (j in 1:length(table_temp[,1])) {
    tmp_comp = as.character(table_temp[j, "compound"])
    tmp_Freq = table_temp[j, "Freq"]
    all_Freq = table_all_deg_cpds[which(table_all_deg_cpds$compound == tmp_comp), "all_Freq"]
    Freq_ratio = tmp_Freq/all_Freq
    table_temp[j, "all_Freq"] = all_Freq
    table_temp[j, "Freq_ratio"] = Freq_ratio
  }
  table_temp = table_temp[order(table_temp$Freq_ratio,table_temp$Freq, -table_temp$all_Freq, decreasing = TRUE),]
  
  return(table_temp)
}

## function
get_table_up_deg_cpds <- function(deg_inout, temp_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) {
  ## 1.up
  up_rows = temp_pct[which(temp_pct$ifup == "up"), ]
  up_ind_cpds = c()
  up_od_cpds = c()
  for (j in 1:length(up_rows)) {
    tmp_cid = up_rows[j, "cycid"]
    tmp_cin = up_rows[j, "c_in"]
    tmp_cout = up_rows[j, "c_out"]
    #1.1 up_ind_cpds
    tmp_cin_ind = cycle_edgesucs_in[which((cycle_edgesucs_in$cycid == tmp_cid) & (cycle_edgesucs_in$node == tmp_cin)), "ind"]
    tmp_cout_ind = cycle_edgesucs_in[which((cycle_edgesucs_in$cycid == tmp_cid) & (cycle_edgesucs_in$node == tmp_cout)), "ind"]
    up_ind_cpds = append(up_ind_cpds, tmp_cin_ind)
    up_ind_cpds = append(up_ind_cpds, tmp_cout_ind)
    #1.2 up_od_cpds
    tmp_cin_od = cycle_edgesucs_out[which((cycle_edgesucs_out$cycid == tmp_cid) & (cycle_edgesucs_out$node == tmp_cin)), "od"]
    tmp_cout_od = cycle_edgesucs_out[which((cycle_edgesucs_out$cycid == tmp_cid) & (cycle_edgesucs_out$node == tmp_cout)), "od"]
    up_od_cpds = append(up_od_cpds, tmp_cin_od)
    up_od_cpds = append(up_od_cpds, tmp_cout_od)
  }
  
  table_up_ind_cpds = as.data.frame(table(up_ind_cpds))
  colnames(table_up_ind_cpds) = c("compound", "Freq")

  table_up_od_cpds = as.data.frame(table(up_od_cpds))
  colnames(table_up_od_cpds) = c("compound", "Freq")
  
  #与这些化合物总的出现次数对比
  table_up_ind_cpds = compare_with_all_Freq(table_up_ind_cpds)
  table_up_od_cpds = compare_with_all_Freq(table_up_od_cpds)
  
  #删去never_considered_comp 中的化合物
  table_up_ind_cpds = table_up_ind_cpds[-which(table_up_ind_cpds$compound %in% never_considered_comp_names),]
  table_up_od_cpds = table_up_od_cpds[-which(table_up_od_cpds$compound %in% never_considered_comp_names),]
  
  
  if(deg_inout == "ind") {
    return(table_up_ind_cpds)
  }else {
    return(table_up_od_cpds)
  }
  
}

get_table_gap_deg_cpds <- function(deg_inout, temp_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) {
  ## 2.gap
  gap_rows = temp_pct[which(temp_pct$ifgap == "gap"), ]
  gap_ind_cpds = c()
  gap_od_cpds = c()
  for (j in 1:length(gap_rows)) {
    tmp_cid = gap_rows[j, "cycid"]
    tmp_cin = gap_rows[j, "c_in"]
    tmp_cout = gap_rows[j, "c_out"]
    #2.1 gap_ind_cpds
    tmp_cin_ind = cycle_edgesucs_in[which((cycle_edgesucs_in$cycid == tmp_cid) & (cycle_edgesucs_in$node == tmp_cin)), "ind"]
    tmp_cout_ind = cycle_edgesucs_in[which((cycle_edgesucs_in$cycid == tmp_cid) & (cycle_edgesucs_in$node == tmp_cout)), "ind"]
    gap_ind_cpds = append(gap_ind_cpds, tmp_cin_ind)
    gap_ind_cpds = append(gap_ind_cpds, tmp_cout_ind)
    #2.2 gap_od_cpds
    tmp_cin_od = cycle_edgesucs_out[which((cycle_edgesucs_out$cycid == tmp_cid) & (cycle_edgesucs_out$node == tmp_cin)), "od"]
    tmp_cout_od = cycle_edgesucs_out[which((cycle_edgesucs_out$cycid == tmp_cid) & (cycle_edgesucs_out$node == tmp_cout)), "od"]
    gap_od_cpds = append(gap_od_cpds, tmp_cin_od)
    gap_od_cpds = append(gap_od_cpds, tmp_cout_od)
  }
  
  table_gap_ind_cpds = as.data.frame(table(gap_ind_cpds))
  colnames(table_gap_ind_cpds) = c("compound", "Freq")
  
  table_gap_od_cpds = as.data.frame(table(gap_od_cpds))
  colnames(table_gap_od_cpds) = c("compound", "Freq")

  #与这些化合物总的出现次数对比
  table_gap_ind_cpds = compare_with_all_Freq(table_gap_ind_cpds)
  table_gap_od_cpds = compare_with_all_Freq(table_gap_od_cpds)
  
  #删去never_considered_comp 中的化合物
  table_gap_ind_cpds = table_gap_ind_cpds[-which(table_gap_ind_cpds$compound %in% never_considered_comp_names),]
  table_gap_od_cpds = table_gap_od_cpds[-which(table_gap_od_cpds$compound %in% never_considered_comp_names),]
  
  
  if(deg_inout == "ind") {
    return(table_gap_ind_cpds)
  }else {
    return(table_gap_od_cpds)
  }
  
}


#########################################################################
# 1.统计
# 全部环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# 全部环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#########################################################################
## test
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")
## 1.读入环上每个cpd的出度od 和 入度ind 信息
#res_cycle_edgesucs_in res_cycle_edgesucs_out
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_cycle_edgesucs_in.RData")
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_cycle_edgesucs_out.RData")
#never_considered_comp_names
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/never_considered_compounds.RData")

##
#tumor_name = "COAD"
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
table_up_ind_cpds = list()
table_up_od_cpds = list()
table_gap_ind_cpds = list()
table_gap_od_cpds = list()

for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumor_name)
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]
  
  ##tumor
  temp_table_up_ind_cpds = get_table_up_deg_cpds("ind", tmp_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) 
  temp_table_up_od_cpds = get_table_up_deg_cpds("od", tmp_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) 
  temp_table_gap_ind_cpds = get_table_gap_deg_cpds("ind", tmp_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names)
  temp_table_gap_od_cpds = get_table_gap_deg_cpds("od", tmp_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names)
  
  table_up_ind_cpds[[tumor_name]] = temp_table_up_ind_cpds
  table_up_od_cpds[[tumor_name]] = temp_table_up_od_cpds
  table_gap_ind_cpds[[tumor_name]] = temp_table_gap_ind_cpds
  table_gap_od_cpds[[tumor_name]] = temp_table_gap_od_cpds
}

## save
save(table_up_ind_cpds, table_up_od_cpds, table_gap_ind_cpds, table_gap_od_cpds, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/all_degnode_freq.RData")






#########################################################################
# 2.统计
# ug环中 up行的cin&cout化合物 的deg_in deg_out的化合物出现频率
# ug环中 gap行的cin&cout化合物 的deg_in deg_out的化合物出现频率
#########################################################################
#gapup_cycle_chain_list.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData")

##
#tumor_name = "COAD"
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
ug_table_up_ind_cpds = list()
ug_table_up_od_cpds = list()
ug_table_gap_ind_cpds = list()
ug_table_gap_od_cpds = list()

for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumor_name)
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]
  ug_metric_pct = tmp_metric_pct[which(tmp_metric_pct$cycid %in% ug_c),]
  
  ##tumor
  temp_ug_table_up_ind_cpds = get_table_up_deg_cpds("ind", ug_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) 
  temp_ug_table_up_od_cpds = get_table_up_deg_cpds("od", ug_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names) 
  temp_ug_table_gap_ind_cpds = get_table_gap_deg_cpds("ind", ug_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names)
  temp_ug_table_gap_od_cpds = get_table_gap_deg_cpds("od", ug_metric_pct, cycle_edgesucs_in, cycle_edgesucs_out, never_considered_comp_names)
  
  ug_table_up_ind_cpds[[tumor_name]] = temp_ug_table_up_ind_cpds
  ug_table_up_od_cpds[[tumor_name]] = temp_ug_table_up_od_cpds
  ug_table_gap_ind_cpds[[tumor_name]] = temp_ug_table_gap_ind_cpds
  ug_table_gap_od_cpds[[tumor_name]] = temp_ug_table_gap_od_cpds
}

## save
save(ug_table_up_ind_cpds, ug_table_up_od_cpds, ug_table_gap_ind_cpds, ug_table_gap_od_cpds, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/ug_degnode_freq.RData")


