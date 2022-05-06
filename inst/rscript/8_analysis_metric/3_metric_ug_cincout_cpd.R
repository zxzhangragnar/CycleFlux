##################################################################################
##################################################################################
#
# 此文件的结果 和 6_count_gap_compound.R 结果的区别
# (频率统计)
# 1.统计全部环中的全部gap(或up)的反应 的c_in(或c_out)化合物 出现的次数频率
# 从全部环中 筛选出 ug的环
# 2.统计ug环中的全部gap(或up)的反应 的c_in(或c_out)化合物 出现的次数频率
# 统计ug环中的全部化合物 出现的次数频率
# 并进行对比
#
##################################################################################




################################################################################
## 统计ug环中的全部gap(或up)的反应 的c_in(或c_out)化合物 出现的次数频率

get_c_up_metric_table <- function(in_or_out, never_considered_comp_arr, temp_pct) {
  # table_all_cyc_cpds
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/table_all_cyc_cpds.RData")
  
  c_in_up_arr = temp_pct[which(temp_pct$ifup == "up"), in_or_out]
  ug_c_in_up_table = as.data.frame(table(c_in_up_arr))
  ug_c_in_up_table = ug_c_in_up_table[order(ug_c_in_up_table$Freq,decreasing = TRUE),]
  
  colnames(ug_c_in_up_table) = c("compound", "Freq")
  ug_c_in_up_table = ug_c_in_up_table[-which(ug_c_in_up_table$compound %in% never_considered_comp_arr),]
  
  ug_c_in_up_table[,"all_Freq"] = table_all_cyc_cpds[ug_c_in_up_table$compound,"Freq"]
  ug_c_in_up_table[,"ratio_Freq"] = ug_c_in_up_table[,"Freq"]/ug_c_in_up_table[,"all_Freq"]
  ug_c_in_up_table = ug_c_in_up_table[order(ug_c_in_up_table$ratio_Freq,ug_c_in_up_table$Freq,-ug_c_in_up_table$all_Freq,decreasing = T),]
  
  return(ug_c_in_up_table)
}
get_c_gap_metric_table <- function(in_or_out, never_considered_comp_arr, temp_pct) {
  # table_all_cyc_cpds
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/table_all_cyc_cpds.RData")
  
  c_in_gap_arr = temp_pct[which(temp_pct$ifgap == "gap"), in_or_out]
  ug_c_in_gap_table = as.data.frame(table(c_in_gap_arr))
  ug_c_in_gap_table = ug_c_in_gap_table[order(ug_c_in_gap_table$Freq,decreasing = TRUE),]
  
  colnames(ug_c_in_gap_table) = c("compound", "Freq")
  ug_c_in_gap_table = ug_c_in_gap_table[-which(ug_c_in_gap_table$compound %in% never_considered_comp_arr),]
  
  ug_c_in_gap_table[,"all_Freq"] = table_all_cyc_cpds[ug_c_in_gap_table$compound,"Freq"]
  ug_c_in_gap_table[,"ratio_Freq"] = ug_c_in_gap_table[,"Freq"]/ug_c_in_gap_table[,"all_Freq"]
  ug_c_in_gap_table = ug_c_in_gap_table[order(ug_c_in_gap_table$ratio_Freq,ug_c_in_gap_table$Freq,-ug_c_in_gap_table$all_Freq,decreasing = T),]
  
  return(ug_c_in_gap_table)
}


#cycle_edge_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")

#delete never consider cpds
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/never_considered_compounds.RData")

package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
library(readr)
never_considered_comp = read_csv(file.path(package_path, "tool_data/never_considered_comp.csv"))
never_considered_comp_arr = never_considered_comp$x

################################################################################
## 全部环中的 gap(up) cin(cout) 化合物频率
#tumor_name = "COAD"

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
all_c_in_up_table = list()
all_c_out_up_table = list()
all_c_in_gap_table = list()
all_c_out_gap_table = list()

for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumor_name)
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]
  
  ##tumor
  temp_all_c_in_up_table = get_c_up_metric_table("c_in", never_considered_comp_arr, tmp_metric_pct)
  temp_all_c_out_up_table = get_c_up_metric_table("c_out", never_considered_comp_arr, tmp_metric_pct)
  temp_all_c_in_gap_table = get_c_gap_metric_table("c_in", never_considered_comp_arr, tmp_metric_pct)
  temp_all_c_out_gap_table = get_c_gap_metric_table("c_out", never_considered_comp_arr, tmp_metric_pct)
  
  all_c_in_up_table[[tumor_name]] = temp_all_c_in_up_table
  all_c_out_up_table[[tumor_name]] = temp_all_c_out_up_table
  all_c_in_gap_table[[tumor_name]] = temp_all_c_in_gap_table
  all_c_out_gap_table[[tumor_name]] = temp_all_c_out_gap_table
}

## save
save(all_c_in_up_table, all_c_out_up_table, all_c_in_gap_table, all_c_out_gap_table, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/all_cincout_freq.RData")


################################################################################
## ug环中的 gap(up) cin(cout) 化合物频率
###
#gapup_cycle_chain_list.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData")
#tumor_name = "COAD"

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
ug_c_in_up_table = list()
ug_c_out_up_table = list()
ug_c_in_gap_table = list()
ug_c_out_gap_table = list()

for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumor_name)
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]
  ug_metric_pct = tmp_metric_pct[which(tmp_metric_pct$cycid %in% ug_c),]
  
  ##tumor
  temp_ug_c_in_up_table = get_c_up_metric_table("c_in", never_considered_comp_arr, ug_metric_pct)
  temp_ug_c_out_up_table = get_c_up_metric_table("c_out", never_considered_comp_arr, ug_metric_pct)
  temp_ug_c_in_gap_table = get_c_gap_metric_table("c_in", never_considered_comp_arr, ug_metric_pct)
  temp_ug_c_out_gap_table = get_c_gap_metric_table("c_out", never_considered_comp_arr, ug_metric_pct)
  
  ug_c_in_up_table[[tumor_name]] = temp_ug_c_in_up_table
  ug_c_out_up_table[[tumor_name]] = temp_ug_c_out_up_table
  ug_c_in_gap_table[[tumor_name]] = temp_ug_c_in_gap_table
  ug_c_out_gap_table[[tumor_name]] = temp_ug_c_out_gap_table
}

## save
save(ug_c_in_up_table, ug_c_out_up_table, ug_c_in_gap_table, ug_c_out_gap_table, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/ug_cincout_freq.RData")




