
##在执行其它 metric 文件之前 生成这两个结果作为对比

##################################################################################################################################################
## 化合物在全部环中的频率
get_all_cyc_cpds <- function() {
  # cycle_directed
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_cycle_union_directed.RData")
  
  all_cyc_cpds = c()
  for (i in 1:length(cycle_directed[,1])) {
    cpds_str = cycle_directed[i, "cpds"]
    cpds_str = substring(cpds_str, 2,nchar(cpds_str)-1)
    cpds_arr <- unlist(strsplit(cpds_str,split = ", "))
    for (j in 1:length(cpds_arr)) {
      tmp_cpd = substring(cpds_arr[j], 2,nchar(cpds_arr[j])-1)
      all_cyc_cpds = append(all_cyc_cpds, tmp_cpd)
    }
  }
  
  table_all_cyc_cpds = as.data.frame(table(all_cyc_cpds))
  table_all_cyc_cpds = table_all_cyc_cpds[order(table_all_cyc_cpds$Freq,decreasing = TRUE),]
  colnames(table_all_cyc_cpds) = c("compound", "Freq")
  return(table_all_cyc_cpds)
}




##################################################################################################################################################
## 化合物在全部度中的频率
get_all_deg_cpds <- function() {
  #res_cycle_edgesucs_in res_cycle_edgesucs_out
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_cycle_edgesucs_in.RData")
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_cycle_edgesucs_out.RData")
  
  all_deg_cpds = c()
  all_deg_cpds = append(all_deg_cpds, cycle_edgesucs_in$ind)
  all_deg_cpds = append(all_deg_cpds, cycle_edgesucs_out$od)
  table_all_deg_cpds = as.data.frame(table(all_deg_cpds))
  table_all_deg_cpds = table_all_deg_cpds[order(table_all_deg_cpds$Freq, decreasing = TRUE),]
  
  colnames(table_all_deg_cpds) = c("compound", "all_Freq")
  table_all_deg_cpds$all_Freq = table_all_deg_cpds$all_Freq*2
  return(table_all_deg_cpds)  
}


## save
table_all_cyc_cpds = get_all_cyc_cpds()
table_all_deg_cpds = get_all_deg_cpds()

save(table_all_cyc_cpds, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/table_all_cyc_cpds.RData")
save(table_all_deg_cpds, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/table_all_deg_cpds.RData")
