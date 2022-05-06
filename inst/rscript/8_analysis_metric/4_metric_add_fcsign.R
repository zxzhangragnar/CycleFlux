
##################################### part1 flux gap foldchangesign ###########################################
#######################################################################################################

# （数量统计）
# 这些涉及到gap(gap左右端点)的cin和cout的个数
# 及all和gap的table出现数量统计的对比
# 增加如下3项:
# fcsign_sum (某个化合物 foldchangesign 的总和)
# up_count   (某个化合物 foldchangesign>0 的反应的个数)
# down_count (某个化合物 foldchangesign<0 的反应的个数)


## 先运行 3_metric_ug_cincout_cpd.R

##############################################################################################
add_table_freq_fcsign_sum_cin <- function(table_temp, temp_pct) {
  #metric
  table_temp[,"fcsign_sum"] = 0
  table_temp[,"up_count"] = 0
  table_temp[,"down_count"] = 0
  # someC_1->Cin->rid->Cout->someC_2 
  # 原来做Cin的，现在做Cout，为了找someC_1
  for (i in 1:length(table_temp[,1])) {
    tmp_cin = table_temp[i,"compound"]
    tmp_df = temp_pct[which(temp_pct$c_out == tmp_cin),c("rid","foldchangesign", "ifgap", "ifup")]
    tmp_df = unique(tmp_df)
    fcsign_sum = 0
    up_count = 0
    down_count = 0
    
    ##这个化合物 出去的反应 为gap(up)的有多少
    cpd_ifgap_count = 0
    cpd_ifup_count = 0
    
    for (j in 1:length(tmp_df[,1])) {
      fcsign_sum = fcsign_sum + tmp_df[j,"foldchangesign"]
      if(tmp_df[j,"foldchangesign"] == 0) {
        print("NA updown")
      }else if(tmp_df[j,"foldchangesign"] > 0) {
        up_count = up_count + 1
      }else if(tmp_df[j,"foldchangesign"] < 0) {
        down_count = down_count + 1
      }
      
      if(tmp_df[j,"ifgap"] == "gap") {
        cpd_ifgap_count = cpd_ifgap_count + 1
      }
      
      if(tmp_df[j,"ifup"] == "up") {
        cpd_ifup_count = cpd_ifup_count + 1
      }
    }
    #put into df
    table_temp[i,"fcsign_sum"] = fcsign_sum/length(tmp_df[,1])
    table_temp[i,"up_count"] = up_count/length(tmp_df[,1])
    table_temp[i,"down_count"] = down_count/length(tmp_df[,1])
    table_temp[i,"cpd_ifgap_count"] = cpd_ifgap_count/length(tmp_df[,1])
    table_temp[i,"cpd_ifup_count"] = cpd_ifup_count/length(tmp_df[,1])
    
  }
  #put in
  return(table_temp)
}

add_table_freq_fcsign_sum_cout <- function(table_temp, temp_pct) {
  #metric
  table_temp[,"fcsign_sum"] = 0
  table_temp[,"up_count"] = 0
  table_temp[,"down_count"] = 0
  # someC_1->Cin->rid->Cout->someC_2 
  # 原来做Cout的，现在做Cin，为了找someC_2
  for (i in 1:length(table_temp[,1])) {
    tmp_cout = table_temp[i,"compound"]
    tmp_df = temp_pct[which(temp_pct$c_in == tmp_cout),c("rid", "foldchangesign", "ifgap", "ifup")]
    tmp_df = unique(tmp_df)
    fcsign_sum = 0
    up_count = 0
    down_count = 0
    
    ##这个化合物 出去的反应 为gap(up)的有多少
    cpd_ifgap_count = 0
    cpd_ifup_count = 0
    
    for (j in 1:length(tmp_df[,1])) {
      fcsign_sum = fcsign_sum + tmp_df[j,"foldchangesign"]
      if(tmp_df[j,"foldchangesign"] == 0) {
        print("NA updown")
      }else if(tmp_df[j,"foldchangesign"] > 0) {
        up_count = up_count + 1
      }else if(tmp_df[j,"foldchangesign"] < 0) {
        down_count = down_count + 1
      }
      
      if(tmp_df[j,"ifgap"] == "gap") {
        cpd_ifgap_count = cpd_ifgap_count + 1
      }
      
      if(tmp_df[j,"ifup"] == "up") {
        cpd_ifup_count = cpd_ifup_count + 1
      }      
    }
    
    #put into df
    table_temp[i,"fcsign_sum"] = fcsign_sum/length(tmp_df[,1])
    table_temp[i,"up_count"] = up_count/length(tmp_df[,1])
    table_temp[i,"down_count"] = down_count/length(tmp_df[,1])
    table_temp[i,"cpd_ifgap_count"] = cpd_ifgap_count/length(tmp_df[,1])
    table_temp[i,"cpd_ifup_count"] = cpd_ifup_count/length(tmp_df[,1])
  }
  #put in
  return(table_temp)
}

#test
#cycle_edge_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")


####################################################################################
#all_cincout_freq
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/all_cincout_freq.RData")

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumors_array[i])
  ##tumor
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]
  all_c_in_up_table[[tumor_name]] = add_table_freq_fcsign_sum_cin(all_c_in_up_table[[tumor_name]], tmp_metric_pct)
  all_c_out_up_table[[tumor_name]] = add_table_freq_fcsign_sum_cout(all_c_out_up_table[[tumor_name]], tmp_metric_pct)
  all_c_in_gap_table[[tumor_name]] = add_table_freq_fcsign_sum_cin(all_c_in_gap_table[[tumor_name]], tmp_metric_pct)
  all_c_out_gap_table[[tumor_name]] = add_table_freq_fcsign_sum_cout(all_c_out_gap_table[[tumor_name]], tmp_metric_pct)
}

## save
save(all_c_in_up_table, all_c_out_up_table, all_c_in_gap_table, all_c_out_gap_table, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/all_cincout_freq.RData")


####################################################################################

#ug_cincout_freq
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/ug_cincout_freq.RData")

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  print(tumors_array[i])
  ##tumor
  tmp_metric_pct = cycle_edge_flux_list[[tumor_name]]

  ug_c_in_up_table[[tumor_name]] = add_table_freq_fcsign_sum_cin(ug_c_in_up_table[[tumor_name]], tmp_metric_pct)
  ug_c_out_up_table[[tumor_name]] = add_table_freq_fcsign_sum_cout(ug_c_out_up_table[[tumor_name]], tmp_metric_pct)
  ug_c_in_gap_table[[tumor_name]] = add_table_freq_fcsign_sum_cin(ug_c_in_gap_table[[tumor_name]], tmp_metric_pct)
  ug_c_out_gap_table[[tumor_name]] = add_table_freq_fcsign_sum_cout(ug_c_out_gap_table[[tumor_name]], tmp_metric_pct)
}

## save
save(ug_c_in_up_table, ug_c_out_up_table, ug_c_in_gap_table, ug_c_out_gap_table, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/8_analysis_metric/metric_result/ug_cincout_freq.RData")
