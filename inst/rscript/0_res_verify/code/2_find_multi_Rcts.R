
##################################### part1 find_multi rcts ###########################################
#######################################################################################################

# 6.
# 把12.31 joint meeting中
# Reaction中 占多个反应物的情况找出来
# 
# 把Rct占多个反应物的情况找出来
# 如 第248个 C05382->R01641
# 把这种全部找出来
# (多个反应物情况)
# 
# 目的是：cycid找出来. 汇报时remove 掉这些即可
#

#####################################################################################

get_cyc_multi_rct<-function(cycle_directed, cyc_multi_rct){
  for (i in 1:length(cycle_directed[,1])) {
    print(i)
    temp_ord_cpds_str = cycle_directed[i,"ord_cpds_str"]
    temp_ord_cpds_str = substring(temp_ord_cpds_str, 2,nchar(temp_ord_cpds_str)-1)
    temp_ord_cpds_str = unlist(strsplit(temp_ord_cpds_str,split = ","))[1]
    temp_ord_cpds_str = substring(temp_ord_cpds_str, 2,nchar(temp_ord_cpds_str)-1)
    
    temp_ord_cpds_str = unlist(strsplit(temp_ord_cpds_str,split = "->"))
    
    temp_ord_cpds_str = temp_ord_cpds_str[seq(2,length(temp_ord_cpds_str)-1,2)]
    
    temp_ocpds_table_arr = c()
    if (length(temp_ord_cpds_str) != length(unique(temp_ord_cpds_str))) {
      cyc_multi_rct[i,"if_multi_rct"] = 1
      
      temp_ocpds_table = table(temp_ord_cpds_str)
      temp_ocpds_table_df = as.data.frame(temp_ocpds_table)
      for (j in 1:length(temp_ocpds_table_df[,1])) {
        if (temp_ocpds_table_df[j,"Freq"] > 1) {
          ocpds_str = paste0(temp_ocpds_table_df[j,"temp_ord_cpds_str"],":",as.character(temp_ocpds_table_df[j,"Freq"]))
          temp_ocpds_table_arr = append(temp_ocpds_table_arr, ocpds_str) 
        }
      }
      cyc_multi_rct[i,"multi_rcts"] = paste(temp_ocpds_table_arr, collapse  = ";")
    }
    
  }
  return(cyc_multi_rct)
}


###
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_cycle_union_directed.RData")

#['C00101->R00946->C00440->R07168->C00143->R04125->C00101', 'C00101->R04125->C00143->R07168->C00440->R00946->C00101']

cyc_multi_rct = as.data.frame(cycle_directed[,1])
cyc_multi_rct[,"if_multi_rct"] = 0
cyc_multi_rct[,"multi_rcts"] = NA
colnames(cyc_multi_rct) = c("cycid","if_multi_rct","multi_rcts")


##
#test
cyc_multi_rct = get_cyc_multi_rct(cycle_directed, cyc_multi_rct)

save(cyc_multi_rct, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/0_res_verify/cyc_multi_rct.RData")

table(cyc_multi_rct[,"multi_rcts"])
# 0   1 
# 131 117 




###learn table
# temp_ocpds_table = table(temp_ord_cpds_str)
# temp_ocpds_table_df = as.data.frame(temp_ocpds_table)
# 
# names(temp_ocpds_table)
# as.numeric(temp_ocpds_table)


