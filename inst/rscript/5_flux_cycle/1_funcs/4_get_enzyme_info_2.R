# 先执行: 4_get_enzyme_info_1.R


# Task:tumor_cyc_enzyme_count_type
#在ec中  在aup adown中 增加上  多少gene up

##################################### part1 #####################################
#############################################################################################



#在ec中  在aup adown中 增加上  每个 enzyme 多少gene up 多少gene down
add_enzyme_count_gene_num <- function(tumor_cyc_enzyme_count_type, tumor_name) {
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



# test 
# tumor_cyc_enzyme_count_type
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_enzyme_count_type = tumor_cyc_enzyme_count_type

for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  tumor_enzyme_info = as.data.frame(add_enzyme_count_gene_num(tumor_cyc_enzyme_count_type, tumors_array[i]))
  tumor_cyc_enzyme_count_type[[tumors_array[i]]] = tumor_enzyme_info
}

# save
save(tumor_cyc_enzyme_count_type, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")





