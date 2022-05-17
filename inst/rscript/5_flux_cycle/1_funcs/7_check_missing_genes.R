# disuse
#
# genecards.org
# 
# 
# "BCMO1"  "NAPRT1"
# 别名:原名
# BCMO1:BCO1 
# NAPRT1:NAPRT 
# 
# gene_names = rownames(d.matrix)
# "BCO1" %in% gene_names
# "NAPRT" %in% gene_names


#gene_missing_list
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_missing_list.RData")
# 
# library(readr)
# to_find_correct_symbol_v1_1<-read_csv("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/to_find_correct_symbol_v1.1.csv", col_names = FALSE)
# 
# #BCMO1
# to_find_correct_symbol_v1_1 = rbind(to_find_correct_symbol_v1_1,c("BCMO1","BCO1"))
# #NAPRT1
# to_find_correct_symbol_v1_1 = rbind(to_find_correct_symbol_v1_1,c("NAPRT1","NAPRT"))
# 
# to_find_correct_symbol = to_find_correct_symbol_v1_1
# 

#直接在to_find_correct_symbol_v1.1.csv文件中修改并添加2行即可



#########################################################################
#########################################################################

#tumor_cyc_enzyme_genecount_updown
# load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_genecount_and_pure.RData")
# 
# tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
# missing_genes_list = c()
# for (i in 1:length(tumors_array)) {
#   tumor_name = tumors_array[i]
#   temp_tumor_cyc_enzyme_genecount_updown = tumor_cyc_enzyme_genecount_updown[[tumor_name]]
#   #保存到tumor字典
#   for (j in 1:length(temp_tumor_cyc_enzyme_genecount_updown[,1])) {
#     if(is.na(temp_tumor_cyc_enzyme_genecount_updown[j, "gene_count"])) {
#       print(j)
#       temp_missing_gene = temp_tumor_cyc_enzyme_genecount_updown[j, "gene_name"]
#       missing_genes_list = append(missing_genes_list, temp_missing_gene)
#     }
#   }
# }
# 
# missing_genes_list = unique(missing_genes_list)
# 











