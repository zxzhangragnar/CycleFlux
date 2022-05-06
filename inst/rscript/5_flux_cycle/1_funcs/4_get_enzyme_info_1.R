##############################################################
# 目的:
# tumor_normal_gene_count.RData各基因的表达量
##############################################################



##################################### part1 考虑count #####################################
#############################################################################################

#########################################                 
#find the enzyme's gene's express count 
#表达量：count

# 1.get Tumor表达量：count
get_tumor_gene_count<-function(tumor_name)
{
  #tumor_name = "COAD"
  tumor_filename = paste0("TCGA-",tumor_name,".RData") #TCGA-COAD.RData
  load_str = paste0("E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/", tumor_filename)
  load(load_str)
  tcga_tname = d.matrix
  
  gene_count = as.data.frame(tcga_tname[,1])
  gene_count[,1] = 0
  colnames(gene_count) = c("count")
  
  mean_t = as.vector(rowMeans(tcga_tname))
  gene_count[,"count"] = log2(mean_t+1)
    
  print(tumor_name)
  return(gene_count)
}

# 2.get Normal表达量：count
get_normal_gene_count<-function(tumor_name)
{
  #init
  #tumor_name = "COAD"
  tumor_N_filename = paste0("TCGA-",tumor_name,"_N.RData") #TCGA-COAD_N.RData
  load_N_str = paste0("E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/", tumor_N_filename)
  load(load_N_str)
  tcga_tname_N = d.matrix
  
  gene_count = as.data.frame(tcga_tname_N[,1])
  gene_count[,1] = 0
  colnames(gene_count) = c("count")
  
  mean_n = as.vector(rowMeans(tcga_tname_N))
  gene_count[,"count"] = log2(mean_n+1)
  
  print(tumor_name)
  return(gene_count)
}

#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_gene_count = list()
normal_gene_count = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  t_count = as.data.frame(get_tumor_gene_count(tumors_array[i]))
  tumor_gene_count[[tumors_array[i]]] = t_count
  
  #保存到normal字典
  n_count = as.data.frame(get_normal_gene_count(tumors_array[i]))
  normal_gene_count[[tumors_array[i]]] = n_count
}

save(tumor_gene_count, normal_gene_count, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_normal_gene_count.RData")








##################################################################################
##################################################################################


#########################################                 
# enzy 规则 造规则
# enzy找：1.有多少酶中全部基因 up down
# 2.有多少酶表达量最高的gene up down

# 给cycle_edge_expression增加一条 gene_count
# 将tumor_count\normal_count\tumor_FC添加到 tumor_cyc_enzyme_genecount_updown 中
get_cyc_enzyme_genecount<-function(tumor_name, cyc_gene_not_in_tgca_but_in_tofind) {
  #init
  #tumor_gene_count normal_gene_count
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_normal_gene_count.RData")
  #dp_part_net
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData")
  #cycle_edge_expression
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/cycle_edge_expression.RData")
  #gene_foldchange_list
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_foldchange_list.RData")
  
  tumor_gene_count = tumor_gene_count
  normal_gene_count = normal_gene_count
  dp_part_net = dp_part_net
  cycle_edge_expression = cycle_edge_expression
  tumor_gene_FC = gene_foldchange_list[[tumor_name]]
  rownames_tumor_gene_FC = rownames(tumor_gene_FC)
  
  # gene_count_df: tumor_count
  gene_count_df = tumor_gene_count[[tumor_name]]
  # normal_gene_df:normal_count  
  normal_gene_df = normal_gene_count[[tumor_name]]
  # gene_FC_df: tumor_gene_FC
  gene_FC_df = as.data.frame(tumor_gene_FC[,"foldchange"])
  colnames(gene_FC_df) = c("foldchange")
  rownames(gene_FC_df) = rownames_tumor_gene_FC
  
  cyc_enzyme_genecount = data.frame()
  for (i in 1:length(cycle_edge_expression[,1])) {
    #这个环的这个反应(相同cin,cout)所对应的所有酶
    temp_cin = cycle_edge_expression[i,"c_in"]
    temp_cout = cycle_edge_expression[i,"c_out"]
    
    temp_enzyme_array = cycle_edge_expression[i,"enzyme"]
    temp_enzyme_array <- unlist(strsplit(temp_enzyme_array,split = ";"))
    
    for (j in 1:length(temp_enzyme_array)) {
      temp_enzyme = temp_enzyme_array[j]
      #这个酶所对应的所有基因
      temp_gene_array = as.character(unlist(dp_part_net[which((dp_part_net$EC == temp_enzyme) & ((dp_part_net$C_in == temp_cin)|(dp_part_net$C_in == temp_cout)) & ((dp_part_net$C_out == temp_cin)|(dp_part_net$C_out == temp_cout))), "Gene_symbol"]))

      temp_gene_array = unique(temp_gene_array)
      temp_gene_array <- unlist(strsplit(temp_gene_array,split = ";"))
      for (k in 1:length(temp_gene_array)) {
        temp_gene = temp_gene_array[k]
        temp_gene_old = temp_gene
        if (temp_gene %in% cyc_gene_not_in_tgca_but_in_tofind) {
          temp_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_gene),2][[1]]
        }
        temp_row = data.frame()
        temp_row[1,"cycid"] = cycle_edge_expression[i,1]
        temp_row[1,"rid"] = cycle_edge_expression[i,2]
        temp_row[1,"enzyme"] = temp_enzyme
        temp_row[1,"gene_name"] = temp_gene_old
        #gene_count(tumor_count)
        temp_row[1,"gene_count"] = gene_count_df[c(temp_gene),]
        #normal_count
        temp_row[1,"normal_gene_count"] = normal_gene_df[c(temp_gene),]
        #FoldChange
        temp_row[1,"FoldChange"] = gene_FC_df[c(temp_gene),]
        
        cyc_enzyme_genecount = rbind(cyc_enzyme_genecount, temp_row)  
      }
      
    }
  }
  print(tumor_name)
  return(cyc_enzyme_genecount)
}

# test
library(readr)
to_find_correct_symbol_v1_1 <- read_csv("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/to_find_correct_symbol_v1.1.csv", col_names = FALSE)

load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_missing_list.RData")
cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_enzyme_genecount = list()
tumor_cyc_enzyme_genecount_pure = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  enzyme_genecount = as.data.frame(get_cyc_enzyme_genecount(tumors_array[i], cyc_gene_not_in_tgca_but_in_tofind))
  enzyme_genecount_pure = unique(enzyme_genecount[,c("enzyme","gene_name","gene_count","normal_gene_count","FoldChange")])
  
  tumor_cyc_enzyme_genecount[[tumors_array[i]]] = enzyme_genecount
  tumor_cyc_enzyme_genecount_pure[[tumors_array[i]]] = enzyme_genecount_pure
}





#############################################################                 
# 计算环中某个酶ec的某个基因gene的sign 上调下调值
get_tumor_gap_enzyme_genecount_updowntcga<-function(tumor_name, datatype, cyc_gene_not_in_tgca_but_in_tofind)
{
  load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_foldchange_list.RData")
  gene_foldchange_list = gene_foldchange_list
  
  #tumor_name = "COAD"
  print(tumor_name)
  stat_r = as.data.frame(gene_foldchange_list[tumor_name])
  names(stat_r) = c("gene_id","log2FoldChange","up_down","sign_up_down","pvalue")
  et = stat_r
  etSig = et
  
  
  #based on tumor_cyc_enzyme_genecount
  if (datatype == "pure") {
    cyc_enzyme_genecount_df = tumor_cyc_enzyme_genecount_pure[[tumor_name]]
  }else {
    cyc_enzyme_genecount_df = tumor_cyc_enzyme_genecount[[tumor_name]]
  }
  cyc_up_down = cyc_enzyme_genecount_df
  cyc_up_down[,"sign"] = 0
  for (i in 1:length(cyc_up_down[,1])) {
    # 此处加上判断，否则因找不到基因而报错  Error in x[[jj]][iseq] <- vjj : 更换参数长度为零 
    temp_str_gene = cyc_up_down[i, "gene_name"] 
    
    #此处应该为 cyc_gene_not_in_tgca_but_in_tofind 
    #否则会报错 Error in x[[jj]][iseq] <- vjj : 更换参数长度为零
    if (temp_str_gene %in% cyc_gene_not_in_tgca_but_in_tofind) {
      temp_str_gene = to_find_correct_symbol_v1_1[which(to_find_correct_symbol_v1_1[,1] == temp_str_gene),2][[1]]
    }
    genesign = as.double(etSig[temp_str_gene,"sign_up_down"])
    #cyc_up_down[i,"sign"] = cyc_up_down[i,"sign"] + genesign
    cyc_up_down[i,"sign"] = genesign
    
  } 
  cyc_up_down[which(cyc_up_down$sign > 0),"updown"] = "up"
  cyc_up_down[which(cyc_up_down$sign < 0),"updown"] = "down"
  # write gap
  #cyc_up_down[which(abs(cyc_up_down$sign) >= 3),"ifgap"] = "gap"
  
  return(cyc_up_down)
}

#########################################
# test
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_missing_list.RData")
cyc_gene_not_in_tgca_but_in_tofind = gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]]

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_enzyme_genecount_updown = list()
tumor_cyc_enzyme_genecount_updown_pure = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  ec_count_updown = get_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i], "normal", cyc_gene_not_in_tgca_but_in_tofind)
  tumor_cyc_enzyme_genecount_updown[[tumors_array[i]]] = ec_count_updown
  
  ec_count_updown_pure = get_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i], "pure", cyc_gene_not_in_tgca_but_in_tofind)
  tumor_cyc_enzyme_genecount_updown_pure[[tumors_array[i]]] = ec_count_updown_pure
}

#save
save(tumor_cyc_enzyme_genecount, tumor_cyc_enzyme_genecount_pure, tumor_cyc_enzyme_genecount_updown, tumor_cyc_enzyme_genecount_updown_pure, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_genecount_and_pure.RData")




#############################################################                 
# 将每个酶的各个gene按count大小排序
order_tumor_gap_enzyme_genecount_updowntcga<-function(tumor_name)
{
  #new data frame
  count_pure_temp_ordered = data.frame()
  
  count_pure_temp = tumor_cyc_enzyme_genecount_updown_pure[[tumor_name]]
  
  #enzyme_table = as.character(as.data.frame(table(count_pure_temp$enzyme))[,1])
  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i] 
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
    temp_enzyme_df = temp_enzyme_df[order(temp_enzyme_df$gene_count, decreasing = T),]
    #add to new data frame
    count_pure_temp_ordered = rbind(count_pure_temp_ordered, temp_enzyme_df)  
  }
  
  return(count_pure_temp_ordered)
}

#test
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_genecount_and_pure.RData")
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_count_pure_temp_ordered = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  count_pure_temp_ordered = order_tumor_gap_enzyme_genecount_updowntcga(tumors_array[i])
  tumor_cyc_count_pure_temp_ordered[[tumors_array[i]]] = count_pure_temp_ordered
}

#do not need to be saved





#############################################################                 
# 分成几类type:
# 1.酶中所有基因gene全上调，全下调
# 2.酶中表达量最高的基因gene上调，下调

#get type function
get_tumor_gap_enzyme_genecount_updowntcga_type<-function(tumor_name) {
  
  ## highstcount_gene
  highstcount_gene_up = data.frame()
  highstcount_gene_down = data.frame()
  count_pure_temp = tumor_cyc_count_pure_temp_ordered[[tumor_name]]
  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i] 
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
    if (!is.na(temp_enzyme_df[1,"updown"])) { #BCMO1 NA
      if (temp_enzyme_df[1,"updown"] == "up") { # row 1 is highstcount
        #add to new data frame
        highstcount_gene_up = rbind(highstcount_gene_up, temp_enzyme_df)      
      }else if(temp_enzyme_df[1,"updown"] == "down") {
        #add to new data frame
        highstcount_gene_down = rbind(highstcount_gene_down, temp_enzyme_df)      
      }    
    }
  }
  
  ## all_gene
  all_gene_up = data.frame()
  all_gene_down = data.frame()
  count_pure_temp = tumor_cyc_count_pure_temp_ordered[[tumor_name]]
  enzyme_table = unique(count_pure_temp$enzyme)
  for (i in 1:length(enzyme_table)) {
    temp_enzyme = enzyme_table[i] 
    temp_enzyme_df = count_pure_temp[which(count_pure_temp$enzyme == temp_enzyme),]
    is_all_up = TRUE
    is_all_down = TRUE
    for (j in 1:length(temp_enzyme_df)) {
      if (!is.na(temp_enzyme_df[j,"updown"])) { #BCMO1 NA
        if (temp_enzyme_df[j,"updown"] != "up") { # row 1 is highstcount
          is_all_up = FALSE
        }else if(temp_enzyme_df[j,"updown"] != "down") {
          is_all_down = FALSE
        }    
      }    
    }
    if (is_all_up) {
      all_gene_up = rbind(all_gene_up, temp_enzyme_df)
    }
    if (is_all_down) {
      all_gene_down = rbind(all_gene_down, temp_enzyme_df)
    }
  }
  
  type_result = list()
  type_result[["highstcount_gene_up"]] = unique(highstcount_gene_up$enzyme)
  type_result[["highstcount_gene_down"]] = unique(highstcount_gene_down$enzyme)
  type_result[["all_gene_up"]] = unique(all_gene_up$enzyme)
  type_result[["all_gene_down"]] = unique(all_gene_down$enzyme)
  
  return(type_result)
  
}



#test
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_count_pure_type = list()# 二级数组list
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  type_result = get_tumor_gap_enzyme_genecount_updowntcga_type(tumors_array[i])
  tumor_cyc_count_pure_type[[tumors_array[i]]] = type_result #二维数组
}


#########################################   
# save part2
save(tumor_cyc_count_pure_type, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_count_pure_type.RData")




################################### merge part2 ############################################                 
###############################################################################
#结合通过count计算的4种情况分类
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_count_pure_type.RData")

#tumor_cyc_enzyme_genecount_updown
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_genecount_and_pure.RData")

# type_result[["highstcount_gene_up"]] = unique(highstcount_gene_up$enzyme)
# type_result[["highstcount_gene_down"]] = unique(highstcount_gene_down$enzyme)
# type_result[["all_gene_up"]] = unique(all_gene_up$enzyme)
# type_result[["all_gene_down"]] = unique(all_gene_down$enzyme)



merge_enzyme_count_ptype<-function(tumor_name) {
  temp_ec_count = tumor_cyc_enzyme_genecount_updown[[tumor_name]]
  temp_ec_ptype = tumor_cyc_count_pure_type[[tumor_name]]
  
  temp_ec_count[,"ec_situation"] = 0
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["highstcount_gene_up"]]),"ec_situation"] = "hup"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["highstcount_gene_down"]]),"ec_situation"] = "hdown"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["all_gene_up"]]),"ec_situation"] = "aup"
  temp_ec_count[which(temp_ec_count$enzyme %in% temp_ec_ptype[["all_gene_down"]]),"ec_situation"] = "adown"
  
  return(temp_ec_count)
}


#merge
tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
tumor_cyc_enzyme_count_type = list()
for (i in 1:length(tumors_array)) {
  #保存到tumor字典
  temp_ec_count = as.data.frame(merge_enzyme_count_ptype(tumors_array[i]))
  tumor_cyc_enzyme_count_type[[tumors_array[i]]] = temp_ec_count
}



# save merge part2
save(tumor_cyc_enzyme_count_type, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_cyc_enzyme_count_type.RData")

















# debug
# d.matrix[c("GSTT2B"),]
# GSTT2B
# 
# # PPAP2A not in cyc_gene_expressinfo and not in to_find
# PPAP2A
# d.matrix[c("PPAP2A"),]
# "PPAP2A" %in% cyc_gene_expressinfo
# 
# test_vector = c(1,2,3,4,5)




