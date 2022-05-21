#disuse
pvalue_wilcox_test_all<-function(data_t, data_n) {
  pvalue_statr<-matrix(1,nrow(data_t),1)
  for(i in 1:nrow(data_t))
  {
    a<-data_t[i,]
    b<-data_n[i,]
    pvalue_statr[i,1]<-wilcox.test(a,b)$p.value
    if(i%%5000==1)
    {
      print(i)
    }
  }
  colnames(pvalue_statr)<-c("p.value")
  rownames(pvalue_statr)<-rownames(data_t)
  return(pvalue_statr)
}


# 详见：TCGA_tutorial_a.R
my_FoldChange_func <- function(data_t, data_n) {
  stat_r = matrix(1,nrow(data_t),1)
  mean_n = as.vector(rowMeans(data_n))+0.001
  mean_t = as.vector(rowMeans(data_t))+0.001
  
  library(gtools)
  fc <- foldchange(mean_t,mean_n)
  log_Fold_change = foldchange2logratio(fc, base = 2)
  
  stat_r[,1]<-log_Fold_change
  colnames(stat_r)<-c("foldchange")
  rownames(stat_r)<-rownames(data_t)
  stat_r = as.data.frame(stat_r)
  stat_r[which(stat_r > 0), "up_down"] = "Up"
  stat_r[which(stat_r$foldchange < 0), "up_down"] = "Down"
  
  stat_r[which(stat_r$foldchange > 0), "sign_up_down"] = 1
  stat_r[which(stat_r$foldchange < 0), "sign_up_down"] = -1
  
  
  #### pvalue ####
  pvalue_res = pvalue_wilcox_test_all(data_t, data_n)
  pvalue_statr = pvalue_res[,"p.value"]
  #### pvalue ####
  
  #### gene_id
  stat_r = cbind(rownames(stat_r), stat_r, pvalue_statr)
  ####
  colnames(stat_r) = c("gene_id", "foldchange", "up_down", "sign_up_down", "pvalue")
  
  return(stat_r)
}

#load and assign to variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



#########################################################################
get_gene_foldchange_main <- function(output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  
  #tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
  tumors_array = c(input_tumor_name)
  
  gene_foldchange_list = list()
  for (i in 1:length(tumors_array)) {
    data_tumor = loadRData(input_tumor_data[i])
    data_normal = loadRData(input_normal_data[i])
    #保存到tumor字典
    deg_foldchange = as.data.frame(my_FoldChange_func(data_tumor, data_normal))
    gene_foldchange_list[[tumors_array[i]]] = deg_foldchange
  }
  
  
  #save
  res_sub_path = "3_flux_subnet/result_tool"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "gene_foldchange_list.RData")
  
  #save(gene_foldchange_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_foldchange_list.RData")
  save(gene_foldchange_list, file=res_file_path)
  
}
 




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# 
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# input_tumor_name = "COAD"
# 
# get_gene_foldchange_main(output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)



