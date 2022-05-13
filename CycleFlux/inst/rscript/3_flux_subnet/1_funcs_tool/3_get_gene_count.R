
#######################################################################################################
# 目标：为tumor_cyc_gap_foldchangsign.RData添加c_in和c_out两列

#######################################################################################################




##
get_pct_statr_info<-function(data_tumor, data_normal)
{
  ##
  data_t = data_tumor
  data_n = data_normal
  #3.t_gene_mean_val
  t_gene_mean_val = as.vector(rowMeans(data_t))
  #4.n_gene_mean_val
  n_gene_mean_val = as.vector(rowMeans(data_n))
  #5.t_gene_0_pct
  f_count_zero<-function(x) sum(x==0)
  t_gene_0_num = as.vector(apply(data_t,1,f_count_zero))
  t_gene_0_pct = t_gene_0_num/ncol(data_t)

  #6.n_gene_0_pct
  n_gene_0_num = as.vector(apply(data_n,1,f_count_zero))
  n_gene_0_pct = n_gene_0_num/ncol(data_n)

  #7.t_gene_log_lss1_pct
  f_count_log_lss1<-function(x) sum(log(x+1)<1)
  t_gene_log_lss1_num = as.vector(apply(data_t,1,f_count_log_lss1))
  t_gene_log_lss1_pct = t_gene_log_lss1_num/ncol(data_t)

  #8.n_gene_log_lss1_pct
  n_gene_log_lss1_num = as.vector(apply(data_n,1,f_count_log_lss1))
  n_gene_log_lss1_pct = n_gene_log_lss1_num/ncol(data_n)

  #put
  pct_statr_info = as.data.frame(matrix(1,nrow(data_t),1))
  rownames(pct_statr_info)<-rownames(data_t)
  pct_statr_info[,1]<-t_gene_mean_val
  colnames(pct_statr_info)<-c("t_gene_mean_val")
  pct_statr_info[, "n_gene_mean_val"] = n_gene_mean_val
  pct_statr_info[, "t_gene_0_pct"] = t_gene_0_pct
  pct_statr_info[, "n_gene_0_pct"] = n_gene_0_pct
  pct_statr_info[, "t_gene_log_lss1_pct"] = t_gene_log_lss1_pct
  pct_statr_info[, "n_gene_log_lss1_pct"] = n_gene_log_lss1_pct

  return(pct_statr_info)
}


#load and assign to variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}




#########################################################################
get_gene_count_main <- function(output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {

  data_tumor = loadRData(input_tumor_data)
  data_normal = loadRData(input_normal_data)

  tumors_array = c(input_tumor_name)

  gene_meanval_list = list()
  # old statr
  for (i in 1:length(tumors_array)) {
    #保存到tumor字典
    statr = as.data.frame(get_pct_statr_info(data_tumor, data_normal))
    gene_meanval_list[[tumors_array[i]]] = statr
  }

  #save
  res_sub_path = "3_flux_subnet/result_tool"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "gene_meanval_list.RData")

  #save(gene_meanval_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/3_flux_subnet/result_tool/gene_meanval_list.RData")
  save(gene_meanval_list, file=res_file_path)

}



# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
#
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# input_tumor_name = "COAD"
#
# get_gene_count_main(output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)







