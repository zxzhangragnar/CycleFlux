#表达量：count


#load and assign to variable
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# 1.get Tumor表达量：count
get_tumor_gene_count<-function(tumor_name, input_tumor_data)
{
  data_tumor = loadRData(input_tumor_data)

  gene_count = as.data.frame(data_tumor[,1])
  gene_count[,1] = 0
  colnames(gene_count) = c("count")

  mean_t = as.vector(rowMeans(data_tumor))
  gene_count[,"count"] = log2(mean_t+1)

  print(tumor_name)
  return(gene_count)
}

# 2.get Normal表达量：count
get_normal_gene_count<-function(tumor_name, input_normal_data)
{
  data_normal = loadRData(input_normal_data)

  gene_count = as.data.frame(data_normal[,1])
  gene_count[,1] = 0
  colnames(gene_count) = c("count")

  mean_n = as.vector(rowMeans(data_normal))
  gene_count[,"count"] = log2(mean_n+1)

  print(tumor_name)
  return(gene_count)
}



get_enzyme_info_0_main <- function(res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  #test
  #tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
  tumors_array = c(input_tumor_name)
  tumor_gene_count = list()
  normal_gene_count = list()
  for (i in 1:length(tumors_array)) {
    #保存到tumor字典
    t_count = as.data.frame(get_tumor_gene_count(tumors_array[i], input_tumor_data))
    tumor_gene_count[[tumors_array[i]]] = t_count

    #保存到normal字典
    n_count = as.data.frame(get_normal_gene_count(tumors_array[i], input_normal_data))
    normal_gene_count[[tumors_array[i]]] = n_count
  }


  #save
  res_sub_path = "5_flux_cycle/result_enzyme"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "tumor_normal_gene_count.RData")

  #save(tumor_gene_count, normal_gene_count, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_enzyme/tumor_normal_gene_count.RData")
  save(tumor_gene_count, normal_gene_count, file=res_file_path)

}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
#
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# input_tumor_name = "COAD"
#
# get_enzyme_info_0_main(res_path, input_tumor_name, input_tumor_data, input_normal_data)


