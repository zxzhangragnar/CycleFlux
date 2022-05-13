
#######################################################################################################
# 目标：为cycle_flux_list
# 更新原来的ifgap项的内容
# 并增加一列ifup

#######################################################################################################




#####################################################
get_edge_max_gene_foldchange<-function(tumor_name, cycle_edge_flux_list) {

  gene_fc_arr<-cycle_edge_flux_list[[tumor_name]][,"gene_fold_change"] # gene_fold_change
  max_gene_fc_arr<-c()
  for(i in 1:length(gene_fc_arr))
  {
    temp_gene_fcs<-unlist(strsplit(gene_fc_arr[i],";"))
    temp_val<-c()
    for(j in 1:length(temp_gene_fcs))
    {
      temp_val<-c(temp_val,unlist(strsplit(temp_gene_fcs[j],":"))[2])
    }
    temp_val<-as.numeric(temp_val)
    max_gene_fc_arr<-c(max_gene_fc_arr,max(temp_val,na.rm=T)) #添加temp_val中的最大者
  }

  return(max_gene_fc_arr)
}

get_edge_max_gene_meanval<-function(tumor_name, cycle_edge_flux_list) {

  gene_meanval_arr<-cycle_edge_flux_list[[tumor_name]][,"t_gene_mean_val"] # t_gene_mean_val
  max_gene_meanval_arr<-c()
  for(i in 1:length(gene_meanval_arr))
  {
    temp_gene_meanvals<-unlist(strsplit(gene_meanval_arr[i],";"))
    temp_val<-c()
    for(j in 1:length(temp_gene_meanvals))
    {
      temp_val<-c(temp_val,unlist(strsplit(temp_gene_meanvals[j],":"))[2])
    }
    temp_val<-as.numeric(temp_val)
    max_gene_meanval_arr<-c(max_gene_meanval_arr,max(temp_val,na.rm=T)) #添加temp_val中的最大者
  }

  return(max_gene_meanval_arr)
}


get_gap_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  edge_max_gene_foldchange = get_edge_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  edge_max_gene_meanval = get_edge_max_gene_meanval(tumor_name, cycle_edge_flux_list)

  fc_val = edge_max_gene_foldchange
  gn_val = edge_max_gene_meanval

  gap_c<-unique(cycle_edge_flux_list[[tumor_name]][which((fc_val<(-2))|(gn_val<1.71)),1])
  return(gap_c)
}



get_up_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  edge_max_gene_foldchange = get_edge_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  edge_max_gene_meanval = get_edge_max_gene_meanval(tumor_name, cycle_edge_flux_list)

  fc_val = edge_max_gene_foldchange
  gn_val = edge_max_gene_meanval

  up_c<-unique(cycle_edge_flux_list[[tumor_name]][which((fc_val>(1))&(gn_val>10)),1])#p.value of diff<0.001, log(fc)>0.5, and the gene show significant up regulation are larger than 10
  return(up_c)
}



##################################################################
refesh_cycle_flux_list<-function(tumor_name, cycle_flux_list, cycle_edge_flux_list){

  gap_c = get_gap_c_cycid(tumor_name, cycle_edge_flux_list)
  up_c = get_up_c_cycid(tumor_name, cycle_edge_flux_list)

  temp_tumor_df = cycle_edge_flux_list[[tumor_name]]
  temp_tumor_updownsign_df = temp_tumor_df[,c("cycid","foldchangesign")]
  #updownsign_df_sum = aggregate(temp_tumor_updownsign_df$foldchangesign, list(temp_tumor_updownsign_df$cycid), sum)
  updownsign_df_mean = aggregate(temp_tumor_updownsign_df$foldchangesign, list(temp_tumor_updownsign_df$cycid), mean)
  colnames(updownsign_df_mean) = c("cycid", "updownsign")
  cycle_flux_list[[tumor_name]][, "updownsign"] = updownsign_df_mean$updownsign

  cycle_flux_list[[tumor_name]][, "ifup"] = "normal"
  cycle_flux_list[[tumor_name]][, "ifgap"] = "normal"
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% up_c), "ifup"] = "up"
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% gap_c), "ifgap"] = "gap"

  return(cycle_flux_list)
}


cyc_info_gap_main <- function(res_path, input_tumor_name) {

  load(file.path(res_path, "/5_flux_cycle/result_final/cycle_flux_list.RData"))
  load(file.path(res_path, "/4_flux_edge/result_final/cycle_edge_flux_list.RData"))

  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_flux_list = refesh_cycle_flux_list(tumor_name, cycle_flux_list, cycle_edge_flux_list)
  }

  #save
  res_sub_path = "5_flux_cycle/result_final"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_flux_list.RData")
  save(cycle_flux_list, file=res_file_path)

}





# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
#
# cyc_info_gap_main(res_path, input_tumor_name)
#






