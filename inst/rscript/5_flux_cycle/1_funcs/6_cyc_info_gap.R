
#######################################################################################################
# 目标：为cycle_flux_list
# 更新原来的ifgap项的内容
# 并增加一列ifup

#######################################################################################################




####################
# functions

get_gap_c_max_gene_foldchange<-function(tumor_name, cycle_edge_flux_list) {
  
  aaa<-cycle_edge_flux_list[[tumor_name]][,"gene_fold_change"]
  ddd<-c()
  for(i in 1:length(aaa))
  {
    bbb<-unlist(strsplit(aaa[i],";"))
    ccc<-c()
    ccc2<-c()
    for(j in 1:length(bbb))
    {
      ccc<-c(ccc,unlist(strsplit(bbb[j],":"))[1])
      ccc2<-c(ccc2,unlist(strsplit(bbb[j],":"))[2])
    }
    ccc2<-as.numeric(ccc2)
    names(ccc2)<-ccc
    ddd<-c(ddd,max(ccc2,na.rm=T)) #添加ccc2中的最大者
  }
  
  return(ddd)
}

get_gap_c_max_gene_meanval<-function(tumor_name, cycle_edge_flux_list) {
  
  ddd2<-c()
  aaa<-cycle_edge_flux_list[[tumor_name]][,"t_gene_mean_val"]
  for(i in 1:length(aaa))
  {
    bbb<-unlist(strsplit(aaa[i],";"))
    ccc<-c()
    ccc2<-c()
    for(j in 1:length(bbb))
    {
      ccc<-c(ccc,unlist(strsplit(bbb[j],":"))[1])
      ccc2<-c(ccc2,unlist(strsplit(bbb[j],":"))[2])
    }
    ccc2<-as.numeric(ccc2)
    names(ccc2)<-ccc
    ddd2<-c(ddd2,max(ccc2,na.rm=T))
  }
  
  return(ddd2)
}


get_gap_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  c_max_gene_foldchange = get_gap_c_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  c_max_gene_meanval = get_gap_c_max_gene_meanval(tumor_name, cycle_edge_flux_list)
  
  ddd = c_max_gene_foldchange
  ddd2 = c_max_gene_meanval
  
  gap_c<-unique(cycle_edge_flux_list[[tumor_name]][which((ddd<(-2))|(ddd2<1.71)),1])
  return(gap_c)
}



get_up_c_cycid<-function(tumor_name, cycle_edge_flux_list) {
  c_max_gene_foldchange = get_gap_c_max_gene_foldchange(tumor_name, cycle_edge_flux_list)
  c_max_gene_meanval = get_gap_c_max_gene_meanval(tumor_name, cycle_edge_flux_list)
  
  ddd = c_max_gene_foldchange
  ddd2 = c_max_gene_meanval
  
  up_c<-unique(cycle_edge_flux_list[[tumor_name]][which((ddd>(1))&(ddd2>10)),1])#p.value of diff<0.001, log(fc)>0.5, and the gene show significant up regulation are larger than 10
  return(up_c)
}



##################################################################
#test


###
refesh_cycle_flux_list<-function(tumor_name, cycle_flux_list, cycle_edge_flux_list){
  print(tumor_name)
  
  #有gap现象的环的cycid
  gap_c = get_gap_c_cycid(tumor_name, cycle_edge_flux_list)
  #有up现象的环的cycid
  up_c = get_up_c_cycid(tumor_name, cycle_edge_flux_list)
  
  #################################
  #updownsign
  temp_tumor_df = cycle_edge_flux_list[[tumor_name]]
  temp_tumor_updownsign_df = temp_tumor_df[,c("cycid","foldchangesign")]
  #分类求和函数, 相同cycid的求和, 并mean()取均值
  #updownsign_df_sum = aggregate(temp_tumor_updownsign_df$foldchangesign, list(temp_tumor_updownsign_df$cycid), sum)
  updownsign_df_mean = aggregate(temp_tumor_updownsign_df$foldchangesign, list(temp_tumor_updownsign_df$cycid), mean)
  colnames(updownsign_df_mean) = c("cycid", "updownsign")
  cycle_flux_list[[tumor_name]][, "updownsign"] = updownsign_df_mean$updownsign
  
  
  #################################
  #ifup
  cycle_flux_list[[tumor_name]][, "ifup"] = "normal"
  #ifgap
  cycle_flux_list[[tumor_name]][, "ifgap"] = "normal"
  #若既是up又是gap,则视为gap
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% up_c), "ifup"] = "up"
  cycle_flux_list[[tumor_name]][which(cycle_flux_list[[tumor_name]]$cycle_id %in% gap_c), "ifgap"] = "gap"


  ##调整顺序 原 8:ifup 9:ifgap 10:updownsign 移动到前面
  #cycle_flux_list[[tumor_name]] = cycle_flux_list[[tumor_name]][,c(1,2,3,8,9,10,4:7)]
  
  return(cycle_flux_list)
}


##2.每一列为1个环 (cycle_flux_list)
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_final/cycle_flux_list.RData")
cycle_flux_list = cycle_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")
cycle_edge_flux_list = cycle_edge_flux_list

tumors_array = c("COAD", "ESCA","BLCA","BRCA","CESC","LUAD","LUSC","PRAD","STAD","THCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC")
for (i in 1:length(tumors_array)) {
  tumor_name = tumors_array[i]
  cycle_flux_list = refesh_cycle_flux_list(tumor_name, cycle_flux_list, cycle_edge_flux_list)
}


save(cycle_flux_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_final/cycle_flux_list.RData")







