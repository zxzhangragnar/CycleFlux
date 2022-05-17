
##1.
get_hsa_net <- function(hsa_net, pathway_name) {
  hsa_net = hsa_net[which(hsa_net$Pathway == pathway_name),]
  rownames(hsa_net) = c(1:length(rownames(hsa_net)))
  return(hsa_net)
}



##2.此部分在python 中处理
# get_all_pathway_partnet_rowidx_df <- function(hsa_net) {
#   #1.得到hsa_net中的全部pathway种类并去重 并删去pathway为空""的值
#   pathway_arr = as.character(hsa_net$Pathway)
#   pathway_arr = unique(pathway_arr)
#   #pathway_arr = pathway_arr[-which(pathway_arr == "")] #83个
#   #2.标记每个pathway对应hsa_net行号 然后在python中再提取
#   all_pathway_partnet_rowidx_df = data.frame()
#   
#   for(i in 1:length(pathway_arr)) {
#     spnidx = which(grepl(pathway_arr[i], hsa_net$Pathway))
#     this_pth_str = paste(spnidx,collapse=";")
#     all_pathway_partnet_rowidx_df[i,1] = this_pth_str
#   }
#   #为pathway添加行名
#   rownames(all_pathway_partnet_rowidx_df) = pathway_arr
#   return(all_pathway_partnet_rowidx_df)
# }


load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData")

##1.
hsa_net = get_hsa_net(hsa_net, "hsa01100")
#hsa_net = get_hsa_net(hsa_net, "hsa00350")
#hsa_net = get_hsa_net(hsa_net, "hsa00280")
#hsa_net = get_hsa_net(hsa_net, "hsa00260")




###############################################################################
# save
save(hsa_net, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData")






################################################################################
### debug 用于测试时选择 某个子pathway
# pathwaylength_df = all_pathway_partnet_rowidx_df
# for(i in 1:length(all_pathway_partnet_rowidx_df[,1])) {
#   tmp_list = unlist(strsplit(all_pathway_partnet_rowidx_df[i,1],split = ";"))
#   pathwaylength_df[i,"length"] = length(tmp_list)
# }


