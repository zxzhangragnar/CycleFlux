
########################################### 1.part_net #################################################
#######################################################################################################

load("E:/scFEA_universal/Data/zc513/Module_gene_map_SCMP_v1.1_functions.RData")
rm(list=setdiff(ls(), "F1_template_hsa_new_update"))
part_net = as.data.frame(F1_template_hsa_new_update)


# load("E:/scFEA_universal/my_R/aimA/original_InputData/input_data/KEGG_functions_FEA.RData")
# rm(list=setdiff(ls(), "F1_template_hsa"))
# part_net = F1_template_hsa


######################################### 2.merge_net_pathway ###########################################
#######################################################################################################
#此步不执行

######################################### 3.split_net_cincout ###########################################
#######################################################################################################
# 将net中有多cin多cout的一行拆成多行
# 注意：拆分后Rid将不唯一
#################################################
# 如 R00001
# C_in:C00404;C00001
# C_out:有一个
#
# 解决方案：
# 分成2行数据
# 保证每行只有一个C_in   （每行的Rid都是R00001）
# 
# ----------
# 如 R00002
# C_in 和 C_out 都有多个值
# 
# 解决方案：
# 分成多行n*m行数据(排列组合)   
# （每行的Rid都是R00002）
#################################################

#拆分多个c_in
library(tidyr)
hsa_net <- part_net %>% separate_rows(C_in, sep = ";")
hsa_net <- hsa_net %>% separate_rows(C_out, sep = ";")
rm(part_net)


#################################### 4.union_compound_tobe_rm.R ########################################
#######################################################################################################
library(readr)
CoFactorCID_d <- read_csv("E:/scFEA_universal/my_R/aimA/original_InputData/compound_to_be_delete/CoFactorCID_d.csv")
Compound_info_status_1_ <- read_csv("E:/scFEA_universal/my_R/aimA/original_InputData/compound_to_be_delete/Compound_info_status(1).csv")
never_considered_comp <- read_csv("E:/scFEA_universal/my_R/aimA/original_InputData/compound_to_be_delete/never_considered_comp.csv")

#生成：要从大网络中删去的全部化合物union
#来自zc628的三个csv文件

# union_compound_tobe_rm <- function() {
# 
#   #构建完整网络
#   part_com_status_1_ = Compound_info_status_1_
#   part_cofactor = CoFactorCID_d
#   part_never_con = never_considered_comp
#   
#   #308.正解一：C16254是tca主要环的一个必要节点，615行
#   part_never_con = part_never_con[-615,]
#   
#   #筛选符合条件的行（满足Consider_Stat==0的行）
#   part_com_status_0 = part_com_status_1_[which(part_com_status_1_$Consider_Stat==0),]
#   
#   rm(part_com_status_1_)
#   
#   #取3个part的compound集合的cid这列的值的并集 union
#   part_com_status_0_cid = as.list(part_com_status_0$CID)
#   part_cofactor_cid = as.list(part_cofactor$x)
#   part_never_con_cid = as.list(part_never_con$x)
#   
#   rm_comp_union = union(part_com_status_0_cid,part_cofactor_cid)
#   rm_comp_union = union(rm_comp_union,part_never_con_cid)
#   
#   rm_comp_union = unlist(rm_comp_union)
#   
#   return(rm_comp_union)
# }


union_compound_tobe_rm <- function() {

  #构建完整网络
  part_com_status_1_ = Compound_info_status_1_
  part_cofactor = CoFactorCID_d
  part_never_con = never_considered_comp

  #保留全部以C开头的化合物,删去全部以G开头的化合物
  part_never_con = part_never_con[-which(substr(part_never_con$x,1,1) == "C"),]
  

  #筛选符合条件的行（满足Consider_Stat==0的行）
  part_com_status_0 = part_com_status_1_[which(part_com_status_1_$Consider_Stat==0),]

  rm(part_com_status_1_)

  #取3个part的compound集合的cid这列的值的并集 union
  part_com_status_0_cid = as.list(part_com_status_0$CID)
  part_cofactor_cid = as.list(part_cofactor$x)
  part_never_con_cid = as.list(part_never_con$x)

  rm_comp_union = union(part_com_status_0_cid,part_cofactor_cid)
  rm_comp_union = union(rm_comp_union,part_never_con_cid)

  rm_comp_union = unlist(rm_comp_union)

  return(rm_comp_union)
}



#生成需要去掉的compound的并集
rm_comp_union = union_compound_tobe_rm()


#删去原始数据
rm(list=setdiff(ls(), c("hsa_net","rm_comp_union")))

#################################### 5.rm_com_from_net.R ########################################
#######################################################################################################
#注:union_compound_tobe_rm 在 union_compound_tobe_rm.R中
#先载入要删掉的化合物集合
# rm_comp_union = union_compound_tobe_rm()
# 
# #载入整理好的part_net
# part_net = get_part_net()
# part_net = merge_net_pathway()
# hsa_net = split_net_cin()
# hsa_net = split_net_cout()

###########################################################
#找到在hsa_net中要删掉的行的序号
#先遍历cin
#记录在for循环之后要删去的行
find_del_comp_index <- function() {
  row_tobe_del_index = c()
  
  for(i in 1:nrow(hsa_net)) {
    if(hsa_net[i,]$C_in %in% rm_comp_union) {
      row_tobe_del_index = append(row_tobe_del_index,i)
    }
  }
  
  #再遍历cout
  for(i in 1:nrow(hsa_net)) {
    if(hsa_net[i,]$C_out %in% rm_comp_union) {
      row_tobe_del_index = append(row_tobe_del_index,i)
    }
  }
  row_tobe_del_index = unique(row_tobe_del_index)
  row_tobe_del_index = sort(row_tobe_del_index)
  return(row_tobe_del_index)
}

row_tobe_del_index = find_del_comp_index()
###########################################################
#从hsa_net中删去这些行
# 这是最终结果
#part_net_after_delet = hsa_net[-row_tobe_del_index,]
hsa_net = hsa_net[-row_tobe_del_index,]

rm(list=setdiff(ls(), "hsa_net"))



# save 
# save(hsa_net, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData")
save(hsa_net, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData")








