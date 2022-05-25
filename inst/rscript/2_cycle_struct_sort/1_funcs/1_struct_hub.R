# 
# 1.找到这些枢纽节点 这些节点同时属于很多环
# 



##################################### part1 找到枢纽节点 这些节点同时属于很多环 #################################

get_compound_directed <- function(cycle_directed) {
  compound_directed = data.frame()
  for (i in 1:length(cycle_directed[,1])) {
    cycle_id = cycle_directed[i, "cycle_id"]
    compound_chain = cycle_directed[i, "compound_chain"]
    
    compound_chain = unlist(strsplit(compound_chain, split = ";"))
    cycle_node = unique(unlist(strsplit(compound_chain[1], split = "->")))
    
    for (j in 1:length(cycle_node)) {
      cpdname = cycle_node[j]
      temp_row = data.frame()
      temp_row[1, "cpdname"] = cpdname
      temp_row[1, "cycle"] = cycle_id
      
      compound_directed = rbind(compound_directed, temp_row)
    }
  }
  
  return(compound_directed)
}


get_cpd_hubnode <- function(frequency, cycle_directed) {
  compound_directed = get_compound_directed(cycle_directed)
  
  cpd_table = as.data.frame(table(compound_directed[,c("cpdname")]))
  colnames(cpd_table) = c("cpd","fre")
  #function in function
  merge_cpd_cycle <- function(compound_directed) {
    cpd_part_df = compound_directed[,c("cpdname","cycle")]
    library(dplyr)
    cpd_part_df = unique(cpd_part_df)
    cpd_part_df = aggregate(cpd_part_df, list(cpd_part_df$cpdname), paste, collapse = " ; ")
    cpd_part_df = cpd_part_df[,-2]
    colnames(cpd_part_df)[1] <- 'cpdname' 
    return(cpd_part_df)
  }
  
  cpd_part_df = merge_cpd_cycle(compound_directed)
  cpd_cyc_fre = cbind(cpd_table, cpd_part_df[,"cycle"])
  colnames(cpd_cyc_fre) = c("cpdname", "fre", "cycle")
  

  cpd_hubnode = cpd_cyc_fre[which(cpd_cyc_fre$fre >= frequency),]

  
  return(cpd_hubnode)
}





##################################### part2 去掉枢纽节点h后，判断这些环是否仍然相连？ #################################
# 如何判断这些环是否仍然相连？
# 对环的shared_cycle内的所有环逐一判断，是否仍然相连？
# 若仍然都相连，则h枢纽性不大
# 若很多都不再相连了，则h枢纽性大
# 判断2个环是否相连：看 2个环的compound是否有交集，若有交集 即这2个环仍然相连 即满足shared_cycle

get_cyc_shared_cyc <- function(cycle_directed) {
  cyc_shared_cyc = cycle_directed[,c("cycle_id", "cpds", "shared_cycle", "shared_cycles")]
  cyc_shared_cyc[, "old_hubdeg"] = cyc_shared_cyc[, "shared_cycle"]
  #wash "cpds"
  for (i in 1:length(rownames(cyc_shared_cyc))) {
    temp_cyc_compounds = cyc_shared_cyc[i, "cpds"]
    temp_cyc_compounds = substring(temp_cyc_compounds, 2, nchar(temp_cyc_compounds)-1)
    temp_cyc_compounds <- unlist(strsplit(temp_cyc_compounds, split = ", "))
    for (j in 1:length(temp_cyc_compounds)) {
      temp_cyc_compounds[j] = substring(temp_cyc_compounds[j], 2, nchar(temp_cyc_compounds[j])-1)
    }
    temp_cyc_compounds = paste(temp_cyc_compounds[1:length(temp_cyc_compounds)], collapse = ";")
    cyc_shared_cyc[i, "cpds"] = temp_cyc_compounds
  }
  
  #wash shared_cycles
  for (i in 1:length(rownames(cyc_shared_cyc))) {
    temp_cyc_compounds = cyc_shared_cyc[i, "shared_cycles"]
    temp_cyc_compounds = substring(temp_cyc_compounds, 2, nchar(temp_cyc_compounds)-1)
    temp_cyc_compounds <- unlist(strsplit(temp_cyc_compounds, split = ", "))
    temp_cyc_compounds = paste(temp_cyc_compounds[1:length(temp_cyc_compounds)], collapse = ";")
    cyc_shared_cyc[i, "shared_cycles"] = temp_cyc_compounds
  }
  cyc_shared_cyc[, "old_hubdegs"] = cyc_shared_cyc[, "shared_cycles"]
  return(cyc_shared_cyc)
}





###
#remove a compound from all the cycles
remove_cpd_in_cyc_shared_cyc <- function(cpd_to_be_rm, cyc_shared_cyc) {
  #cpd_to_be_rm = "C05345"
  for (i in 1:length(rownames(cyc_shared_cyc))) {
    temp_cyc_compounds = cyc_shared_cyc[i, "cpds"]
    temp_cyc_compounds <- unlist(strsplit(temp_cyc_compounds, split = ";"))
    after_rm_temp_cyc_compounds = c()
    for (j in 1:length(temp_cyc_compounds)) {
      if(temp_cyc_compounds[j] == cpd_to_be_rm) {
        #pass
      }else {
        after_rm_temp_cyc_compounds = append(after_rm_temp_cyc_compounds, temp_cyc_compounds[j])
      } 
    }
    after_rm_temp_cyc_compounds = paste(after_rm_temp_cyc_compounds[1:length(after_rm_temp_cyc_compounds)], collapse = ";")
    cyc_shared_cyc[i, "aft_rm_cpds"] = after_rm_temp_cyc_compounds    
  }
  
  return(cyc_shared_cyc)
}






# if cycle_A and cycle_B have Intersection?
judge_cycle_have_intersection <- function(cyc_shared_cyc, cycle_id_A, cycle_id_B) {
  # A:line_1 B:line_2
  cycle_A_cpd = cyc_shared_cyc[which(cyc_shared_cyc$cycle_id == cycle_id_A), "aft_rm_cpds"]
  cycle_B_cpd = cyc_shared_cyc[which(cyc_shared_cyc$cycle_id == cycle_id_B), "aft_rm_cpds"]
  cycle_A_cpd <- unlist(strsplit(cycle_A_cpd, split = ";"))
  cycle_B_cpd <- unlist(strsplit(cycle_B_cpd, split = ";"))
  
  interaction_bet_cycles = intersect(cycle_A_cpd, cycle_B_cpd)
  has_interaction = FALSE
  if(length(interaction_bet_cycles) >= 1) {
    has_interaction = TRUE
  }
  return(has_interaction)  
}






# new hubdeg of cyc_shared_cyc
calculate_new_hubdeg <- function(cyc_shared_cyc) {
  for (i in 1:length(rownames(cyc_shared_cyc))) {
    #print(i)
    cycle_id_slf = cyc_shared_cyc[i, "cycle_id"]
    old_hubdegs = cyc_shared_cyc[i, "old_hubdegs"]
    old_hubdegs <- unlist(strsplit(old_hubdegs, split = ";"))
    new_hubdegs = c()
    for (j in 1:length(old_hubdegs)) {
      cycle_id_ngh = old_hubdegs[j]
      if_intersect = judge_cycle_have_intersection(cyc_shared_cyc, cycle_id_slf, cycle_id_ngh)
      if(if_intersect) {
        new_hubdegs = append(new_hubdegs, cycle_id_ngh)
      }
    }
    cyc_shared_cyc[i, "new_hubdeg"] = length(new_hubdegs)
    new_hubdegs = paste(new_hubdegs[1:length(new_hubdegs)], collapse = ";")
    cyc_shared_cyc[i, "new_hubdegs"] = new_hubdegs
  }
  return(cyc_shared_cyc)
}

# wilcox_test
deg_wilcox_test <- function(data_new, data_old) {
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test
  wilcox_res = wilcox.test(data_new, data_old, exact = FALSE, alternative = "less")
  # if alternative = "greater" pvalue = 1
  pvalue = wilcox_res[["p.value"]]
  
  if(pvalue < 0.05){
    if_obvious_different = TRUE
  }else {
    if_obvious_different = FALSE
  }
  
  res_if_obvious_different = list()
  res_if_obvious_different[["pvalue"]] = pvalue
  res_if_obvious_different[["judge"]] = if_obvious_different
  return(res_if_obvious_different)
}


#main func
get_cpd_to_be_rm_judge_list_tool <- function(cpd_to_be_rm, cpd_hubnode, cycle_directed) {
  cyc_shared_cyc = get_cyc_shared_cyc(cycle_directed)
  cyc_shared_cyc = remove_cpd_in_cyc_shared_cyc(cpd_to_be_rm, cyc_shared_cyc)
  cyc_shared_cyc = calculate_new_hubdeg(cyc_shared_cyc)
  #compare
  #cyc_shared_cyc[, "old_hubdegs"] == cyc_shared_cyc[, "new_hubdegs"]
  
  ### 看看这些环在去掉枢纽节点后，连通性是否变化很大
  temp_hubcpd_cycs = cpd_hubnode[which(cpd_hubnode$cpdname == cpd_to_be_rm), "cycle"]
  temp_hubcpd_cycs <- unlist(strsplit(temp_hubcpd_cycs, split = " ; "))
  #去掉枢纽节点后，连通性对比
  temp_hubcpd_changed_compared = cyc_shared_cyc[which(cyc_shared_cyc$cycle_id %in% temp_hubcpd_cycs), c("cycle_id", "old_hubdeg", "new_hubdeg")]
  cyc_shared_df = cyc_shared_cyc[which(cyc_shared_cyc$cycle_id %in% temp_hubcpd_cycs), ]
  
  # hub_node: 50% cycle changed amplitude > 50%
  # test diff
  # using wilcox.test
  data_old = temp_hubcpd_changed_compared[,"old_hubdeg"]
  data_new = temp_hubcpd_changed_compared[,"new_hubdeg"]
  res_if_obvious_different = deg_wilcox_test(data_new, data_old) 
  
  temp_cpd_to_be_rm_judge_list = list()
  temp_cpd_to_be_rm_judge_list[["cyc_shared_df"]] = as.data.frame(cyc_shared_df)
  temp_cpd_to_be_rm_judge_list[["compared"]] = as.data.frame(temp_hubcpd_changed_compared)
  temp_cpd_to_be_rm_judge_list[["judge"]] = res_if_obvious_different
  
  return(temp_cpd_to_be_rm_judge_list)
}

#final result
get_cpd_to_be_rm_judge_list <- function(cpd_hubnode, cycle_directed) {
  cpd_to_be_rm_judge_list = list()
  if(length(rownames(cpd_hubnode)) > 0) {
    for (c in 1:length(rownames(cpd_hubnode))) {
      print(c)
      cpd_to_be_rm = as.character(cpd_hubnode[c,"cpdname"])
      temp_list = get_cpd_to_be_rm_judge_list_tool(cpd_to_be_rm, cpd_hubnode, cycle_directed)
      cpd_to_be_rm_judge_list[[cpd_to_be_rm]] = temp_list
    }    
  }
  return(cpd_to_be_rm_judge_list)
}




############################################################
#枢纽节点h 所属的环中 哪些环属于 hub-linking结构？
#对于cpd_hubnode某个cpd所属的某个环cycn 其 old_hubdeg/new_hubdeg > 3 
#则这个环为 这个枢纽节点h的 hub-linking 环

############################################################
#对于 枢纽节点 h1和h2其所属的环 有很多重合的部分 怎么算？

##
#res
get_hubcpd_cycno_list_df <- function(cpd_to_be_rm_judge_list, fc_threshold) {
  library(gtools)
  hubcpd_cycno_list = list()
  hubcpd_cycno_list_df = data.frame()
  hubcpd_names = names(cpd_to_be_rm_judge_list)
  for (i in 1:length(cpd_to_be_rm_judge_list)) {
    temp_hubcpd = hubcpd_names[i] #"C00019"
    temp_hubcpd_compared_df = cpd_to_be_rm_judge_list[[temp_hubcpd]][["compared"]]
    temp_part_df = data.frame()
    for (j in 1:length(rownames(temp_hubcpd_compared_df))) {
      old_hubdeg = temp_hubcpd_compared_df[j, "old_hubdeg"]
      new_hubdeg = temp_hubcpd_compared_df[j, "new_hubdeg"]
      fc <- foldchange(old_hubdeg, new_hubdeg)
      if(fc > fc_threshold) {
        temp_hubcpd_compared_df[j, "hubcpd"] = temp_hubcpd
        temp_part_df = rbind(temp_part_df, temp_hubcpd_compared_df[j, ])
      }
    }
    hubcpd_cycno_list[[temp_hubcpd ]] = temp_part_df
    hubcpd_cycno_list_df = rbind(hubcpd_cycno_list_df, temp_part_df)
  }  
  
  return(hubcpd_cycno_list_df)
}

###
# 一个枢纽节点h 最终确定的属于hub-linking结构的环的数量一定要大于 10 
# 否则 这个枢纽节点和其所确定的环 无效
filter_hub_num_under10 <- function(cpd_to_be_rm_judge_list, frequency, fc_threshold) {
  hub_struct_cycle_id = c()
  if (length(cpd_to_be_rm_judge_list) > 0) {
    hubcpd_cycno_list_df = get_hubcpd_cycno_list_df(cpd_to_be_rm_judge_list, fc_threshold)
    
    hub_cycle_ids = hubcpd_cycno_list_df[,"cycle_id"]
    #final result: hub linking struct
    hub_cycle_ids = unique(hub_cycle_ids)
    #filter
    hub_node_belongs_num = as.list(table(hubcpd_cycno_list_df[,"hubcpd"]))
    for (hn in names(hub_node_belongs_num)) {
      if(hub_node_belongs_num[[hn]] < frequency) {
        #remove
        hubcpd_cycno_list_df = hubcpd_cycno_list_df[-which(hubcpd_cycno_list_df$hubcpd == hn),]
      }
    }
    hub_struct_cycle_id = hubcpd_cycno_list_df[, "cycle_id"]
  }

  return(hub_struct_cycle_id)
}




#############################################################################################
# test
struct_hub_main <- function(output_path, res_path, prm_1=5, prm_2=5, prm_3=3) {
  
  load(file.path(output_path, "cycle_directed.RData"))
  
  ##可疑hub-node
  cpd_hubnode = get_cpd_hubnode(prm_1, cycle_directed)
  
  cpd_to_be_rm_judge_list = get_cpd_to_be_rm_judge_list(cpd_hubnode, cycle_directed)
  
  hub_struct_cycle_id = filter_hub_num_under10(cpd_to_be_rm_judge_list, prm_2, prm_3)
  
  hub_struct_cycle_id = as.numeric(hub_struct_cycle_id)
  
  #save
  res_sub_path = "2_cycle_struct_sort/result_struct"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "hub_struct_cycle_id.RData")
  
  save(hub_struct_cycle_id, file=res_file_path)
}



#######################################################################################################################                 
# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# #struct_hub_main(10, 10, 3)
# struct_hub_main(output_path, res_path, 5, 5, 3)
