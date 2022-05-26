# 
# "COAD" 
# ## each edge
# 
# 1.Up仅头入尾出：7,58,60,106,149,484,211,213,261,265,268,272,278,376
# 2.Up中间有up出度: 405,406,414,278,500, 336(M),348(M)
# 3.Gap中间有up的出度:426,430,469,470,472,476,497,469,225,227
# 4.Up另辟蹊径up:36,37,473,478,123,311(M),319(M), 320(M),325(M),348,487(M)
# 
# 这4种情况分别代表了什么样的生物现象(生物故事)?
# 


####################
# 1.Up仅头入尾出
#
# # 1.Up仅头入尾出:
# 环中仅2个节点node有 up的度
# 且这2个节点，在这个环中作为c_in或c_out时，ifup==up

####################
# 2.Up中间有up出度:
#   环中有>3个节点node有 up的度
# 且这些节点，在这个环中作为c_in或c_out时，两端的ifgap不全为gap


####################
# 3.另辟蹊径:
#   存在一个度节点（ind或od）
# 这个节点同时是这个环中2个节点的度
# 且这2个节点，在这个环中作为c_in或c_out时，ifup==up
# 且这2个节点之间的边，存在gap的边




#################################################################################
#################################################################################

# tumor_name = "COAD"
# indeg_df = cycle_edge_flux_list_in[[tumor_name]]
# outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
# cyc_df = cycle_edge_flux_list[[tumor_name]]
# ug_c = names(gapup_cycle_chain_list[[tumor_name]])


##################################################################################
##################### 1.Up仅头入尾出 #############################################
get_class_1_c <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ##
  class_1_c = c()
  for (i in 1:length(ug_c)) {
    class_1_judge = TRUE
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycle_id == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycle_id == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycle_id == ug_c[i]),]
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)
    
    if (length(up_degnode_arr) == 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(!("up" %in% tmp_node_row[,"ifup"])){
          class_1_judge = FALSE
        }
        
      }
    }else {
      class_1_judge = FALSE
    }
    
    if (class_1_judge) {
      class_1_c = append(class_1_c, ug_c[i])
    }
    
  }
  return(class_1_c)
  
}



##################################################################################
##################### 2.Up中间有up出度 #############################################
get_class_2_c <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ##
  class_2_c = c()
  for (i in 1:length(ug_c)) {
    class_2_judge = TRUE
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycle_id == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycle_id == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycle_id == ug_c[i]),]
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)
    
    if (length(up_degnode_arr) > 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(("gap" %in% tmp_node_row[,"ifgap"]) & !("up" %in% tmp_node_row[,"ifup"])){
          class_2_judge = FALSE
        }
        
      }
    }else {
      class_2_judge = FALSE
    }
    
    if (class_2_judge) {
      class_2_c = append(class_2_c, ug_c[i])
    }
    
  }
  return(class_2_c)
  
}




##################################################################################
##################### 3.另辟蹊径 #############################################
get_class_3_c <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ##
  class_3_c = c()
  for (i in 1:length(ug_c)) {
    class_3_judge = FALSE
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycle_id == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycle_id == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycle_id == ug_c[i]),]
    
    up_indeg = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),]
    up_outdeg = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),]
    
    up_indeg_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"ind"]
    up_outdeg_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"od"]
    
    ist_deg_arr = intersect(up_indeg_arr, up_outdeg_arr)
    if(length(ist_deg_arr) != 0) {
      #shared_deg_arr = c()
      for (j in 1:length(ist_deg_arr)) {
        indegnode = up_indeg[which(up_indeg$ind == ist_deg_arr[j]), "node"] #hnode
        outdegnode = up_outdeg[which(up_outdeg$od == ist_deg_arr[j]), "node"] #tnode
        if(length(indegnode)>0 & length(outdegnode)>0){
          cond1 = !identical(indegnode, outdegnode)
          cond2 = (identical(indegnode, outdegnode) & !((length(indegnode)==1) & (length(outdegnode)==1)))
          if ( cond1 | cond2 ) {
            class_3_judge = TRUE
            #shared_deg_arr = append(shared_deg_arr, ist_deg_arr[j])
          } 
        }
        
      }    
    }
    
    if (class_3_judge) {
      class_3_c = append(class_3_c, ug_c[i])
    }
    
  }
  
  return(class_3_c)
}






##################################################################################
##################### 4.从1点爆发 #############################################
get_class_4_c <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ##
  class_4_c = c()
  for (i in 1:length(ug_c)) {
    class_4_judge = TRUE
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycle_id == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycle_id == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycle_id == ug_c[i]),]
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)
    
    if (length(up_degnode_arr) == 1) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(!("up" %in% tmp_node_row[,"ifup"])){
          class_4_judge = FALSE
        }
        
      }
    }else {
      class_4_judge = FALSE
    }
    
    if (class_4_judge) {
      class_4_c = append(class_4_c, ug_c[i])
    }
    
  }
  return(class_4_c)
  
}


##################################################################################
##################### 5.其他类(Gap中间有up的出度) ###################################################
get_ug_class_df <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  cyc_df = as.data.frame(cycle_edge_flux_list[[tumor_name]])
  
  #分类结果
  ug_class_df = as.data.frame(unique(cyc_df[which(cyc_df$cycle_id %in% ug_c), "cycle_id"]))
  ug_class_df[,"classid"] = 0
  colnames(ug_class_df) = c("cycle_id", "classid")
  
  class_1_c = get_class_1_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_2_c = get_class_2_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_3_c = get_class_3_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_4_c = get_class_4_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  
  ## 属于第3类的环, 从第1类和第2类中去掉
  class_1_c = setdiff(class_1_c, class_3_c)
  class_2_c = setdiff(class_2_c, class_3_c)
  
  class_ist = c(class_1_c, class_2_c, class_3_c, class_4_c)
  class_ist = unique(class_ist)
  
  ## 5.其他类(Gap中间有up的出度)
  class_5_c = setdiff(ug_c, class_ist)
  
  ug_class_df[which(ug_class_df$cycle_id %in% class_1_c), "classid"] = 1
  ug_class_df[which(ug_class_df$cycle_id %in% class_2_c), "classid"] = 2
  ug_class_df[which(ug_class_df$cycle_id %in% class_3_c), "classid"] = 3
  ug_class_df[which(ug_class_df$cycle_id %in% class_4_c), "classid"] = 4
  ug_class_df[which(ug_class_df$cycle_id %in% class_5_c), "classid"] = 5
  
  return(ug_class_df)
}





#################################################################################
cyc_classification_main <- function(res_path, input_tumor_name) {
  
  # init
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_in.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_out.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  
  ##
  cycle_upgap_class_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    ug_class_df = get_ug_class_df(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) 
    cycle_upgap_class_list[[tumor_name]] = ug_class_df  
  }
  
  # save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_upgap_class_list.RData")
  
  save(cycle_upgap_class_list, file=res_file_path)
}




# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# cyc_classification_main(res_path, input_tumor_name)




