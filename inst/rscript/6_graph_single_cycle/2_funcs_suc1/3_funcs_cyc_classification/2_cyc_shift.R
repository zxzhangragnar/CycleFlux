#
# 3.
# 某个地方 环的loss
# 环的shift: 1.断 2.新边
# 
# 写成的格式   直接给出更深的结果(那些边)
# 
# 4. formula
# 原来环 -> 哪些断了 -> 留下的东西 -> 新的边
# 
# 5.                                       ok
# 一共有多少环
# 符合上面类型 (shift) (即另辟蹊径)
# 给出statistical统计结果
# 
# 6.
# 变的意义 —— 生物学意义
#
#







##################### 3.另辟蹊径 #############################################
get_cyc_shift_formula <- function(tumor_name, cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_chain_list = gapup_cycle_chain_list[[tumor_name]]
  ##
  # result
  new_old_path_df = data.frame()

  
  for (i in 1:length(ug_c)) {
    cycid = ug_c[i]
    print(cycid)
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycid == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycid == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycid == ug_c[i]),]
    
    up_indeg = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),]
    up_outdeg = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),]
    
    up_indeg_arr = up_indeg[,"ind"]
    up_outdeg_arr = up_outdeg[,"od"]
    
    ist_deg_arr = intersect(up_indeg_arr, up_outdeg_arr)
    if(length(ist_deg_arr) != 0) {
      #shared_deg_arr = c()
      for (j in 1:length(ist_deg_arr)) {
        indegnode = up_indeg[which(up_indeg$ind == ist_deg_arr[j]), "node"] #hnode
        outdegnode = up_outdeg[which(up_outdeg$od == ist_deg_arr[j]), "node"] #tnode
        #老路的全部node节点即为 留下的东西(留下的node节点)
        #观察 留下的node节点 中 哪些edge断了gap 
        
        cond1 = !identical(indegnode, outdegnode)
        cond2 = (identical(indegnode, outdegnode) & !((length(indegnode)==1) & (length(outdegnode)==1)))
        if ( cond1 | cond2 ) {
          permutation_node = data.frame()
          for (od in outdegnode) {
            for (ind in indegnode) {
              temp_node_df = data.frame()
              temp_node_df[1, "od"] = od
              temp_node_df[1, "ind"] = ind
              
              permutation_node = rbind(permutation_node, temp_node_df)
            }
          }
          
          if (length(which(permutation_node$od == permutation_node$ind)) > 0) {
            permutation_node = permutation_node[-which(permutation_node$od == permutation_node$ind),]
          }
          
          
          # 1.[新路]另辟的这条道路上 只有1个节点ist_deg_arr
          new_path_str = paste0(permutation_node[1,"od"], " -> U -> ", ist_deg_arr[j], " -> U -> ", permutation_node[1,"ind"])
          
          # 2.[老路]原道路的 起始点indegnode和 终止点outdegnode
          # 哪些断了 -> 留下的东西
          #在 gapup_cycle_chain_list 中寻找
          
          old_path_arr = c()
          for (chn in ug_chain_list[[i]]) {
            temp_chain_str = chn
            temp_chain_arr = unlist(strsplit(temp_chain_str,split = " -> "))
            
            od_index = which(temp_chain_arr == permutation_node[1,"od"])
            ind_index = which(temp_chain_arr == permutation_node[1,"ind"])
            length(od_index)
            length(ind_index)
            
            if(length(od_index) > 1) {
              od_index = od_index[1]
            }
            
            if(length(ind_index) > 1) {
              ind_index = ind_index[1]
            }
            
            if(od_index > ind_index) {
              index_1 = ind_index
              index_2 = od_index
            }else {
              index_1 = od_index
              index_2 = ind_index
            }
            
            
            
            temp_old_path_arr_1 = temp_chain_arr[index_1:index_2]
            if("G" %in% temp_old_path_arr_1) {
              temp_old_path_str_1 = paste(temp_old_path_arr_1, collapse  = " -> ")
              old_path_arr = append(old_path_arr, temp_old_path_str_1)
            }
            
            #逆向的环
            chain_part_1 = temp_chain_arr[1:index_1]
            chain_part_2 = temp_chain_arr[index_2:(length(temp_chain_arr)-1)]
            temp_old_path_arr_2 = c(chain_part_2, chain_part_1)
            
            if("G" %in% temp_old_path_arr_2) {
              temp_old_path_str_2 = paste(temp_old_path_arr_2, collapse  = " -> ")
              old_path_arr = append(old_path_arr, temp_old_path_str_2)
            }
          }
          
          old_path_arr = unique(old_path_arr)
          
          #找老路的全部node节点即可
          old_path_str = paste(old_path_arr, collapse  = " ; ")
          
          
          ## 为df添加行
          temp_old_path_df = data.frame()
          temp_old_path_df[1, "cycid"] = cycid
          temp_old_path_df[1, "new_path"] = new_path_str
          temp_old_path_df[1, "old_path"] = old_path_str
          temp_old_path_df[1, "endpoint_node"] = paste0(permutation_node[1,"od"], ";", permutation_node[1,"ind"])
          temp_old_path_df[1, "shift_node"] = ist_deg_arr[j]
          
          new_old_path_df = rbind(new_old_path_df, temp_old_path_df)
          
          
        }
      }
    }

  }
  
  return(new_old_path_df)
}



#################################################################################
cyc_shift_main <- function(res_path, input_tumor_name) {
  
  # init
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/3_funcs_cyc_classification/cycle_upgap_class_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_in.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data/cycle_edge_flux_list_out.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))

  
  ##
  cycle_shift_path_df_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    new_old_path_df = get_cyc_shift_formula(tumor_name, cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    cycle_shift_path_df_list[[tumor_name]] = new_old_path_df
  }
  
  # save
  res_sub_path = "6_graph_single_cycle/2_funcs_suc1/3_funcs_cyc_classification"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_shift_path_df_list.RData")
  
  save(cycle_shift_path_df_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/2_funcs_suc1/3_funcs_cyc_classification/cycle_shift_path_df_list.RData")
  
}





# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# cyc_shift_main(res_path, input_tumor_name)




















