


#情况3: 合并情况1和情况2 (_e1)
#e1<--up--(c2<--up--c1)<--up--e1
###################################################################

# st2_cycid = select_st2_1_gapinfo_e1e2$cycid
# st3_cycid = select_st3_1_gapinfo_e1e2$cycid
# intersect_cycid = intersect(st2_cycid, st3_cycid)

## 不要求 某个环同时有 od & ind 才计算在内
## 而是 做并集 只要这个环 有 od | ind 就计算在内 



###########################################################################
###########################################################################
# 1.只包括 up

plot_up_cycle<-function(tumor_name, plot_name, select_in_upinfo, select_out_upinfo, gapup_cycle_chain_list, never_considered_comp_names) {
  tmp_chain_str = gapup_cycle_chain_list[[tumor_name]]
  if (length(tmp_chain_str) > 0) {
    cyc_id_arr = names(tmp_chain_str)
    
    for (k in 1:length(cyc_id_arr)) {
      tmp_cyc_id = cyc_id_arr[k]
      tmp_cycle_chain_str = tmp_chain_str[[tmp_cyc_id]] #[1]C->U->c->G->c [2]c->G->c->U->c
      tmp_cycle_chain_arr = unlist(strsplit(tmp_cycle_chain_str[1],split = " -> ")) #'c' 'U' 'c' 'G' 'c'
      
      ## 添加环中的节点
      cyc_node_vector = c()
      for (i in 1:length(tmp_cycle_chain_arr)) {
        if((i%%2 == 1) & (i != length(tmp_cycle_chain_arr))) {
          cyc_node_vector = append(cyc_node_vector, tmp_cycle_chain_arr[i])
        }
      }
      
      
      ##
      ## successors部分
      #第k个环
      #有入度的结点
      up_out_node = select_out_upinfo[which(select_out_upinfo$cycid == tmp_cyc_id), "node"]
      up_out_od = select_out_upinfo[which(select_out_upinfo$cycid == tmp_cyc_id), "od"]
      #有出度的结点
      up_in_node = select_in_upinfo[which(select_in_upinfo$cycid == tmp_cyc_id), "node"]
      up_in_ind = select_in_upinfo[which(select_in_upinfo$cycid == tmp_cyc_id), "ind"]
      
      #degree_node_arr = c(rbind(up_out_od, up_in_ind))
      degree_node_arr = c(up_out_od, up_in_ind)
      degree_node_arr = unique(degree_node_arr)
      
      ## successors部分
      library(igraph)
      all_node_num = length(cyc_node_vector) + length(degree_node_arr)
      all_node_vector = append(cyc_node_vector, degree_node_arr)
      
      g <- make_empty_graph(n = all_node_num)
      
      #分2段给graph中的vertex添加名字
      #1.给cyc上的vertex添加名字 2.给cyc的degree的vertex添加名字
      g <- set.vertex.attribute(g, "name", value=all_node_vector)
      
      
      #######################################################
      #标记 never_considered_comp 在图中这些化合物节点为grey色
      for (i in 1:length(all_node_vector)) {
        if(all_node_vector[i] %in% never_considered_comp_names) {
          V(g)[name==all_node_vector[i]]$color <- "grey"
        }
      }
      #######################################################
      
      ###########
      ## 添加边
      ## 新的建环方式
      for (i in 1:length(tmp_cycle_chain_str)) {
        chain_cycle_temp = tmp_cycle_chain_str[i]
        chain_cycle_temp_arr = unlist(strsplit(chain_cycle_temp,split = " -> "))
        
        for (j in 1:length(chain_cycle_temp_arr)) {
          if(j%%2 == 0) {
            edge_info = chain_cycle_temp_arr[j]
            tmp_cyc_hnode = V(g)[name==chain_cycle_temp_arr[j-1]]
            tmp_cyc_tnode = V(g)[name==chain_cycle_temp_arr[j+1]]
            if(edge_info == "G") {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "red")
            }else if(edge_info == "U") {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "green")
            }else {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "black")
            }
          }
        }
      }
      #######################################################
      
      ####################################################
      ## successors部分
      # 根据name查找某个点 V(g)[name=="C00370"]
      #g <- add_edges(g, c(V(g)[name=="C00370"],V(g)[name=="C01958"],V(g)[name=="C00028"]), color = "black")
      
      
      #有入度的结点
      # up_out_node
      # up_out_od
      #有出度的结点
      # up_in_node
      # up_in_ind
      
      if (!identical(up_out_node == 0, logical(0))) {
        for (i in 1:length(up_out_node)) {
          tmp_cyc_node = V(g)[name==up_out_node[i]]
          tmp_cyc_od = V(g)[name==up_out_od[i]]
          
          if (!are.connected(g, tmp_cyc_node, tmp_cyc_od)) {
            g <- add_edges(g, c(tmp_cyc_node, tmp_cyc_od), color = "cyan")  
          }
          
        }
      }
      
      
      if (!identical(up_in_node == 0, logical(0))) {
        for (i in 1:length(up_in_node)) {
          tmp_cyc_node = V(g)[name==up_in_node[i]]
          tmp_cyc_ind = V(g)[name==up_in_ind[i]]
          
          if (!are.connected(g, tmp_cyc_ind, tmp_cyc_node)) {
            g <- add_edges(g, c(tmp_cyc_ind, tmp_cyc_node), color = "cyan")
          }
          
        }
      }
      
      ## successors部分
      ####################################################
      
      tmp_graph_name = paste0(tumor_name, "_C", cyc_id_arr[k], ".png")
      tmp_plot_name = file.path(plot_name, tmp_graph_name)
      
      png(tmp_plot_name, 500, 500)
      plot(g)
      dev.off()
    }
  }

  
}


###########################################################################
###########################################################################
# 2.既包括 up 也包括 gap

plot_ug_cycle<-function(tumor_name, plot_name, select_in_upinfo, select_out_upinfo, select_in_gapinfo, select_out_gapinfo, gapup_cycle_chain_list, never_considered_comp_names) {
  tmp_chain_str = gapup_cycle_chain_list[[tumor_name]]
  if (length(tmp_chain_str) > 0) {
    cyc_id_arr = names(tmp_chain_str)
    
    for (k in 1:length(cyc_id_arr)) {
      tmp_cyc_id = cyc_id_arr[k]
      tmp_cycle_chain_str = tmp_chain_str[[tmp_cyc_id]] #[1]C->U->c->G->c [2]c->G->c->U->c
      tmp_cycle_chain_arr = unlist(strsplit(tmp_cycle_chain_str[1],split = " -> ")) #'c' 'U' 'c' 'G' 'c'
      
      ## 添加环中的节点
      cyc_node_vector = c()
      for (i in 1:length(tmp_cycle_chain_arr)) {
        if((i%%2 == 1) & (i != length(tmp_cycle_chain_arr))) {
          cyc_node_vector = append(cyc_node_vector, tmp_cycle_chain_arr[i])
        }
      }
      
      
      ##
      ## successors部分
      #第k个环
      #有入度的结点
      up_out_node = select_out_upinfo[which(select_out_upinfo$cycid == tmp_cyc_id), "node"]
      up_out_od = select_out_upinfo[which(select_out_upinfo$cycid == tmp_cyc_id), "od"]
      #有出度的结点
      up_in_node = select_in_upinfo[which(select_in_upinfo$cycid == tmp_cyc_id), "node"]
      up_in_ind = select_in_upinfo[which(select_in_upinfo$cycid == tmp_cyc_id), "ind"]
      
      #有入度的结点
      gap_out_node = select_out_gapinfo[which(select_out_gapinfo$cycid == tmp_cyc_id), "node"]
      gap_out_od = select_out_gapinfo[which(select_out_gapinfo$cycid == tmp_cyc_id), "od"]
      #有出度的结点
      gap_in_node = select_in_gapinfo[which(select_in_gapinfo$cycid == tmp_cyc_id), "node"]
      gap_in_ind = select_in_gapinfo[which(select_in_gapinfo$cycid == tmp_cyc_id), "ind"]
      
      
      #degree_node_arr = c(rbind(up_out_od, up_in_ind, gap_out_od, gap_in_ind))
      degree_node_arr = c(up_out_od, up_in_ind, gap_out_od, gap_in_ind)
      degree_node_arr = unique(degree_node_arr)
      
      ## successors部分
      library(igraph)
      all_node_num = length(cyc_node_vector) + length(degree_node_arr)
      all_node_vector = append(cyc_node_vector, degree_node_arr)
      
      g <- make_empty_graph(n = all_node_num)
      
      #分2段给graph中的vertex添加名字
      #1.给cyc上的vertex添加名字 2.给cyc的degree的vertex添加名字
      g <- set.vertex.attribute(g, "name", value=all_node_vector)
      
      
      #######################################################
      #标记 never_considered_comp 在图中这些化合物节点为绿色
      for (i in 1:length(all_node_vector)) {
        if(all_node_vector[i] %in% never_considered_comp_names) {
          V(g)[name==all_node_vector[i]]$color <- "grey"
        }
      }
      
      #######################################################
      ## 新的建环方式
      for (i in 1:length(tmp_cycle_chain_str)) {
        chain_cycle_temp = tmp_cycle_chain_str[i]
        chain_cycle_temp_arr = unlist(strsplit(chain_cycle_temp,split = " -> "))
        
        for (j in 1:length(chain_cycle_temp_arr)) {
          if(j%%2 == 0) {
            edge_info = chain_cycle_temp_arr[j]
            tmp_cyc_hnode = V(g)[name==chain_cycle_temp_arr[j-1]]
            tmp_cyc_tnode = V(g)[name==chain_cycle_temp_arr[j+1]]
            if(edge_info == "G") {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "red")
            }else if(edge_info == "U") {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "green")
            }else {
              g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "black")
            }
          }
        }
      }
      #######################################################
      
      ####################################################
      ## successors部分
      # 根据name查找某个点 V(g)[name=="C00370"]
      #g <- add_edges(g, c(V(g)[name=="C00370"],V(g)[name=="C01958"],V(g)[name=="C00028"]), color = "black")
      
      
      #有入度的结点
      # up_out_node
      # up_out_od
      #有出度的结点
      # up_in_node
      # up_in_ind
      
      if (!identical(up_out_node == 0, logical(0))) {
        for (i in 1:length(up_out_node)) {
          tmp_cyc_node = V(g)[name==up_out_node[i]]
          tmp_cyc_od = V(g)[name==up_out_od[i]]
          if (!are.connected(g, tmp_cyc_node, tmp_cyc_od)) {
            g <- add_edges(g, c(tmp_cyc_node, tmp_cyc_od), color = "cyan")      
          }
      
        }
      }
      
      
      if (!identical(up_in_node == 0, logical(0))) {
        for (i in 1:length(up_in_node)) {
          tmp_cyc_node = V(g)[name==up_in_node[i]]
          tmp_cyc_ind = V(g)[name==up_in_ind[i]]
          if (!are.connected(g, tmp_cyc_ind, tmp_cyc_node)) {
            g <- add_edges(g, c(tmp_cyc_ind, tmp_cyc_node), color = "cyan")  
          }
          
        }
      }
      
      #有入度的结点
      # gap_out_node
      # gap_out_od
      #有出度的结点
      # gap_in_node
      # gap_in_ind
      if (!identical(gap_out_node == 0, logical(0))) {
        for (i in 1:length(gap_out_node)) {
          tmp_cyc_node = V(g)[name==gap_out_node[i]]
          tmp_cyc_od = V(g)[name==gap_out_od[i]]
          ## 注意: 某条边若已经是 up 了, 则不可能再是 gap 了
          if (!are.connected(g, tmp_cyc_node, tmp_cyc_od)) {
            g <- add_edges(g, c(tmp_cyc_node, tmp_cyc_od), color = "orange")
          }
          
        }
      }
      
      if (!identical(gap_in_node == 0, logical(0))) {
        for (i in 1:length(gap_in_node)) {
          tmp_cyc_node = V(g)[name==gap_in_node[i]]
          tmp_cyc_ind = V(g)[name==gap_in_ind[i]]
          ## 注意: 某条边若已经是 up 了, 则不可能再是 gap 了
          if (!are.connected(g, tmp_cyc_ind, tmp_cyc_node)) {
            g <- add_edges(g, c(tmp_cyc_ind, tmp_cyc_node), color = "orange")  
          }
          
        }
      }
      
      ## successors部分
      ####################################################
      
      tmp_graph_name = paste0(tumor_name, "_C", cyc_id_arr[k], ".png")
      tmp_plot_name = file.path(plot_name, tmp_graph_name)
      
      png(tmp_plot_name, 500, 500)
      plot(g)
      dev.off()
    }    
  }

  
}





##############################################################################
my_ug_degnode_3_main <- function(res_path, graph_path, input_tumor_name) {
  
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  
  ##2.each_edge
  ##select_in_upinfo_eachedge select_out_upinfo_eachedge
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_result/select_in_upinfo_eachedge.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_result/select_out_upinfo_eachedge.RData"))
  
  tumors_array = c(input_tumor_name)
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    ##2.each_edge
    res_sub_path = "degree_each_edge"
    dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
    res_file_path = file.path(graph_path, res_sub_path)
    
    plot_up_cycle(tumor_name, res_file_path, select_in_upinfo_eachedge[[tumor_name]], select_out_upinfo_eachedge[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)
    
  }
  
}






# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
# input_tumor_name = "COAD"
# my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)


















##################################################################################
## 暂时用不到的
## 1.up
# plot_name_1 = paste0("graph_cycle/", tumor_name, "/1_up/", tumor_name, "_C")
# plot_up_cycle(tumor_name, plot_name_1, select_in_upinfo[[tumor_name]], select_out_upinfo[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

##3.both_inout
# plot_name_2 = paste0("graph_cycle/", tumor_name, "/3_both_inout/", tumor_name, "_C")
# plot_up_cycle(tumor_name, plot_name_2, select_in_upinfo_both[[tumor_name]], select_out_upinfo_both[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

##4.with_gap
# plot_name_4 = paste0("graph_cycle/", tumor_name, "/4_with_gap/", tumor_name, "_C")
# plot_ug_cycle(tumor_name, plot_name_4, select_in_upinfo[[tumor_name]], select_out_upinfo[[tumor_name]], select_in_gapinfo[[tumor_name]], select_out_gapinfo[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

##5.both_with_gap
# plot_name_5 = paste0("graph_cycle/", tumor_name, "/5_both_with_gap/", tumor_name, "_C")
# plot_ug_cycle(tumor_name, plot_name_5, select_in_upinfo_eachedge[[tumor_name]], select_out_upinfo_eachedge[[tumor_name]], select_in_gapinfo_eachedge[[tumor_name]], select_out_gapinfo_eachedge[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

##6.each_edge_with_gap
# plot_name_6 = paste0("graph_cycle/", tumor_name, "/6_each_edge_with_gap/", tumor_name, "_C")
# plot_ug_cycle(tumor_name, plot_name_6, select_in_upinfo_eachedge[[tumor_name]], select_out_upinfo_eachedge[[tumor_name]], select_in_gapinfo_eachedge[[tumor_name]], select_out_gapinfo_eachedge[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

