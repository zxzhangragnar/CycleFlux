

#情况1:  e2<--up--e1<--up--c2<--up--c1
####################################################################

#情况2:  c2<--up--c1<--up--e1<--up--e2
####################################################################

#情况3: 合并情况1和情况2 (_e1)
#e1<--up--(c2<--up--c1)<--up--e1
###################################################################

#情况4: 合并情况1和情况2 (_e1e2)
#e2<--up--e1<--up--(c2<--up--c1)<--up--e1<--up--e2
###################################################################








###########################################################################
###########################################################################
# 画出这些图形
#参考 1_gap_analysis.R

plot_up_cycle<-function(tumor_name, plot_name, select_st2_1_gapinfo_e1e2, select_st2_2_gapinfo_e1e2, select_st3_1_gapinfo_e1e2, select_st3_2_gapinfo_e1e2, gapup_cycle_chain_list, never_considered_comp_names) {
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
      ############################
      #要求同时有入度和出度?
      st2_cycid = select_st2_1_gapinfo_e1e2$cycid
      st3_cycid = select_st3_1_gapinfo_e1e2$cycid
      intersect_cycid = intersect(st2_cycid, st3_cycid)
      
      ist_2_1_e1e2 = select_st2_1_gapinfo_e1e2[which(select_st2_1_gapinfo_e1e2$cycid %in% intersect_cycid),]
      ist_2_2_e1e2 = select_st2_2_gapinfo_e1e2[which(select_st2_2_gapinfo_e1e2$cycid %in% intersect_cycid),]
      
      ist_3_1_e1e2 = select_st3_1_gapinfo_e1e2[which(select_st3_1_gapinfo_e1e2$cycid %in% intersect_cycid),]
      ist_3_2_e1e2 = select_st3_2_gapinfo_e1e2[which(select_st3_2_gapinfo_e1e2$cycid %in% intersect_cycid),]
      ###########################
      
      ist_2_1_node = ist_2_1_e1e2[which(ist_2_1_e1e2$cycid == tmp_cyc_id), "node"]
      ist_2_1_od = ist_2_1_e1e2[which(ist_2_1_e1e2$cycid == tmp_cyc_id), "od"]
      ist_2_2_od = ist_2_2_e1e2[which(ist_2_2_e1e2$cycid == tmp_cyc_id), "od"]
      ist_2_2_suc = ist_2_2_e1e2[which(ist_2_2_e1e2$cycid == tmp_cyc_id), "suc"]
      
      ist_3_1_node = ist_3_1_e1e2[which(ist_3_1_e1e2$cycid == tmp_cyc_id), "node"]
      ist_3_1_ind = ist_3_1_e1e2[which(ist_3_1_e1e2$cycid == tmp_cyc_id), "ind"]
      ist_3_2_ind = ist_3_2_e1e2[which(ist_3_2_e1e2$cycid == tmp_cyc_id), "ind"]
      ist_3_2_presuc = ist_3_2_e1e2[which(ist_3_2_e1e2$cycid == tmp_cyc_id), "presuc"]
      
      #degree_node_arr1 = c(rbind(ist_2_1_od,ist_2_2_od,ist_2_2_suc,ist_3_1_ind,ist_3_2_ind,ist_3_2_presuc))
      degree_node_arr = c(ist_2_1_od,ist_2_2_od,ist_2_2_suc,ist_3_1_ind,ist_3_2_ind,ist_3_2_presuc)
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
      
      if (!identical(ist_2_1_node == 0, logical(0))) {
        for (i in 1:length(ist_2_1_node)) {
          ist_2_node = as.integer(V(g)[name==ist_2_1_node[i]])
          ist_2_od = as.integer(V(g)[name==ist_2_1_od[i]]) #ist_2_1_od == ist_2_2_od
          ist_2_suc = as.integer(V(g)[name==ist_2_2_suc[i]])
          
          if (!are.connected(g, ist_2_node, ist_2_od)) {
            g <- add_edges(g, c(ist_2_node,ist_2_od), color = "cyan")
          }
          if (!are.connected(g, ist_2_od, ist_2_suc)) {
            g <- add_edges(g, c(ist_2_od,ist_2_suc), color = "cyan")
          }
          
        }      
      }
      
      
      if (!identical(ist_3_1_node == 0, logical(0))) {
        for (i in 1:length(ist_3_1_node)) {
          ist_3_node = as.integer(V(g)[name==ist_3_1_node[i]])
          ist_3_ind = as.integer(V(g)[name==ist_3_1_ind[i]]) #ist_3_1_od == ist_3_2_od
          ist_3_presuc = as.integer(V(g)[name==ist_3_2_presuc[i]])
          
          if (!are.connected(g, ist_3_ind, ist_3_node)) {
            g <- add_edges(g, c(ist_3_ind,ist_3_node), color = "cyan")
          }
          
          if (!are.connected(g, ist_3_presuc, ist_3_ind)) {
            g <- add_edges(g, c(ist_3_presuc,ist_3_ind), color = "cyan")
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
suc2_my_ug_degnode_3_main <- function(res_path, graph_path, input_tumor_name) {
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_result/select_st2_1_gapinfo_e1e2.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_result/select_st2_2_gapinfo_e1e2.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_result/select_st3_1_gapinfo_e1e2.RData"))
  load(file.path(res_path, "6_graph_single_cycle/2_funcs_suc2/cyc_successors_result/select_st3_2_gapinfo_e1e2.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  
  tumor_name = input_tumor_name
  res_sub_path = "degree2_each_edge"
  dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(graph_path, res_sub_path)

  plot_up_cycle(tumor_name, res_file_path, select_st2_1_gapinfo_e1e2, select_st2_2_gapinfo_e1e2, select_st3_1_gapinfo_e1e2, select_st3_2_gapinfo_e1e2, gapup_cycle_chain_list, never_considered_comp_names)

}



# test
res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
input_tumor_name = "COAD"
suc2_my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)






