
##############################################################################
##############################################################################
# 画出这些环型结构
#igraph

plot_ug_cycle<-function(tumor_name, plot_name, gapup_cycle_chain_list, never_considered_comp_names) {
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
      library(igraph)
      all_node_num = length(cyc_node_vector)
      all_node_vector = cyc_node_vector
      
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
      
      tmp_graph_name = paste0(tumor_name, "_C", cyc_id_arr[k], ".png")
      tmp_plot_name = file.path(plot_name, tmp_graph_name)
      
      png(tmp_plot_name, 500, 500)
      plot(g)
      dev.off()
    }    
  }

}



plot_ug_cycle_cname<-function(tumor_name, plot_name, compounds_dict, gapup_cycle_chain_list, never_considered_comp_names) {
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
          tmp_cname = compounds_dict[which(compounds_dict$ENTRY == tmp_cycle_chain_arr[i]), "...1"]
          cyc_node_vector = append(cyc_node_vector, tmp_cname)
          # cyc_node_vector = append(cyc_node_vector, tmp_cycle_chain_arr[i])
        }
      }
      
      
      ##
      library(igraph)
      all_node_num = length(cyc_node_vector)
      all_node_vector = cyc_node_vector
      
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
      ## 新的建环方式
      for (i in 1:length(tmp_cycle_chain_str)) {
        chain_cycle_temp = tmp_cycle_chain_str[i]
        chain_cycle_temp_arr = unlist(strsplit(chain_cycle_temp,split = " -> "))
        
        for (j in 1:length(chain_cycle_temp_arr)) {
          hnode_cname = compounds_dict[which(compounds_dict$ENTRY == chain_cycle_temp_arr[j-1]), "...1"]
          tnode_cname = compounds_dict[which(compounds_dict$ENTRY == chain_cycle_temp_arr[j+1]), "...1"]
          
          if(j%%2 == 0) {
            edge_info = chain_cycle_temp_arr[j]
            tmp_cyc_hnode = V(g)[name==hnode_cname]
            tmp_cyc_tnode = V(g)[name==tnode_cname]
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
      
      tmp_graph_name = paste0(tumor_name, "_C", cyc_id_arr[k], ".png")
      tmp_plot_name = file.path(plot_name, tmp_graph_name)
      
      png(tmp_plot_name, 500, 500)
      plot(g)
      dev.off()
    }    
  }

}
#关于water -> C00001
#在all_chain_list 和 all_chain_list_cid的部分


############################################## plot ###################################################################


##############################################################################
cycle_draw_cname_main <- function(res_path, graph_path, input_tumor_name) {
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/compounds_dict.RData"))
  
  tumors_array = c(input_tumor_name)
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    ## cname
    res_sub_path = "cname"
    dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
    res_file_path = file.path(graph_path, res_sub_path)
    plot_ug_cycle_cname(tumor_name, res_file_path, compounds_dict, gapup_cycle_chain_list, never_considered_comp_names)
  }
  
}


cycle_draw_cid_main <- function(res_path, graph_path, input_tumor_name) {
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/compounds_dict.RData"))
  
  tumors_array = c(input_tumor_name)
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    ## cid
    res_sub_path = "cid"
    dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
    res_file_path = file.path(graph_path, res_sub_path)
    plot_ug_cycle(tumor_name, res_file_path, gapup_cycle_chain_list, never_considered_comp_names)
    
  }
  
}






# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
# input_tumor_name = "COAD"
# #cycle_draw_cname_main(res_path, graph_path, input_tumor_name)
# cycle_draw_cid_main(res_path, graph_path, input_tumor_name)

