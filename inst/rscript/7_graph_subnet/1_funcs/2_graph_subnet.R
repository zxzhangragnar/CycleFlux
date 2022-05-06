# 
# 约简后的代谢网络
# 
# 标注出来:
#   (1)在原输入文件dp_part_net中 成环的边(cycle的边)
# (2)在原输入文件dp_part_net中 可去掉的边(gap的边)
# 
# 绘制2个层面的图像:
#   (3)表达express层面图像(edge层面的图像)(切换模式1)
# 处理后的代谢网络 每条边 红:gap 绿:up 黑:cycle 灰:never_consider 白:正常
# 
# "点击某条edge边" 还可以查看来自subnet_edge_flux_list.RData的
# gene_fold_change...等信息 
# 
# 【全部edge的up/gap】
# dp_part_net中【全部的边的gap/up】都表示出来
# 这样就包括了 环的 degnode successor的那些边
# 
# 普通环的边和点也标注出来颜色
#

add_normal_cycles <- function(cycle_directed, g) {
  #添加普通环
  for (k in 1:length(cycle_directed[,1])) {
    ord_cpds_str = cycle_directed[k,"ord_cpds_str"]
    ord_cpds_str = substring(ord_cpds_str, 2,nchar(ord_cpds_str)-1)
    ord_cpds_str = unlist(strsplit(ord_cpds_str,split = ", "))
    for (i in 1:length(ord_cpds_str)) {
      cpds_str = substring(ord_cpds_str[i], 2,nchar(ord_cpds_str[i])-1)
      cpds_str = unlist(strsplit(cpds_str,split = "->")) #'c' 'U' 'c' 'G' 'c'
      for (j in 1:length(cpds_str)) {
        if(j%%2 == 0) {
          tmp_cyc_hnode = V(g)[name==cpds_str[j-1]]
          tmp_cyc_tnode = V(g)[name==cpds_str[j+1]]
          V(g)[name==cpds_str[j-1]]$color <- "grey"
          V(g)[name==cpds_str[j+1]]$color <- "grey"
          E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$size = 1
          E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "grey"
        }
      }
    }
  }
  
  return(g)
}



add_upgap_edges <- function(part_subnet, g, judge) {
  # up gap 的边
  for (i in 1:length(part_subnet[,1])) {
    c_in = as.character(part_subnet[i, "c_in"])
    c_out = as.character(part_subnet[i, "c_out"])
    tmp_cyc_hnode = as.integer(V(g)[name==c_in])
    tmp_cyc_tnode = as.integer(V(g)[name==c_out])
    #size和箭头大小配套使用
    E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$size = 3
    if(part_subnet[i, "ifup"] == "up") {
      if ((judge == "up") | (judge == "all")) {
        E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "cyan" 
      }
    }else if(part_subnet[i, "ifgap"] == "gap") {
      if ((judge == "gap") | (judge == "all")) {
        E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "orange" 
      }
    }
  }
  return(g)
}




add_ug_cycles <- function(tumor_name, gapup_cycle_chain_list, g) {
  
  tmp_chain_str = gapup_cycle_chain_list[[tumor_name]]
  if (length(tmp_chain_str) > 0) {
    cyc_id_arr = names(tmp_chain_str)
    for (k in 1:length(cyc_id_arr)) {
      tmp_cyc_id = cyc_id_arr[k]
      tmp_cycle_chain_str = tmp_chain_str[[tmp_cyc_id]] #[1]C->U->c->G->c [2]c->G->c->U->c
      #######################################################
      ## 新的建环方式
      for (i in 1:length(tmp_cycle_chain_str)) {
        chain_cycle_temp = tmp_cycle_chain_str[i]
        chain_cycle_temp_arr = unlist(strsplit(chain_cycle_temp,split = " -> ")) #'c' 'U' 'c' 'G' 'c'
        
        for (j in 1:length(chain_cycle_temp_arr)) {
          if(j%%2 == 0) {
            edge_info = chain_cycle_temp_arr[j]
            tmp_cyc_hnode = V(g)[name==chain_cycle_temp_arr[j-1]]
            tmp_cyc_tnode = V(g)[name==chain_cycle_temp_arr[j+1]]
            V(g)[name==chain_cycle_temp_arr[j-1]]$color <- "yellow"
            V(g)[name==chain_cycle_temp_arr[j+1]]$color <- "yellow"
            
            
            E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$size = 6
            
            if(edge_info == "G") {
              E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "red"
            }else if(edge_info == "U") {
              E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "green"
            }else {
              E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = "black"
            }
          }
        }
      }
      #######################################################
    }    
  }  
  
  return(g)
}






plot_ug_cycle<-function(tumor_name, input_pathway_name, plot_name, subnet_edge_flux_list, cycle_directed, gapup_cycle_chain_list, never_considered_comp_names) {
  
  part_subnet = subnet_edge_flux_list[[tumor_name]]
  
  #1.画子网
  ## 添加全部点
  all_node_vector = union(part_subnet$c_in, part_subnet$c_out)
  ##
  library(igraph)
  g <- make_empty_graph(n = length(all_node_vector))
  g <- set.vertex.attribute(g, "name", value = all_node_vector)
  #######################################################
  #标记 never_considered_comp 在图中这些化合物节点为cornsilk色
  for (i in 1:length(all_node_vector)) {
    if(all_node_vector[i] %in% never_considered_comp_names) {
      V(g)[name==all_node_vector[i]]$color <- "cornsilk"
    }
  }
  #######################################################
  ## 添加全部边
  for (i in 1:length(part_subnet[,1])) {
    c_in = as.character(part_subnet[i, "c_in"])
    c_out = as.character(part_subnet[i, "c_out"])
    tmp_cyc_hnode = as.integer(V(g)[name==c_in])
    tmp_cyc_tnode = as.integer(V(g)[name==c_out])
    #size和箭头大小配套使用
    
    if (!are.connected(g, tmp_cyc_hnode, tmp_cyc_tnode)) {
      g <- add_edges(g, c(tmp_cyc_hnode,tmp_cyc_tnode), color = "cornsilk", size = 1)
    }

  }
  
  #添加普通环
  g = add_normal_cycles(cycle_directed, g)
  
  ###
  library(qgraph)
  tmp_graph_name = paste0(tumor_name, "_", input_pathway_name, ".png")
  tmp_plot_name = file.path(plot_name, tmp_graph_name)
  
  png(tmp_plot_name, height=60, width=120, units="in", res=250)
  par(mfrow=c(1, 3))
  e <- get.edgelist(g,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                         area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
  ###
  
  #添加ug环
  g = add_ug_cycles(tumor_name, gapup_cycle_chain_list, g)

  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("ug cycle", side=1)   
  
  ## 2.with up
  g = add_upgap_edges(part_subnet, g, "up")
  #添加ug环
  g = add_ug_cycles(tumor_name, gapup_cycle_chain_list, g)

  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("with up", side=1)  
  
  ## 2.with up gap
  g = add_upgap_edges(part_subnet, g, "all")
  #添加ug环
  g = add_ug_cycles(tumor_name, gapup_cycle_chain_list, g)
  
  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("with up gap", side=1)    

  dev.off()
}















##############################################################################
graph_subnet_main <- function(output_path, res_path, graph_path, input_tumor_name, input_pathway_name="hsa_all") {
  # init
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData"))
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  load(file.path(res_path, "3_flux_subnet/result_final/subnet_edge_flux_list.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/compounds_dict.RData"))
  load(file.path(output_path, "res_allpathway_cycle_union_directed.RData"))
  
  ##
  tumors_array = c(input_tumor_name)
  
  # save
  res_sub_path = "graph_subnet"
  dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(graph_path, res_sub_path)
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    plot_ug_cycle(tumor_name, input_pathway_name, res_file_path, subnet_edge_flux_list, cycle_directed, gapup_cycle_chain_list, never_considered_comp_names)
  }
  
}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
# input_tumor_name = "COAD"
# input_pathway_name = "hsa00350"
# graph_subnet_main(output_path, res_path, graph_path, input_tumor_name, input_pathway_name)











