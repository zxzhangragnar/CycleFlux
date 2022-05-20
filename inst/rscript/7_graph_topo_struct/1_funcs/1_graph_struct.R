

add_struct_id_cycles <- function(judge, struct_cycle_sort_list, cycle_directed, g) {
  
  #4.all
  for (k in 1:length(cycle_directed[,1])) {
    tmp_cyc_id = cycle_directed[k,"cycid"]
    
    struct_color = "grey"
    struct_size = 3
    
    if (tmp_cyc_id %in% struct_cycle_sort_list$hub) {
      if ((judge == "hub") | (judge == "all")) {
        struct_color = "yellow"
        struct_size = 6 
      }
    }else if (tmp_cyc_id %in% struct_cycle_sort_list$net) {
      if ((judge == "net") | (judge == "all")) {
        struct_color = "green"
        struct_size = 6
      }
    }else if (tmp_cyc_id %in% struct_cycle_sort_list$idv) {
      if ((judge == "idv") | (judge == "all")) {
        struct_color = "red"
        struct_size = 6 
      }
    }
    
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
          V(g)[name==cpds_str[j-1]]$color <- struct_color
          V(g)[name==cpds_str[j+1]]$color <- struct_color
          E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$size = struct_size
          E(g)[get.edge.ids(g, c(tmp_cyc_hnode, tmp_cyc_tnode))]$color = struct_color
        }
      }
    }
  }
  
  return(g)
}





##################################################################################################
plot_ug_cycle<-function(tumor_name, input_pathway_name, plot_name, subnet_edge_flux_list, cycle_directed, struct_cycle_sort_list, never_considered_comp_names) {
  
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
  
  ###
  library(qgraph)
  tmp_graph_name = paste0(tumor_name, "_", input_pathway_name, ".png")
  tmp_plot_name = file.path(plot_name, tmp_graph_name)
  
  png(tmp_plot_name, height=60, width=120, units="in", res=250)
  par(mfrow=c(1, 4))
  e <- get.edgelist(g,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                         area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
  ###
  #1. hub
  g = add_struct_id_cycles("hub", struct_cycle_sort_list, cycle_directed, g) 
  
  # #显示名字
  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("hub forming", side=1)   

  #2. net
  g = add_struct_id_cycles("net", struct_cycle_sort_list, cycle_directed, g) 
  
  # 显示名字
  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("net linking", side=1)   

  #3. idv
  g = add_struct_id_cycles("idv", struct_cycle_sort_list, cycle_directed, g) 
  
  # #显示名字
  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("individual", side=1)   
  
  
  #4.all
  g = add_struct_id_cycles("all", struct_cycle_sort_list, cycle_directed, g) 

  # #显示名字
  plot(g,layout=l,vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
  mtext("all struct", side=1)    

  dev.off()
}














##############################################################################
graph_struct_main <- function(output_path, res_path, graph_path, input_tumor_name, input_pathway_name="hsa_all") {
  # init
  load(file.path(res_path, "6_graph_single_cycle/result_analysis/never_considered_compounds.RData"))
  load(file.path(res_path, "3_flux_subnet/result_final/subnet_edge_flux_list.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/compounds_dict.RData"))
  load(file.path(res_path, "2_cycle_struct_sort/result_struct/struct_cycle_sort_list.RData"))
  load(file.path(output_path, "cycle_directed.RData"))
  
  tumors_array = c(input_tumor_name)

  # save
  res_sub_path = "graph_topo_struct"
  dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(graph_path, res_sub_path)
  
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    plot_ug_cycle(tumor_name, input_pathway_name, res_file_path, subnet_edge_flux_list, cycle_directed, struct_cycle_sort_list, never_considered_comp_names)
  }
  
}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
# input_tumor_name = "COAD"
# input_pathway_name = "hsa00350"
# graph_struct_main(output_path, res_path, graph_path, input_tumor_name, input_pathway_name)
