#在执行这个函数之前 要先执行 0_rela_module_cyc.R


##################################### target ################################################
#############################################################################################
#
#
# 
# index
# 将节点的H-指数定义为最大值H，使得至少存在H个度不小于H的邻居。
# 求h-index算法
# 对某个节点，将其所有邻居的度从小到大排列，从度最大(设度为d)的邻居开始遍历，每次遍历其它所有邻居，
# 看是否存在d以上个度大于等于d的邻居节点，若是则h-index=d，若不存在，则d=d-1继续遍历直到d=1
# 
#  (degree, H-index and coreness) 和 (closeness and betweenness) 
#
#


########################### 在以下计算中，将"环"看成网络中的节点 ####################################



##################################### part1 H-index ################################################
#############################################################################################

# index: 将节点的H-指数定义为最大值H，使得至少存在H个度不小于H的邻居。
# 求h-index算法
# 对某个节点，将其所有邻居的度从小到大排列，
# 从度最大(设度为d)的邻居开始遍历，每次遍历其它所有邻居，
# 看是否存在d以上个度大于等于d的邻居节点，若是则h-index=d，若不存在，则d=d-1继续遍历直到d=1


#计算h-index的函数
get_h_index_tool_func <- function(cites) {
  if(max(cites) == 0) return(0) # assuming this is reasonable
  cites = cites[order(cites, decreasing = TRUE)]
  tail(which(cites >= seq_along(cites)), 1)
}

# 计算每个环的h-index
get_cyc_hindex <- function(cycle_directed) {
  cyc_hindex_df = as.data.frame(cycle_directed[,1])
  names(cyc_hindex_df) = "cycid"
  for (r in 1:length(rownames(cycle_directed))) {
    print(r)
    cycid = cycle_directed[r, "cycid"]
    temp_cyc_nghnode = cycle_directed[which(cycle_directed$cycid == cycid), "neighcycs"]
    
    if(temp_cyc_nghnode == "[]") {
      temp_hindex = 0
    }else {
      temp_cyc_nghnode = substring(temp_cyc_nghnode, 2,nchar(temp_cyc_nghnode)-1)
      temp_cyc_nghnode <- unlist(strsplit(temp_cyc_nghnode,split = ", "))
      
      #neighbor degree list
      temp_cyc_ngbnode_degarr = list()
      for (i in 1:length(temp_cyc_nghnode)) {
        node_name = temp_cyc_nghnode[i]
        temp_cyc_ngbnode_degarr[[node_name]] = cycle_directed[which(cycle_directed$cycid == node_name), "neighcyc"]
      }
      
      # sort the list
      temp_cyc_ngbnode_degarr = unlist(temp_cyc_ngbnode_degarr)
      temp_cyc_ngbnode_degarr = sort(temp_cyc_ngbnode_degarr, decreasing = TRUE)
      temp_cyc_ngbnode_degarr = as.list(temp_cyc_ngbnode_degarr)
      
      # put list into vec
      temp_cyc_ngbnode_degvec = c()
      for (i in 1:length(temp_cyc_ngbnode_degarr)) {
        temp_cyc_ngbnode_degvec = append(temp_cyc_ngbnode_degvec, temp_cyc_ngbnode_degarr[[i]])
      }
      temp_hindex = get_h_index_tool_func(temp_cyc_ngbnode_degvec)
    }
    
    
    cyc_hindex_df[which(cyc_hindex_df$cycid == cycid), "hindex"] = temp_hindex
  }
  
  return(cyc_hindex_df)
}




##################################### part2 coreness ################################################
#############################################################################################

# 求coreness算法
# https://www.rdocumentation.org/packages/igraph/versions/1.2.6/topics/coreness
# The k-core of graph is a maximal subgraph in which each vertex has at least degree k. 
# The coreness of a vertex is k if it belongs to the k-core but not to the (k+1)-core.

# 图的 k-core 是一个每个顶点至少有 k 度的极大子图。
# 如果顶点属于 k-core 但不属于 (k+1)-core，则顶点的核度为 coreness。

get_graph_cyc_df <- function(cycle_directed) {
  graph_cyc_df = data.frame()
  for (i in 1:length(rownames(cycle_directed))) {
    print(i)
    cycid = cycle_directed[i, "cycid"]
    temp_cyc_nghnode = cycle_directed[i, "neighcycs"]
    if(temp_cyc_nghnode == "[]") {
      # pass
      temp_df = data.frame(
        from = c(cycid), 
        to = c(cycid)
      )
      graph_cyc_df = rbind(graph_cyc_df, temp_df)
    }else {
      temp_cyc_nghnode = substring(temp_cyc_nghnode, 2,nchar(temp_cyc_nghnode)-1)
      temp_cyc_nghnode <- unlist(strsplit(temp_cyc_nghnode,split = ", "))
      
      for (j in 1:length(temp_cyc_nghnode)) {
        temp_df = data.frame(
          from = c(cycid), 
          to = c(temp_cyc_nghnode[j])
        )
        graph_cyc_df = rbind(graph_cyc_df, temp_df)
      }  
    }
    
  }
  return(graph_cyc_df)
}

get_cyc_coreness <- function(graph_cyc_df) {
  #build graph
  library(igraph)
  cycle_graph <- graph.data.frame(d = graph_cyc_df, directed = FALSE)
  cyc_coreness = coreness(cycle_graph) 		# small core triangle in a ring
  cyc_coreness = as.data.frame(cyc_coreness)
  #cyc_coreness = as.data.frame(cyc_coreness[-length(rownames(cyc_coreness)),])
  names(cyc_coreness) = "coreness"  
  return(cyc_coreness)
}



#############################################################################################
# test
cycle_stcid_topology_main <- function(output_path, res_path) {
  load(file.path(output_path, "res_allpathway_cycle_union_directed.RData"))
  load(file.path(res_path, "1_cycle_topology/result_topo/cycle_rela_module.RData"))
  
  cyc_hindex_df = get_cyc_hindex(cycle_directed)
  graph_cyc_df = get_graph_cyc_df(cycle_directed)
  cyc_coreness_df = get_cyc_coreness(graph_cyc_df)
  
  cyc_hindex_df = as.data.frame(cyc_hindex_df[,"hindex"])
  names(cyc_hindex_df) = c("h-index")
  cycle_topology_info = cbind(cycle_rela_module, cyc_hindex_df, cyc_coreness_df)
  
  cycle_directed_topoinfo = cycle_directed[,c("neighcyc", "deg", "indeg", "outdeg", "sharedcyc")]
  cycle_topology_info = cbind(cycle_topology_info, cycle_directed_topoinfo)
  
  # merge distance
  distance_col = as.data.frame(cycle_directed[,"distance"])
  names(distance_col) = c("distance")
  cycle_topology_info = cbind(cycle_topology_info, distance_col)
  
  # save
  res_sub_path = "1_cycle_topology/result_topo"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_topology_info.RData")
  
  #save(cycle_topology_info, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cycle_topology_info.RData")
  save(cycle_topology_info, file=res_file_path)
  
}



# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# cycle_stcid_topology_main(output_path, res_path) 
# 































