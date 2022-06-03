


#######################################################################################
# source("1_cycle_directed.R")


separate_irreversible_net <- function(net) {
  for (i in 1:length(rownames(net))) {
    x = net[i,]
    if(is.na(x$Irreversible) | (x$Irreversible==TRUE)) {
      new_row = x
      new_row$C_out = x$C_in
      new_row$C_in = x$C_out
      net = rbind(net, new_row)
    }
  }
  net = unique(net)
  return(net)
}

build_net <- function(hsa_net) {
  print("graph build")
  all_compounds = union(hsa_net$C_in, hsa_net$C_out)
  hsa_net = separate_irreversible_net(hsa_net)

  ## vertex
  library(igraph)
  g <- make_empty_graph(n = length(all_compounds))
  g <- set.vertex.attribute(g, "name", value = all_compounds)

  ## edges
  for (i in 1:length(rownames(hsa_net))) {
    c_in = V(g)[name==as.character(hsa_net[i, "C_in"])]
    c_out = V(g)[name==as.character(hsa_net[i, "C_out"])]

    rid = as.character(hsa_net[i, "Rid"])
    enzyme = as.character(hsa_net[i, "EC"])
    pathway = as.character(hsa_net[i, "Pathway"])
    gene_symbol = as.character(hsa_net[i, "Gene_symbol"])

    if (!are.connected(g, c_in, c_out)) {
      g <- add_edges(g, c(c_in, c_out),rid = rid,enzyme = enzyme,pathway = pathway,gene_symbol = gene_symbol)
    }else {
      old_rid = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$rid, split = ";"))
      old_enzyme = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$enzyme, split = ";"))
      old_pathway = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$pathway, split = ";"))
      old_gene_symbol = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$gene_symbol, split = ";"))

      new_rid = union(old_rid, rid)
      new_enzyme = union(old_enzyme, enzyme)
      new_pathway = union(old_pathway, pathway)
      new_gene_symbol = union(old_gene_symbol, gene_symbol)

      E(g)[get.edge.ids(g, c(c_in, c_out))]$rid = paste(new_rid, collapse = ";")
      E(g)[get.edge.ids(g, c(c_in, c_out))]$enzyme = paste(new_enzyme, collapse = ";")
      E(g)[get.edge.ids(g, c(c_in, c_out))]$pathway = paste(new_pathway, collapse = ";")
      E(g)[get.edge.ids(g, c(c_in, c_out))]$gene_symbol = paste(new_gene_symbol, collapse = ";")
    }
  }
  return(g)
}

get_template_cm <- function(hsa_net, connected_compounds) {

  template_cm = matrix(0,
                       nrow = length(connected_compounds),
                       ncol = length(connected_compounds),
                       dimnames = list(connected_compounds, connected_compounds))

  connected_net = hsa_net[which((hsa_net$C_in %in% connected_compounds) & (hsa_net$C_out %in% connected_compounds)),]
  for (i in 1:length(rownames(connected_net))) {
    temp_cin = as.character(connected_net[i, "C_in"])
    temp_cout = as.character(connected_net[i, "C_out"])

    template_cm[temp_cin, temp_cout] = 1
    template_cm[temp_cout, temp_cin] = 1
  }

  return(template_cm)

}


find_directed_cycle <- function(hsa_net, g) {
  print("find basic cycles")
  library(igraph)
  compound_graph <- as.undirected(g)
  compound_components = igraph::groups(components(compound_graph))
  connected_compounds = as.vector(compound_components[[1]])

  template_cm = get_template_cm(hsa_net, connected_compounds)

  library('ggm')
  fund_cycles = fundCycles(template_cm)
  cycles = fund_cycles[lapply(fund_cycles, length) > 4]

  replace_cname <- function(x) {
    x <- connected_compounds[x]
    y = c(x[1:(length(x)/2)], x[1])
    return(y)
  }

  cycles_cname <- lapply(cycles, replace_cname)

  # judge directed cycle
  judge_directed_cycle <- function(x) {
    pos_judge = TRUE
    neg_judge = TRUE

    for (i in 1:(length(x)-1)) {
      edge_c = c(x[i], x[i+1])
      ud_edge = hsa_net[which( ((hsa_net$C_in %in% edge_c) & (hsa_net$C_out %in% edge_c)) ),]

      d_edge = hsa_net[which( ((hsa_net$C_in == x[i]) & (hsa_net$C_out == x[i+1])) ),]
      if (!(TRUE %in% ud_edge$Irreversible) & length(rownames(d_edge)) == 0) {
        pos_judge = FALSE
      }
    }

    neg_x = rev(x)
    for (i in 1:(length(neg_x)-1)) {
      edge_c = c(neg_x[i], neg_x[i+1])
      ud_edge = hsa_net[which( ((hsa_net$C_in %in% edge_c) & (hsa_net$C_out %in% edge_c)) ),]

      d_edge = hsa_net[which( ((hsa_net$C_in == neg_x[i]) & (hsa_net$C_out == neg_x[i+1])) ),]
      if (!(TRUE %in% ud_edge$Irreversible) & length(rownames(d_edge)) == 0) {
        neg_judge = FALSE
      }

    }

    res = c()
    if (pos_judge) {
      res =  append(res, paste(x, collapse  = "->"))
    }
    if (neg_judge) {
      res = append(res, paste(neg_x, collapse  = "->"))
    }

    return(res)

  }

  cycles_cname_directed = lapply(cycles_cname, judge_directed_cycle)
  directed_cycle = cycles_cname_directed[lapply(cycles_cname_directed, length) > 0]

  return(directed_cycle)
}


#cycle_directed OK
get_cycle_directed <- function(hsa_net, g) {
  print("cycle_directed")
  directed_cycle = find_directed_cycle(hsa_net, g)

  directed_cycle_node = list()
  cycle_directed = data.frame()

  if(length(directed_cycle) > 0) {
    for (i in 1:length(directed_cycle)) {
      compound_chain = directed_cycle[[i]]
      cycle_node = compound_chain[1]
      cycle_node = unlist(strsplit(cycle_node, split = "->"))
      directed_cycle_node[[i]] = cycle_node
    }

    for(i in 1:length(directed_cycle)) {
      print(paste0("cyc_", i))
      cycle_directed[i, "cycle_id"] = i
      cycle_directed[i, "compound_chain"] = paste(directed_cycle[[i]], collapse = ";")

    }
  }

  return(cycle_directed)
}


##
#cycle_edge_expression OK
get_cycle_edge_expression <- function(cycle_directed, g) {
  print("cycle_edge")
  cycle_edge_expression = data.frame()
  for(i in 1:length(rownames(cycle_directed))) {
    cycle_id = as.integer(cycle_directed[i, "cycle_id"])
    compound_chain = as.character(cycle_directed[i, "compound_chain"])
    compound_chain = unlist(strsplit(compound_chain, split = ";"))
    is_used = FALSE
    temp_cycle_df = data.frame()
    for (j in 1:length(compound_chain)) {
      cycle_node = unlist(strsplit(compound_chain[j], split = "->"))
      temp_df = data.frame()
      for (k in 1:(length(cycle_node)-1)) {
        c_in = cycle_node[k]
        c_out = cycle_node[k+1]

        temp_df[k, "cycle_id"] = cycle_id
        temp_df[k, "c_in"] = c_in
        temp_df[k, "c_out"] = c_out
        temp_df[k, "rid"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$rid
        temp_df[k, "enzyme"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$enzyme
        temp_df[k, "pathway"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$pathway
        temp_df[k, "gene_symbol"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$gene_symbol
      }

      #filter cycle belongs to same ec
      ec_kinds = aggregate(temp_df$enzyme, list(temp_df$cycle_id), paste, collapse = ";")
      ec_kinds = unlist(strsplit(ec_kinds$x, split = ";"))
      ec_kinds = unique(ec_kinds)
      if(length(ec_kinds) > 1) {
        is_used = TRUE
        temp_cycle_df = rbind(temp_cycle_df, temp_df)
      }
    }

    if(is_used) {
      cycle_edge_expression = rbind(cycle_edge_expression, temp_cycle_df)
    }

  }

  return(cycle_edge_expression)
}



rm_disuse_cycle_directed <- function(cycle_edge_expression, cycle_directed) {
  used_cycle_id = unique(cycle_edge_expression$cycle_id)
  cycle_directed = cycle_directed[which(cycle_directed$cycle_id %in% used_cycle_id),]
  cycle_directed[, "cycle_id"] = c(1:length(rownames(cycle_directed)))
  rownames(cycle_directed) = cycle_directed[, "cycle_id"]
  return(cycle_directed)
}

reid_cycle_edge_expression <- function(cycle_edge_expression) {
  used_cycle_id = unique(cycle_edge_expression$cycle_id)
  new_cycle_id = c(1:length(used_cycle_id))
  for (i in 1:length(used_cycle_id)) {
    cycle_edge_expression[which(cycle_edge_expression$cycle_id == used_cycle_id[i]), "cycle_id"] = new_cycle_id[i]
  }
  return(cycle_edge_expression)
}


##
#subnet_edge_expression OK
get_subnet_edge_expression <- function(g) {
  print("subnet_edge")
  subnet_edge_expression = data.frame()
  edges = E(g)
  for(i in 1:length(edges)) {
    subnet_edge_expression[i, "c_in"] = head_of(g, edges[[i]])$name
    subnet_edge_expression[i, "c_out"] = tail_of(g, edges[[i]])$name
    subnet_edge_expression[i, "rid"] = edges[[i]]$rid
    subnet_edge_expression[i, "enzyme"] = edges[[i]]$enzyme
    subnet_edge_expression[i, "pathway"] = edges[[i]]$pathway
    subnet_edge_expression[i, "gene_symbol"] = edges[[i]]$gene_symbol
  }

  return(subnet_edge_expression)
}


##
#cycle_edge_expression OK
get_express_string <- function(cycle_var) {
  cycle_var_table = as.data.frame(table(cycle_var))
  cycle_var_arr = c()
  for (j in 1:length(rownames(cycle_var_table))) {
    stat_name = cycle_var_table[j,1]
    stat_val = cycle_var_table[j,2]
    stat_str = paste0(stat_name, ":", stat_val)
    cycle_var_arr = c(cycle_var_arr, stat_str)
  }
  cycle_var_str = paste(cycle_var_arr, collapse = ";")

  return(cycle_var_str)
}

get_cycle_expression <- function(cycle_directed, g) {
  print("cycle")
  cycle_expression = data.frame()
  for(i in 1:length(rownames(cycle_directed))) {
    cycle_id = as.integer(cycle_directed[i, "cycle_id"])
    compound_chain = as.character(cycle_directed[i, "compound_chain"])
    compound_chain = unlist(strsplit(compound_chain, split = ";"))

    cycle_enzyme = c()
    cycle_pathway = c()
    cycle_gene_symbol = c()
    for (j in 1:length(compound_chain)) {
      cycle_node = unlist(strsplit(compound_chain[j], split = "->"))
      for (k in 1:(length(cycle_node)-1)) {
        c_in = cycle_node[k]
        c_out = cycle_node[k+1]

        edge_enzyme = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$enzyme, split = ";"))
        edge_pathway = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$pathway, split = ";"))
        edge_gene_symbol = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$gene_symbol, split = ";"))

        cycle_enzyme = c(cycle_enzyme, edge_enzyme)
        cycle_pathway = c(cycle_pathway, edge_pathway)
        cycle_gene_symbol = c(cycle_gene_symbol, edge_gene_symbol)
      }
    }

    cycle_enzyme_str = get_express_string(cycle_enzyme)
    cycle_pathway_str = get_express_string(cycle_pathway)
    cycle_gene_symbol_str = get_express_string(cycle_gene_symbol)

    cycle_expression[i, "cycle_id"] = cycle_id
    cycle_expression[i, "ec_express"] = cycle_enzyme_str
    cycle_expression[i, "pathway_express"] = cycle_pathway_str
    cycle_expression[i, "gene_express"] = cycle_gene_symbol_str
  }
  return(cycle_expression)

}





##
#cycle_edgesucs_expression_in OK
get_cycle_edgesucs_expression_in <- function(cycle_directed, g) {
  print("cycle_indegree")
  cycle_edgesucs_expression_in = data.frame()
  for(i in 1:length(rownames(cycle_directed))) {
    if((i %% 10) == 0) {
      print(paste0("cyc_", i))
    }
    cycle_id = as.integer(cycle_directed[i, "cycle_id"])
    compound_chain = as.character(cycle_directed[i, "compound_chain"])
    compound_chain = unlist(strsplit(compound_chain, split = ";"))
    cycle_node = unlist(strsplit(compound_chain[1], split = "->"))

    for (k in 1:(length(cycle_node)-1)) {
      node = as.character(cycle_node[k])
      ind_nodes = neighbors(g, node, mode = "in")$name
      ind_nodes = setdiff(ind_nodes, cycle_node)
      for (m in ind_nodes) {
        temp_row = data.frame()

        ind = as.character(m)
        temp_row[1, "cycle_id"] = cycle_id
        temp_row[1, "node"] = node
        temp_row[1, "ind"] = ind

        temp_row[1, "rid"] = E(g)[get.edge.ids(g, c(ind, node))]$rid
        temp_row[1, "enzyme"] = E(g)[get.edge.ids(g, c(ind, node))]$enzyme
        temp_row[1, "pathway"] = E(g)[get.edge.ids(g, c(ind, node))]$pathway
        temp_row[1, "gene_symbol"] = E(g)[get.edge.ids(g, c(ind, node))]$gene_symbol

        cycle_edgesucs_expression_in = rbind(cycle_edgesucs_expression_in, temp_row)
      }
    }
  }

  return(cycle_edgesucs_expression_in)
}



##
#cycle_edgesucs_expression_out OK
get_cycle_edgesucs_expression_out <- function(cycle_directed, g) {
  print("cycle_outdegree")
  cycle_edgesucs_expression_out = data.frame()
  for(i in 1:length(rownames(cycle_directed))) {
    if((i %% 10) == 0) {
      print(paste0("cyc_", i))
    }
    cycle_id = as.integer(cycle_directed[i, "cycle_id"])
    compound_chain = as.character(cycle_directed[i, "compound_chain"])
    compound_chain = unlist(strsplit(compound_chain, split = ";"))
    cycle_node = unlist(strsplit(compound_chain[1], split = "->"))

    for (k in 1:(length(cycle_node)-1)) {
      node = as.character(cycle_node[k])
      od_nodes = neighbors(g, node, mode = "out")$name
      od_nodes = setdiff(od_nodes, cycle_node)
      for (m in od_nodes) {
        temp_row = data.frame()

        od = as.character(m)
        temp_row[1, "cycle_id"] = cycle_id
        temp_row[1, "node"] = node
        temp_row[1, "od"] = od

        temp_row[1, "rid"] = E(g)[get.edge.ids(g, c(node, od))]$rid
        temp_row[1, "enzyme"] = E(g)[get.edge.ids(g, c(node, od))]$enzyme
        temp_row[1, "pathway"] = E(g)[get.edge.ids(g, c(node, od))]$pathway
        temp_row[1, "gene_symbol"] = E(g)[get.edge.ids(g, c(node, od))]$gene_symbol

        cycle_edgesucs_expression_out = rbind(cycle_edgesucs_expression_out, temp_row)
      }
    }
  }

  return(cycle_edgesucs_expression_out)
}

#######################################################################################
# source("2_get_missing_gene.R")

# all genes in cycles
get_gene_expressinfo<-function(subnet_edge_expression){
  cyc_gene_express_list = c()

  for (i in 1:length(subnet_edge_expression[,1])) {
    gene_num_arr = as.character(subnet_edge_expression[i,"gene_symbol"])
    gene_num_arr = unlist(strsplit(gene_num_arr,split = ";"))
    for (j in 1:length(gene_num_arr)) {
      cyc_gene_express_list = append(cyc_gene_express_list,gene_num_arr[j])
    }
  }
  cyc_gene_express_list = unique(cyc_gene_express_list)
  return(cyc_gene_express_list)
}





# get_cyc_gene_not_in_tgca<-function(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1){
#   cyc_gene_kegg = c()
#   cyc_gene_tcga = c()
#   for (i in 1:length(rownames(to_find_correct_symbol_v1_1))) {
#     kegg_gene = to_find_correct_symbol_v1_1[i,1]
#     tcga_gene = to_find_correct_symbol_v1_1[i,2]
#     if(kegg_gene %in% cyc_gene_expressinfo){
#       cyc_gene_kegg = append(cyc_gene_kegg,kegg_gene)
#     }
#     if(tcga_gene %in% cyc_gene_expressinfo){
#       cyc_gene_tcga = append(cyc_gene_tcga,tcga_gene)
#     }
#   }
#
#   cyc_gene_not_in_tgca = c()
#   cyc_gene_in_tgca = c()
#   tcga_gene_list = rownames(all_gene_stat[[1]])
#   cyc_gene_expressinfo = unique(cyc_gene_expressinfo)
#   for (i in 1:length(cyc_gene_expressinfo)) {
#     gene = cyc_gene_expressinfo[i]
#     if(gene %in% tcga_gene_list){
#       cyc_gene_in_tgca = append(cyc_gene_in_tgca,gene)
#     }else{
#       cyc_gene_not_in_tgca = append(cyc_gene_not_in_tgca,gene)
#     }
#   }
#
#   return(cyc_gene_not_in_tgca)
# }
#
#
# get_gene_missing_list<-function(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1){
#   #test
#   cyc_gene_not_in_tgca = get_cyc_gene_not_in_tgca(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1)
#   cyc_gene_not_in_tgca_but_in_tofind = intersect(cyc_gene_not_in_tgca,to_find_correct_symbol_v1_1$X1)
#   cyc_gene_not_in_tgca_not_in_tofind = setdiff(cyc_gene_not_in_tgca, cyc_gene_not_in_tgca_but_in_tofind)
#
#   #result
#   gene_missing_list = list()
#   gene_missing_list[["cyc_gene_expressinfo"]] = cyc_gene_expressinfo
#   gene_missing_list[["cyc_gene_not_in_tgca"]] = cyc_gene_not_in_tgca
#   gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]] = cyc_gene_not_in_tgca_but_in_tofind
#   gene_missing_list[["cyc_gene_not_in_tgca_not_in_tofind"]] = cyc_gene_not_in_tgca_not_in_tofind
#
#   return(gene_missing_list)
# }
#
#
#
#
#
# get_gene_and_tofind_list_main <- function(subnet_edge_expression, stat_all, deg_all) {
#
#
#   cyc_gene_expressinfo = get_gene_expressinfo(subnet_edge_expression)
#
#   library(readr)
#   to_find_correct_symbol_v1_1 = read_csv(system.file("data/to_find_correct_symbol_v1.1.csv", package = "CycleFlux"), col_names = FALSE)
#
#   gene_missing_list = get_gene_missing_list(cyc_gene_expressinfo, stat_all, to_find_correct_symbol_v1_1)
#
#   return(gene_missing_list)
# }

#######################################################################################
# function for gap
# get_mean_fc<-function(tumor_name, edge_state_list, all_gene_stat) {
#
#   temp_cycle_edge_flux = edge_state_list[[tumor_name]]
#   temp_gene_stat = all_gene_stat[[tumor_name]]
#   temp_cycle_edge_flux[,"mean_fc"] = 0
#
#   for (i in 1:length(temp_cycle_edge_flux[,1])) {
#     temp_edge_sum_fc = 0
#     temp_gene_symbol = temp_cycle_edge_flux[i,"gene_symbol"]
#     temp_gene_symbol_arr = unlist(strsplit(temp_gene_symbol,split = ";"))
#
#     for (j in 1:length(temp_gene_symbol_arr)) {
#       temp_gene = temp_gene_symbol_arr[j]
#       if (temp_gene %in% rownames(temp_gene_stat)) {
#         genesign = as.double(temp_gene_stat[temp_gene,"FC"])
#         temp_edge_sum_fc = temp_edge_sum_fc + genesign
#       }
#     }
#     temp_cycle_edge_flux[i,"mean_fc"] = temp_edge_sum_fc/length(temp_gene_symbol_arr) #more equal
#   }
#
#   return(temp_cycle_edge_flux)
# }

get_state_edges <- function(tumor_name, edge_state_list, gene_gapup_info, all_gene_stat, model) {

  temp_edge_state = edge_state_list[[tumor_name]]
  if (model == 1) {
    temp_edge_flux = get_state_edges_model_1(tumor_name, temp_edge_state, gene_gapup_info)
  }else if (model == 2) {
    temp_edge_flux = get_state_edges_model_2(tumor_name, temp_edge_state, gene_gapup_info, all_gene_stat)
  }else {
    temp_edge_flux = get_state_edges_model_1(tumor_name, temp_edge_state, gene_gapup_info)
  }
  return(temp_edge_flux)
}

get_state_edges_model_1 <- function(tumor_name, temp_edge_state, gene_gapup_info) {
  temp_gene_gapup_info = gene_gapup_info[[tumor_name]]
  temp_edge_flux = temp_edge_state
  temp_edge_flux[, "state"] = "no_change"
  for (i in 1:length(temp_edge_flux[,1])) {
    if((i %% 1000) == 0) {
      print(paste0("edge_", i))
    }
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr = unlist(strsplit(gene_arr, split = ";"))

    temp_gene_state = "no_change"
    temp_gap_num = 0
    temp_up_num = 0
    temp_ne_num = 0
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if (temp_gene %in% names(temp_gene_gapup_info)) {
        temp_gene_state = temp_gene_gapup_info[[temp_gene]]
      }
      if(temp_gene_state == "gap"){temp_gap_num = temp_gap_num + 1}
      if(temp_gene_state == "up"){temp_up_num = temp_up_num + 1}
      if(temp_gene_state == "no_express"){temp_ne_num = temp_ne_num + 1}
    }

    if(temp_ne_num > 0) {temp_edge_flux[i, "state"] = "no_express"}
    if(temp_up_num > 0) {temp_edge_flux[i, "state"] = "up"}
    if(temp_gap_num > 0) {temp_edge_flux[i, "state"] = "gap"}

  }
  return(temp_edge_flux)
}

get_state_edges_model_2 <- function(tumor_name, temp_edge_state, gene_gapup_info, all_gene_stat) {

  temp_gene_gapup_info = gene_gapup_info[[tumor_name]]
  temp_all_gene_stat = all_gene_stat[[tumor_name]]
  temp_edge_flux = temp_edge_state
  temp_edge_flux[, "state"] = "no_change"
  for (i in 1:length(temp_edge_flux[,1])) {
    if((i %% 1000) == 0) {
      print(paste0("edge_", i))
    }
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr <- unlist(strsplit(gene_arr, split = ";"))

    gene_mt_df = data.frame()
    ## find mt(mean tumor value) gene
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if (temp_gene %in% rownames(temp_all_gene_stat)) {
        temp_mt = temp_all_gene_stat[temp_gene, "mt"]
      }else {
        temp_mt = 0
      }
      gene_mt_df[j, "gene"] = temp_gene
      gene_mt_df[j, "mt"] = temp_mt
    }

    gene_mt_df = gene_mt_df[order(gene_mt_df$mt,decreasing = TRUE),]
    max_mt_gene = gene_mt_df[1, "gene"]
    if ((max_mt_gene %in% names(temp_gene_gapup_info))) {
      temp_edge_flux[i, "state"] = temp_gene_gapup_info[[max_mt_gene]]
    }
  }
  return(temp_edge_flux)
}

# get_cycle_edge_obvs <- function(tumor_name, all_gene_stat, cycle_edge_expression, prm_4, prm_5) {
#   temp_all_gene_stat = all_gene_stat[[tumor_name]]
#
#   cycle_edge_obvs = data.frame()
#
#   for (i in 1:length(cycle_edge_expression[,1])) {
#     gene_str = cycle_edge_expression[i,"gene_symbol"]
#     gene_arr = unlist(strsplit(gene_str,split = ";"))
#     gene_num = 0
#     obvs_num = 0
#     for (j in 1:length(gene_arr)) {
#       temp_gene = gene_arr[j]
#       temp_pvalue = 0
#       temp_fc = 0
#       if (temp_gene %in% rownames(temp_all_gene_stat)) {
#         temp_pvalue = temp_all_gene_stat[temp_gene,"p.value"]
#         temp_fc = temp_all_gene_stat[temp_gene,"FC"]
#         gene_num = gene_num + 1
#       }
#
#       if ((temp_pvalue < prm_4) & (abs(temp_fc) > prm_5)) {
#         obvs_num = obvs_num + 1
#       }
#     }
#     obvs_rat = obvs_num/gene_num
#     cycle_edge_obvs[i,"DE_cof"] = round(obvs_rat, 2)
#   }
#   return(cycle_edge_obvs)
# }


## main function
get_edge_state_list <- function(edge_expression, tumors_array, all_gene_stat, gene_gapup_info, model) {
  edge_state_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    print(tumor_name)
    edge_state_list[[tumor_name]] = edge_expression

    edge_state_list[[tumor_name]] = get_state_edges(tumor_name, edge_state_list, gene_gapup_info, all_gene_stat, model)
    #edge_state_list[[tumor_name]] = get_mean_fc(tumor_name, edge_state_list, all_gene_stat)

    # temp_cycle_edge_obvs = get_cycle_edge_obvs(tumor_name, all_gene_stat, edge_expression, prm_4, prm_5)
    # edge_state_list[[tumor_name]][,"DE_cof"] = temp_cycle_edge_obvs$DE_cof

  }
  return(edge_state_list)
}


#######################################################################################
# source("3_flux_subnet_add_gap.R")
subnet_edge_flux_list_main <- function(subnet_edge_expression, stat_all, deg_all, input_tumor_name, model) {
  print("subnet_edge_flux")
  #init

  tumors_array = c(input_tumor_name)
  subnet_edge_flux_list = get_edge_state_list(subnet_edge_expression, tumors_array, stat_all, deg_all, model)

  return(subnet_edge_flux_list)
}


#######################################################################################
# source("4_flux_cycle_add_gap.R")
cycle_edge_flux_list_main <- function(cycle_edge_expression, stat_all, deg_all, input_tumor_name, model) {
  print("cycle_edge_flux")
  #init


  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list = get_edge_state_list(cycle_edge_expression, tumors_array, stat_all, deg_all, model)

  return(cycle_edge_flux_list)
}



#######################################################################################
# source("5_flux_cycle_chain_list.R")

get_compoundsVer <- function(compoundsVer) {
  compoundsVer120921<-as.matrix(compoundsVer)
  for(i in 1:nrow(compoundsVer120921))
  {
    compoundsVer120921[i,1]<-unlist(strsplit(compoundsVer120921[i,4],";",fixed=T))[1]
  }
  rownames(compoundsVer120921)<-compoundsVer120921[,1]
  compoundsVer120921 = as.data.frame(compoundsVer120921)

  return(compoundsVer120921)
}

get_all_chain_list_cid<-function(tumor_name, cycle_edge_flux_list, cycle_directed) {
  #tumor_name = "COAD"
  cyc_metric = cycle_edge_flux_list[[tumor_name]]
  all_chain_cyc = list()
  for (i in 1:length(cycle_directed[,1])) {
    cycle_id = cycle_directed[i, "cycle_id"]
    tmp_cyc_metric = cyc_metric[which(cyc_metric$cycle_id == cycle_id),]
    tmp_cyc_chain = c()

    compound_chain = cycle_directed[i, "compound_chain"]
    compound_chain = unlist(strsplit(compound_chain, split = ";"))
    for (j in 1:length(compound_chain)) {
      cycle_node = unlist(strsplit(compound_chain[j], split = "->"))
      tmp_str = c(cycle_node[1])
      for (k in 1:(length(cycle_node)-1)) {
        tmp_cin = cycle_node[k]
        tmp_cout = cycle_node[k+1]
        tmp_state = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "state"]

        if(tmp_state == "gap") {
          tmp_str = c(tmp_str, "G")
        }else if(tmp_state == "up") {
          tmp_str = c(tmp_str, "U")
        }else {
          tmp_str = c(tmp_str, "")
        }

        tmp_str = c(tmp_str, cycle_node[k+1])
      }
      tmp_chain_str = paste(tmp_str, collapse = " -> ")
      tmp_cyc_chain = append(tmp_cyc_chain, tmp_chain_str)
    }

    all_chain_cyc[[i]] = tmp_cyc_chain
  }

  return(all_chain_cyc)
}

get_cycle_chain_list_main <- function(cycle_directed, cycle_edge_flux_list, input_tumor_name) {
  #init
  tumors_array = c(input_tumor_name)

  ## cid
  all_chain_list_cid = list()
  for (i in 1:length(tumors_array)) {
    temp_chain_list_cid =  get_all_chain_list_cid(tumors_array[i], cycle_edge_flux_list, cycle_directed)
    all_chain_list_cid[[tumors_array[i]]] = temp_chain_list_cid
  }

  return(all_chain_list_cid)
}

get_compounds_dict_main <- function() {
  #init
  library(readr)
  compoundsVer = read_csv(system.file("data/compoundsVer120921.csv", package = "CycleFlux"))

  compounds_dict = get_compoundsVer(compoundsVer)

  return(compounds_dict)
}


#######################################################################################
# source("6_degnode_gap_analysis.R")

get_never_considered_comp<-function(never_considered_comp, compounds_dict) {
  never_considered_comp_arr = never_considered_comp$x
  compounds_dict_cids = compounds_dict$ENTRY
  never_considered_comp_arr_cids = intersect(never_considered_comp_arr, compounds_dict_cids)
  never_considered_comp_names = never_considered_comp_arr_cids

  return(never_considered_comp_names)
}



get_ug_chain_list<-function(tumors_array, cycle_edge_flux_list, all_chain_list_cid) {
  ug_chain_list_dict = list()
  for (j in 1:length(tumors_array)) {
    tumor_name = tumors_array[j]
    temp_tumor_df = cycle_edge_flux_list[[tumor_name]]
    gap_c = unique(temp_tumor_df[which(temp_tumor_df$state == "gap"),'cycle_id'])
    up_c = unique(temp_tumor_df[which(temp_tumor_df$state == "up"),'cycle_id'])
    ug_c = intersect(gap_c, up_c)

    temp_chain_list = all_chain_list_cid[[tumor_name]]
    cycle_len = length(temp_chain_list)
    #names(temp_chain_list) = c(0:247)
    names(temp_chain_list) = c(1:cycle_len)
    ug_c = as.integer(ug_c)
    ug_chain_list = temp_chain_list[ug_c]

    ug_chain_list_dict[[tumor_name]] = ug_chain_list
  }
  return(ug_chain_list_dict)
}


single_cycle_gap_analysis_main <- function(cycle_edge_flux_list, all_chain_list_cid, input_tumor_name) {

  ##
  tumors_array = c(input_tumor_name)
  gapup_cycle_chain_list = get_ug_chain_list(tumors_array, cycle_edge_flux_list, all_chain_list_cid)
  return(gapup_cycle_chain_list)
}

never_considered_comp_names_main <- function(compounds_dict, input_tumor_name) {
  # init
  library(readr)
  never_considered_comp = read_csv(system.file("data/never_considered_comp.csv", package = "CycleFlux"))

  ##
  never_considered_comp_names = get_never_considered_comp(never_considered_comp, compounds_dict)
  return(never_considered_comp_names)
}


#######################################################################################
# source("7_degnode_add_gap.R")
cycle_edge_flux_list_in_main <- function(cycle_edgesucs_expression_in, stat_all, deg_all, input_tumor_name, model) {
  print("cycle_indegree_flux")

  #init
  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list_in = get_edge_state_list(cycle_edgesucs_expression_in, tumors_array, stat_all, deg_all, model)

  return(cycle_edge_flux_list_in)

}


cycle_edge_flux_list_out_main <- function(cycle_edgesucs_expression_out, stat_all, deg_all, input_tumor_name, model) {
  print("cycle_outdegree_flux")

  #init
  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list_out = get_edge_state_list(cycle_edgesucs_expression_out, tumors_array, stat_all, deg_all, model)

  return(cycle_edge_flux_list_out)
}



#######################################################################################

get_state_rids <- function(sps_net, gene_gapup_info, all_gene_stat, model, input_tumor_name) {
  rid_net = as.data.frame(sps_net[,c("Rid", "Gene_symbol")])
  colnames(rid_net) = c("rid", "gene_symbol")
  state_rids_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    if (model == 1) {
      temp_edge_flux = get_state_edges_model_1(tumor_name, rid_net, gene_gapup_info)
    }else if (model == 2) {
      temp_edge_flux = get_state_edges_model_2(tumor_name, rid_net, gene_gapup_info, all_gene_stat)
    }else {
      temp_edge_flux = get_state_edges_model_1(tumor_name, rid_net, gene_gapup_info)
    }
    state_rids_list[[tumor_name]] = temp_edge_flux
  }
  return(state_rids_list)
}

#######################################################################################
# source("8_my_ug_degnode_1.R")

#e2<--up--c2<--up--c1
get_ug_outdegree <- function(cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_ug_outdegree_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ##up
    out_degnode_df = cycle_edge_flux_list_out[[tumor_name]]
    ug_c = names(gapup_cycle_chain_list[[tumor_name]])
    out_degnode_df = out_degnode_df[which((out_degnode_df$state == "up") & (out_degnode_df$cycle_id %in% ug_c)),]

    select_ug_outdegree_list[[tumor_name]] = out_degnode_df
  }

  return(select_ug_outdegree_list)
}



get_ug_indegree <- function(cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_ug_indegree_list = list()

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ## up
    in_degnode_df = cycle_edge_flux_list_in[[tumor_name]]
    ug_c = names(gapup_cycle_chain_list[[tumor_name]])
    in_degnode_df = in_degnode_df[which((in_degnode_df$state == "up") & (in_degnode_df$cycle_id %in% ug_c)),]

    select_ug_indegree_list[[tumor_name]] = in_degnode_df
  }

  return(select_ug_indegree_list)
}

get_ug_cycle <- function(cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_ug_cycle_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ##up
    ug_cycle_df = cycle_edge_flux_list[[tumor_name]]
    ug_c = names(gapup_cycle_chain_list[[tumor_name]])
    ug_cycle_df = ug_cycle_df[which(ug_cycle_df$cycle_id %in% ug_c),]

    select_ug_cycle_list[[tumor_name]] = ug_cycle_df
  }

  return(select_ug_cycle_list)
}

#######################################################################################
# source("9_my_ug_degnode_3.R")

#e1<--up--(c2<--up--c1)<--up--e1

plot_up_cycle<-function(tumor_name, plot_name, cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names) {
  # library(igraph)
  # library(qgraph)
  select_shift_info = cycle_shift_path_df_list[[tumor_name]]
  select_in_upinfo = select_ug_indegree_list[[tumor_name]]
  select_out_upinfo = select_ug_outdegree_list[[tumor_name]]
  select_ug_cycle = select_ug_cycle_list[[tumor_name]]

  shift_cycle_ids = unique(select_shift_info$cycle_id)
  for (i in 1:length(shift_cycle_ids)) {
    cycid = shift_cycle_ids[i]
    ind_df = select_in_upinfo[which(select_in_upinfo$cycle_id == cycid),]
    od_df = select_out_upinfo[which(select_out_upinfo$cycle_id == cycid),]
    cycle_df = select_ug_cycle[which(select_ug_cycle$cycle_id == cycid),]
    shift_df = select_shift_info[which(select_shift_info$cycle_id == cycid),]

    ## nodes
    ind_node = unique(ind_df$ind)
    od_node = unique(od_df$od)
    cycle_node = union(cycle_df$c_in, cycle_df$c_out)
    all_node = unique(c(ind_node, od_node, cycle_node))

    g <- make_empty_graph(n = length(all_node))
    g <- set.vertex.attribute(g, "name", value=all_node)
    for (j in 1:length(all_node)) {
      if(all_node[j] %in% never_considered_comp_names) {
        V(g)[name==all_node[j]]$color <- "grey"
      }
    }

    ##edges
    for (j in 1:length(rownames(cycle_df))) {
      c_in = cycle_df[j, "c_in"]
      c_out = cycle_df[j, "c_out"]
      state = cycle_df[j, "state"]
      c_in_node = V(g)[name==c_in]
      c_out_node = V(g)[name==c_out]
      if(state == "gap") {
        g <- add_edges(g, c(c_in_node,c_out_node), color = "red", size = 3)
      }else if(state == "up") {
        g <- add_edges(g, c(c_in_node,c_out_node), color = "green", size = 3)
      }else {
        g <- add_edges(g, c(c_in_node,c_out_node), color = "black", size = 3)
      }
    }

    for (j in 1:length(rownames(od_df))) {
      node = od_df[j, "node"]
      od = od_df[j, "od"]
      v_node = V(g)[name==node]
      v_od = V(g)[name==od]
      if (!are.connected(g, v_node, v_od)) {
        g <- add_edges(g, c(v_node, v_od), color = "cyan", size = 1)
      }
    }

    for (j in 1:length(rownames(ind_df))) {
      node = ind_df[j, "node"]
      ind = ind_df[j, "ind"]
      v_node = V(g)[name==node]
      v_ind = V(g)[name==ind]
      if (!are.connected(g, v_ind, v_node)) {
        g <- add_edges(g, c(v_ind, v_node), color = "cyan", size = 1)
      }
    }

    ##shift
    if (length(rownames(shift_df)) > 0) {
      for (j in 1:length(rownames(shift_df))) {
        new_path = shift_df[j, "new_path"]
        new_path_nodes = unlist(strsplit(new_path, split = " -> "))
        if(length(new_path_nodes) == 5) {
          start_node = new_path_nodes[1]
          shift_node = new_path_nodes[3]
          end_node = new_path_nodes[5]

          V(g)[name==shift_node]$color <- "orchid"
          if (are.connected(g, start_node, shift_node)) {
            E(g)[get.edge.ids(g, c(start_node, shift_node))]$color = "orchid"
            E(g)[get.edge.ids(g, c(start_node, shift_node))]$size = 6
          }
          if (are.connected(g, shift_node, end_node)) {
            E(g)[get.edge.ids(g, c(shift_node, end_node))]$color = "orchid"
            E(g)[get.edge.ids(g, c(shift_node, end_node))]$size = 6
          }
        }


      }

    }



    #plot
    tmp_graph_name = paste0(tumor_name, "_C", cycid, ".png")
    tmp_plot_name = file.path(plot_name, tmp_graph_name)

    png_size = ceiling(vcount(g)/8)
    png_size[png_size < 15] = 15

    png(tmp_plot_name, height=png_size, width=png_size, units="in", res=100)
    e <- get.edgelist(g,names=FALSE)
    l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                           area=(vcount(g)^2),repulse.rad=(vcount(g)^3))
    plot(g, layout=l, vertex.size=5, edge.arrow.size=0.8, edge.width=E(g)$size)

    dev.off()
  }

}

plot_single_cycle <- function(cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names, input_tumor_name) {
  library(igraph)
  library(qgraph)
  ##select_ug_indegree_list select_ug_outdegree_list
  tumors_array = c(input_tumor_name)

  ##2.each_edge
  res_sub_path = "single_graphs"
  dir.create(file.path(res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_sub_path)

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    plot_up_cycle(tumor_name, res_file_path, cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names)
  }

}




#######################################################################################
# plot net cycle new
plot_net_up_cycle<-function(tumor_name, plot_name, cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names) {
  # library(igraph)
  # library(qgraph)
  select_shift_info = cycle_shift_path_df_list[[tumor_name]]
  shift_cycle_ids = unique(select_shift_info$cycle_id)

  select_ug_cycle = select_ug_cycle_list[[tumor_name]]
  select_ug_cycle = select_ug_cycle[which(select_ug_cycle$cycle_id %in% shift_cycle_ids),]

  ## nodes
  shf_node = unique(select_shift_info$shift_node)
  cycle_node = union(select_ug_cycle$c_in, select_ug_cycle$c_out)
  all_node = unique(c(shf_node, cycle_node))

  g <- make_empty_graph(n = length(all_node))
  g <- set.vertex.attribute(g, "name", value=all_node)
  for (j in 1:length(all_node)) {
    if(all_node[j] %in% never_considered_comp_names) {
      V(g)[name==all_node[j]]$color <- "grey"
    }
  }

  for (i in 1:length(shift_cycle_ids)) {
    cycid = shift_cycle_ids[i]
    cycle_df = select_ug_cycle[which(select_ug_cycle$cycle_id == cycid),]
    shift_df = select_shift_info[which(select_shift_info$cycle_id == cycid),]

    ##edges
    for (j in 1:length(rownames(cycle_df))) {
      c_in = cycle_df[j, "c_in"]
      c_out = cycle_df[j, "c_out"]
      state = cycle_df[j, "state"]
      c_in_node = V(g)[name==c_in]
      c_out_node = V(g)[name==c_out]
      if (!are.connected(g, c_in_node, c_out_node)) {
        if(state == "gap") {
          g <- add_edges(g, c(c_in_node,c_out_node), color = "red", size = 3)
        }else if(state == "up") {
          g <- add_edges(g, c(c_in_node,c_out_node), color = "green", size = 3)
        }else {
          g <- add_edges(g, c(c_in_node,c_out_node), color = "black", size = 3)
        }
      }

    }

    ##shift
    if (length(rownames(shift_df)) > 0) {
      for (j in 1:length(rownames(shift_df))) {
        new_path = shift_df[j, "new_path"]
        new_path_nodes = unlist(strsplit(new_path, split = " -> "))
        start_node = new_path_nodes[1]
        shift_node = new_path_nodes[3]
        end_node = new_path_nodes[5]

        v_start_node = as.integer(V(g)[name==start_node])
        v_shift_node = as.integer(V(g)[name==shift_node])
        v_end_node = as.integer(V(g)[name==end_node])

        V(g)[name==shift_node]$color <- "orchid"
        if (!are.connected(g, v_start_node, v_shift_node)) {
          g <- add_edges(g, c(v_start_node, v_shift_node), color = "orchid", size = 6)
        }
        if (!are.connected(g, v_shift_node, v_end_node)) {
          g <- add_edges(g, c(v_shift_node, v_end_node), color = "orchid", size = 6)
        }

      }

    }

  }

  #plot
  tmp_graph_name = paste0(tumor_name, "_net_shift", ".png")
  tmp_plot_name = file.path(plot_name, tmp_graph_name)

  png_size = ceiling(vcount(g)/8)
  png_size[png_size < 15] = 15

  png(tmp_plot_name, height=png_size, width=png_size, units="in", res=100)
  e <- get.edgelist(g,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                          area=(18*vcount(g)^2),repulse.rad=(vcount(g)^3.2))
  #l = layout_nicely(g)
  #l = layout_with_kk(g)
  #l = layout_with_graphopt(g)

  plot(g, layout=l, vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)

  dev.off()
}



plot_net_cycle <- function(cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names, input_tumor_name) {

  ##select_ug_indegree_list select_ug_outdegree_list
  tumors_array = c(input_tumor_name)
  library(igraph)
  library(qgraph)
  ##2.each_edge
  res_sub_path = "net_graphs"
  dir.create(file.path(res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_sub_path)

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    plot_net_up_cycle(tumor_name, res_file_path, cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names)
  }

}



#######################################################################################
# source("11_cyc_shift.R")
##
#Another way

get_tmp_comm_degnode_df <- function(cycle_id, permutation_node, ug_chain, comm_degnode) {
  tmp_comm_degnode_df = data.frame()

  for (p in 1:length(rownames(permutation_node))) {
    new_path_str = paste0(permutation_node[p,"od"], " -> U -> ", comm_degnode, " -> U -> ", permutation_node[p,"ind"])
    old_path_arr = c()
    for (c in 1:length(ug_chain)) {
      temp_chain_str = ug_chain[c]
      temp_chain_arr = unlist(strsplit(temp_chain_str,split = " -> "))

      od_index = which(temp_chain_arr == permutation_node[p,"od"])
      ind_index = which(temp_chain_arr == permutation_node[p,"ind"])
      if(length(od_index) > 1) { od_index = od_index[1] }
      if(length(ind_index) > 1) { ind_index = ind_index[1] }

      index_1 = min(od_index, ind_index)
      index_2 = max(od_index, ind_index)

      temp_old_path_arr_1 = temp_chain_arr[index_1:index_2]
      if("G" %in% temp_old_path_arr_1) {
        temp_old_path_str_1 = paste(temp_old_path_arr_1, collapse  = " -> ")
        old_path_arr = append(old_path_arr, temp_old_path_str_1)
      }

      chain_part_1 = temp_chain_arr[1:index_1]
      chain_part_2 = temp_chain_arr[index_2:(length(temp_chain_arr)-1)]
      temp_old_path_arr_2 = c(chain_part_2, chain_part_1)
      if("G" %in% temp_old_path_arr_2) {
        temp_old_path_str_2 = paste(temp_old_path_arr_2, collapse  = " -> ")
        old_path_arr = append(old_path_arr, temp_old_path_str_2)
      }

    }

    old_path_arr = unique(old_path_arr)
    old_path_str = paste(old_path_arr, collapse  = " ; ")

    tmp_comm_degnode_df[p, "cycle_id"] = cycle_id
    tmp_comm_degnode_df[p, "new_path"] = new_path_str
    tmp_comm_degnode_df[p, "old_path"] = old_path_str
    tmp_comm_degnode_df[p, "endpoint_node"] = paste0(permutation_node[p,"od"], ";", permutation_node[p,"ind"])
    tmp_comm_degnode_df[p, "shift_node"] = comm_degnode
  }

  return(tmp_comm_degnode_df)
}

get_tmp_cycle_shift_df <- function(cycle_id, comm_deg_arr, up_indeg, up_outdeg, ug_chain) {
  tmp_cycle_shift_df = data.frame()

  if(length(comm_deg_arr) != 0) {
    for (j in 1:length(comm_deg_arr)) {
      comm_degnode = comm_deg_arr[j]
      start_cyclenode = up_indeg[which(up_indeg$ind == comm_degnode), "node"] #hnode
      end_cyclenode = up_outdeg[which(up_outdeg$od == comm_degnode), "node"] #tnode

      cond1 = !identical(start_cyclenode, end_cyclenode)
      cond2 = (identical(start_cyclenode, end_cyclenode) & !((length(start_cyclenode)==1) & (length(end_cyclenode)==1)))
      if ( cond1 | cond2 ) {
        permutation_node = data.frame()
        for (od in end_cyclenode) {
          for (ind in start_cyclenode) {
            temp_node_df = data.frame()
            temp_node_df[1, "od"] = od
            temp_node_df[1, "ind"] = ind
            permutation_node = rbind(permutation_node, temp_node_df)
          }
        }

        if (length(which(permutation_node$od == permutation_node$ind)) > 0) {
          permutation_node = permutation_node[-which(permutation_node$od == permutation_node$ind),]
        }

        tmp_comm_degnode_df = get_tmp_comm_degnode_df(cycle_id, permutation_node, ug_chain, comm_degnode)
        tmp_cycle_shift_df = rbind(tmp_cycle_shift_df, tmp_comm_degnode_df)
      }
    }
  }

  return(tmp_cycle_shift_df)
}


get_cyc_shift_formula <- function(tumor_name, select_ug_indegree_list, select_ug_outdegree_list, gapup_cycle_chain_list) {
  ug_indegree = select_ug_indegree_list[[tumor_name]]
  ug_outdegree = select_ug_outdegree_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_chain_list = gapup_cycle_chain_list[[tumor_name]]

  all_cycle_shift_df = data.frame()
  for (i in 1:length(ug_c)) {
    cycle_id = ug_c[i]
    ug_chain = ug_chain_list[[i]]
    up_indeg = ug_indegree[which(ug_indegree$cycle_id == cycle_id),]
    up_outdeg = ug_outdegree[which(ug_outdegree$cycle_id == cycle_id),]

    up_indegnodes = up_indeg[,"ind"]
    up_outdegnodes = up_outdeg[,"od"]

    comm_deg_arr = intersect(up_indegnodes, up_outdegnodes)

    tmp_cycle_shift_df = get_tmp_cycle_shift_df(cycle_id, comm_deg_arr, up_indeg, up_outdeg, ug_chain)
    all_cycle_shift_df = rbind(all_cycle_shift_df, tmp_cycle_shift_df)
  }

  return(all_cycle_shift_df)
}


cyc_shift_main <- function(select_ug_indegree_list, select_ug_outdegree_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  cycle_shift_path_df_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    all_cycle_shift_df = get_cyc_shift_formula(tumor_name, select_ug_indegree_list, select_ug_outdegree_list, gapup_cycle_chain_list)
    cycle_shift_path_df_list[[tumor_name]] = all_cycle_shift_df
  }

  return(cycle_shift_path_df_list)
}



#######################################################################################
get_cycle_stat_list <- function(cycle_edge_flux_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_shift_path_df_list, input_tumor_name) {

  tumors_array = c(input_tumor_name)
  cycle_stat_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_edge_flux = cycle_edge_flux_list[[tumor_name]]
    cycle_edge_flux_in = cycle_edge_flux_list_in[[tumor_name]]
    cycle_edge_flux_out = cycle_edge_flux_list_out[[tumor_name]]
    cycle_shift_path = cycle_shift_path_df_list[[tumor_name]]

    cycids = unique(cycle_edge_flux$cycle_id)
    cycle_stat = data.frame()
    for (j in 1:length(cycids)) {
      cycle_id = cycids[j]
      cycle_df = cycle_edge_flux[which(cycle_edge_flux$cycle_id == cycle_id),]
      cycle_in_df = cycle_edge_flux_in[which(cycle_edge_flux_in$cycle_id == cycle_id),]
      cycle_out_df = cycle_edge_flux_out[which(cycle_edge_flux_out$cycle_id == cycle_id),]
      cycle_shift_df = cycle_shift_path[which(cycle_shift_path$cycle_id == cycle_id),]

      gap_edge = length(which(cycle_df$state == "gap"))
      up_edge = length(which(cycle_df$state == "up"))
      no_express_edge = length(which(cycle_df$state == "no_express"))
      no_change_edge = length(which(cycle_df$state == "no_change"))

      gap_indegree = length(which(cycle_in_df$state == "gap"))
      up_indegree = length(which(cycle_in_df$state == "up"))

      gap_outdegree = length(which(cycle_out_df$state == "gap"))
      up_outdegree = length(which(cycle_out_df$state == "up"))

      shift_edge = length(cycle_shift_df$new_path)
      ##
      cycle_stat[j, "cycle_id"] = cycle_id
      cycle_stat[j, "gap_edge"] = gap_edge
      cycle_stat[j, "up_edge"] = up_edge
      cycle_stat[j, "no_express_edge"] = no_express_edge
      cycle_stat[j, "no_change_edge"] = no_change_edge
      cycle_stat[j, "gap_indegree"] = gap_indegree
      cycle_stat[j, "up_indegree"] = up_indegree
      cycle_stat[j, "gap_outdegree"] = gap_outdegree
      cycle_stat[j, "up_outdegree"] = up_outdegree
      cycle_stat[j, "shift_edge"] = shift_edge
    }

    cycle_stat_list[[tumor_name]] = cycle_stat
  }

  return(cycle_stat_list)
}


#######################################################################################
# source("12_freq_stat.R")
#
# get_gap_cycleid_freq_list <- function(input_tumor_name, cycle_edge_flux_list) {
#   gap_cycleid_freq_list = list()
#   for (i in 1:length(input_tumor_name)) {
#     tumor_name = input_tumor_name[i]
#     gap_cycleid = cycle_edge_flux_list[[tumor_name]][which(cycle_edge_flux_list[[tumor_name]]$state == "gap"),]
#     gap_cycleid = unique(gap_cycleid)
#     gap_cycleid_freq = as.data.frame(table(gap_cycleid$cycle_id))
#
#     gap_cycleid_freq = gap_cycleid_freq[order(gap_cycleid_freq$Freq, decreasing = T),]
#
#     colnames(gap_cycleid_freq) = c("cycle_id", "gap_freq")
#     gap_cycleid_freq_list[[tumor_name]] = gap_cycleid_freq
#   }
#   return(gap_cycleid_freq_list)
# }

get_shift_node_freq_list <- function(input_tumor_name, cycle_shift_path_df_list) {
  shift_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    shift_node_freq = as.data.frame(table(cycle_shift_path_df_list[[tumor_name]]$shift_node))
    shift_node_freq = shift_node_freq[order(shift_node_freq$Freq, decreasing = T),]

    shift_node_df = cycle_shift_path_df_list[[tumor_name]]
    if (length(rownames(shift_node_df)) != 0) {
      cycle_ids =  aggregate(shift_node_df$cycle_id, list(shift_node_df$shift_node), paste, collapse = ";")
      rownames(cycle_ids) = cycle_ids[,1]

      colnames(shift_node_freq) = c("node", "shift_freq")

      for (j in 1:length(shift_node_freq[,1])) {
        shift_node_freq[j,"cycle_id"] = cycle_ids[shift_node_freq[j,"node"],2]
      }
    }

    shift_node_freq_list[[tumor_name]] = shift_node_freq
  }
  return(shift_node_freq_list)
}

get_shift_edge_node_freq_list <- function(input_tumor_name, cycle_shift_path_df_list) {
  shift_edge_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    temp_cycle_shift_path_df = cycle_shift_path_df_list[[tumor_name]]
    shift_edgenode_freq = data.frame()
    if (length(rownames(temp_cycle_shift_path_df)) != 0) {
      shift_edgenode_freq = as.data.frame(table(temp_cycle_shift_path_df[,c("endpoint_node","shift_node")]))
      shift_edgenode_freq = shift_edgenode_freq[order(shift_edgenode_freq$Freq, decreasing = T),]

      colnames(shift_edgenode_freq) = c("edge_node", "node", "shift_freq")

      shift_edgenode_freq = shift_edgenode_freq[which(shift_edgenode_freq$shift_freq != 0),]
    }

    shift_edge_node_freq_list[[tumor_name]] = shift_edgenode_freq
  }
  return(shift_edge_node_freq_list)
}


convert_deg <- function(deg) {
  deg_dict = data.frame(
    "state_value" = c("no_express", "up", "gap", "no_change")
  )
  rownames(deg_dict) = c("-10", "1", "-1", "0")

  deg_val = as.character(deg)
  deg_names = names(deg)
  deg_state_val = deg_dict[deg_val, "state_value"]
  names(deg_state_val) = deg_names

  return(deg_state_val)
}


loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#######################################################################################
#' getCycleFlux function
#'
#' This function for cycle shift situation
#'
#' @param
#' model:
#' IF model=1
#' For an edge, if any gene on the edge is up, the edge is up, if the edge contains both the gap gene and the up gene, it is regarded as a gap.
#' IF model=2
#' For an edge, if the gene with the largest mean of tumor value is up, then the edge is up
#' @param
#' single_graph:
#' If this parameter is TRUE, an image will be generated in this directory for each cycle that matches the shift definition.
#' @param
#' net_graph:
#' If this parameter is TRUE, an image containing all cycles that match the shift definition will be generated in this directory.
#' Default:
#' model=1
#' single_graph=TRUE
#' net_graph=TRUE
#' @keywords getCycleFlux
#' @export
#' @examples
#' getCycleFlux(input_net_file, input_stat_gene_file, input_deg_gene_file, 1, TRUE, FALSE)
#'
getCycleFlux <- function(input_net_file, input_stat_gene_file, input_deg_gene_file, model=1, single_graph=TRUE, net_graph=TRUE) {
  #source("1_cycle_directed.R")

  ##
  hsa_net = loadRData(file.path(input_net_file))
  g = build_net(hsa_net)
  cycle_directed = get_cycle_directed(hsa_net, g)
  cycle_edge_expression = get_cycle_edge_expression(cycle_directed, g)
  cycle_directed = rm_disuse_cycle_directed(cycle_edge_expression, cycle_directed)
  cycle_edge_expression = reid_cycle_edge_expression(cycle_edge_expression)

  subnet_edge_expression = get_subnet_edge_expression(g)
  cycle_expression = get_cycle_expression(cycle_directed, g)
  cycle_edgesucs_expression_in = get_cycle_edgesucs_expression_in(cycle_directed, g)
  cycle_edgesucs_expression_out = get_cycle_edgesucs_expression_out(cycle_directed, g)

  ##
  stat_all = loadRData(file.path(input_stat_gene_file))
  deg_all = loadRData(file.path(input_deg_gene_file))
  deg_all = lapply(deg_all, convert_deg)

  input_tumor_name = names(deg_all)

  #source("2_get_missing_gene.R")
  #gene_missing_list = get_gene_and_tofind_list_main(subnet_edge_expression, stat_all, deg_all)
  #source("3_flux_subnet_add_gap.R")
  subnet_edge_flux_list = subnet_edge_flux_list_main(subnet_edge_expression, stat_all, deg_all, input_tumor_name, model)
  #source("4_flux_cycle_add_gap.R")
  cycle_edge_flux_list = cycle_edge_flux_list_main(cycle_edge_expression, stat_all, deg_all, input_tumor_name, model)

  #source("5_flux_cycle_chain_list.R")
  all_chain_list_cid = get_cycle_chain_list_main(cycle_directed, cycle_edge_flux_list, input_tumor_name)
  compounds_dict = get_compounds_dict_main()

  #source("6_degnode_gap_analysis.R")
  gapup_cycle_chain_list = single_cycle_gap_analysis_main(cycle_edge_flux_list, all_chain_list_cid, input_tumor_name)
  never_considered_comp_names = never_considered_comp_names_main(compounds_dict, input_tumor_name)

  #source("7_degnode_add_gap.R")
  cycle_edge_flux_list_in = cycle_edge_flux_list_in_main(cycle_edgesucs_expression_in, stat_all, deg_all, input_tumor_name, model)
  cycle_edge_flux_list_out = cycle_edge_flux_list_out_main(cycle_edgesucs_expression_out, stat_all, deg_all, input_tumor_name, model)

  #source("8_my_ug_degnode_1.R")
  select_ug_indegree_list = get_ug_indegree(cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)
  select_ug_outdegree_list = get_ug_outdegree(cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)
  select_ug_cycle_list = get_ug_cycle(cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)

  #### result ####
  #source("10_cyc_classification.R")

  #source("11_cyc_shift.R")
  cycle_shift_path_df_list = cyc_shift_main(select_ug_indegree_list, select_ug_outdegree_list, gapup_cycle_chain_list, input_tumor_name)

  #source("9_my_ug_degnode_3.R")
  if (single_graph) {
    plot_single_cycle(cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names, input_tumor_name)
  }

  if (net_graph) {
    plot_net_cycle(cycle_shift_path_df_list, select_ug_indegree_list, select_ug_outdegree_list, select_ug_cycle_list, never_considered_comp_names, input_tumor_name)
  }

  ##
  rids_state_list = get_state_rids(hsa_net, deg_all, stat_all, model, input_tumor_name)
  cycle_stat_list = get_cycle_stat_list(cycle_edge_flux_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_shift_path_df_list, input_tumor_name)

  #source("12_freq_stat.R")
  # gap_cycleid_freq_list = get_gap_cycleid_freq_list(input_tumor_name, cycle_edge_flux_list)
  shift_node_freq_list = get_shift_node_freq_list(input_tumor_name, cycle_shift_path_df_list)
  shift_edge_node_freq_list = get_shift_edge_node_freq_list(input_tumor_name, cycle_shift_path_df_list)

  result_list = list()
  result_list[["cycle_edge"]] = cycle_edge_flux_list
  result_list[["shift_edge"]] = cycle_shift_path_df_list
  result_list[["cycle_chain"]] = all_chain_list_cid
  result_list[["cycle_stat"]] = cycle_stat_list
  result_list[["rids_state"]] = rids_state_list
  result_list[["shift_edge_freq"]] = shift_node_freq_list
  result_list[["shift_edge_nfreq"]] = shift_edge_node_freq_list

  return(result_list)
}







