


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
    print(paste0("cyc_", i))
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
    print(paste0("cyc_", i))
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





get_cyc_gene_not_in_tgca<-function(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1){
  cyc_gene_kegg = c()
  cyc_gene_tcga = c()
  for (i in 1:length(rownames(to_find_correct_symbol_v1_1))) {
    kegg_gene = to_find_correct_symbol_v1_1[i,1]
    tcga_gene = to_find_correct_symbol_v1_1[i,2]
    if(kegg_gene %in% cyc_gene_expressinfo){
      cyc_gene_kegg = append(cyc_gene_kegg,kegg_gene)
    }
    if(tcga_gene %in% cyc_gene_expressinfo){
      cyc_gene_tcga = append(cyc_gene_tcga,tcga_gene)
    }
  }

  cyc_gene_not_in_tgca = c()
  cyc_gene_in_tgca = c()
  tcga_gene_list = rownames(all_gene_stat[[1]])
  cyc_gene_expressinfo = unique(cyc_gene_expressinfo)
  for (i in 1:length(cyc_gene_expressinfo)) {
    gene = cyc_gene_expressinfo[i]
    if(gene %in% tcga_gene_list){
      cyc_gene_in_tgca = append(cyc_gene_in_tgca,gene)
    }else{
      cyc_gene_not_in_tgca = append(cyc_gene_not_in_tgca,gene)
    }
  }

  return(cyc_gene_not_in_tgca)
}


get_gene_missing_list<-function(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1){
  #test
  cyc_gene_not_in_tgca = get_cyc_gene_not_in_tgca(cyc_gene_expressinfo, all_gene_stat, to_find_correct_symbol_v1_1)
  cyc_gene_not_in_tgca_but_in_tofind = intersect(cyc_gene_not_in_tgca,to_find_correct_symbol_v1_1$X1)
  cyc_gene_not_in_tgca_not_in_tofind = setdiff(cyc_gene_not_in_tgca, cyc_gene_not_in_tgca_but_in_tofind)

  #result
  gene_missing_list = list()
  gene_missing_list[["cyc_gene_expressinfo"]] = cyc_gene_expressinfo
  gene_missing_list[["cyc_gene_not_in_tgca"]] = cyc_gene_not_in_tgca
  gene_missing_list[["cyc_gene_not_in_tgca_but_in_tofind"]] = cyc_gene_not_in_tgca_but_in_tofind
  gene_missing_list[["cyc_gene_not_in_tgca_not_in_tofind"]] = cyc_gene_not_in_tgca_not_in_tofind

  return(gene_missing_list)
}





get_gene_and_tofind_list_main <- function(subnet_edge_expression, stat_all, deg_all) {


  cyc_gene_expressinfo = get_gene_expressinfo(subnet_edge_expression)

  library(readr)
  to_find_correct_symbol_v1_1 = read_csv(system.file("data/to_find_correct_symbol_v1.1.csv", package = "CycleFlux"), col_names = FALSE)

  gene_missing_list = get_gene_missing_list(cyc_gene_expressinfo, stat_all, to_find_correct_symbol_v1_1)

  return(gene_missing_list)
}

#######################################################################################
# function for gap
# get_mean_fc<-function(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list) {
#   missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
#
#   temp_cycle_edge_flux = cycle_edge_flux_list[[tumor_name]]
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
#       if (!(temp_gene %in% missing_genes)) {
#         genesign = as.double(temp_gene_stat[temp_gene,"FC"])
#         temp_edge_sum_fc = temp_edge_sum_fc + genesign
#       }
#     }
#     temp_cycle_edge_flux[i,"mean_fc"] = temp_edge_sum_fc/length(temp_gene_symbol_arr) #more equal
#   }
#
#   return(temp_cycle_edge_flux)
# }

get_state_edges <- function(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list, prm_1) {

  if (prm_1 == 1) {
    temp_edge_flux = get_state_edges_model_1(tumor_name, cycle_edge_flux_list, gene_gapup_info, gene_missing_list)
  }else if (prm_1 == 2) {
    temp_edge_flux = get_state_edges_model_2(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list)
  }else {
    temp_edge_flux = get_state_edges_model_1(tumor_name, cycle_edge_flux_list, gene_gapup_info, gene_missing_list)
  }
  return(temp_edge_flux)
}

get_state_edges_model_1 <- function(tumor_name, cycle_edge_flux_list, gene_gapup_info, gene_missing_list) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  temp_gene_gapup_info = gene_gapup_info[[tumor_name]]
  temp_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_edge_flux[, "state"] = "normal"
  for (i in 1:length(temp_edge_flux[,1])) {
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr = unlist(strsplit(gene_arr, split = ";"))
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if (!(temp_gene %in% missing_genes)) {
        temp_edge_flux[i, "state"] = temp_gene_gapup_info[[temp_gene]]
      }
    }
  }
  return(temp_edge_flux)
}

get_state_edges_model_2 <- function(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]

  temp_gene_gapup_info = gene_gapup_info[[tumor_name]]
  temp_all_gene_stat = all_gene_stat[[tumor_name]]
  temp_edge_flux = cycle_edge_flux_list[[tumor_name]]
  temp_edge_flux[, "state"] = "normal"
  for (i in 1:length(temp_edge_flux[,1])) {
    gene_arr = temp_edge_flux[i, "gene_symbol"]
    gene_arr <- unlist(strsplit(gene_arr, split = ";"))

    gene_mt_df = data.frame()
    ## find mt(mean tumor value) gene
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      if (temp_gene %in% missing_genes) {
        temp_mt = 0
      }else {
        temp_mt = temp_all_gene_stat[temp_gene, "mt"]
      }
      gene_mt_df[j, "gene"] = temp_gene
      gene_mt_df[j, "mt"] = temp_mt
    }

    gene_mt_df = gene_mt_df[order(gene_mt_df$mt,decreasing = TRUE),]
    max_mt_gene = gene_mt_df[1, "gene"]
    if (!(max_mt_gene %in% missing_genes)) {
      temp_edge_flux[i, "state"] = temp_gene_gapup_info[[max_mt_gene]]
    }
  }
  return(temp_edge_flux)
}

get_cycle_edge_obvs <- function(tumor_name, all_gene_stat, cycle_edge_expression, gene_missing_list, prm_1, prm_2) {
  missing_genes = gene_missing_list[["cyc_gene_not_in_tgca"]]
  temp_all_gene_stat = all_gene_stat[[tumor_name]]

  cycle_edge_obvs = data.frame()

  for (i in 1:length(cycle_edge_expression[,1])) {
    gene_str = cycle_edge_expression[i,"gene_symbol"]
    gene_arr = unlist(strsplit(gene_str,split = ";"))
    gene_num = 0
    obvs_num = 0
    for (j in 1:length(gene_arr)) {
      temp_gene = gene_arr[j]
      temp_pvalue = 0
      temp_fc = 0
      if (!(temp_gene %in% missing_genes)) {
        temp_pvalue = temp_all_gene_stat[temp_gene,"p.value"]
        temp_fc = temp_all_gene_stat[temp_gene,"FC"]
        gene_num = gene_num + 1
      }

      if ((temp_pvalue < prm_1) & (abs(temp_fc) > prm_2)) {
        obvs_num = obvs_num + 1
      }
    }
    obvs_rat = obvs_num/gene_num
    cycle_edge_obvs[i,"DE_cof"] = round(obvs_rat, 2)
  }
  return(cycle_edge_obvs)
}


## main function
get_cycle_edge_flux_list <- function(edge_expression, tumors_array, all_gene_stat, gene_gapup_info, gene_missing_list, prm_1, prm_2, prm_3) {
  cycle_edge_flux_list = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    cycle_edge_flux_list[[tumor_name]] = edge_expression

    cycle_edge_flux_list[[tumor_name]] = get_state_edges(tumor_name, cycle_edge_flux_list, gene_gapup_info, all_gene_stat, gene_missing_list, prm_3)
    #cycle_edge_flux_list[[tumor_name]] = get_mean_fc(tumor_name, cycle_edge_flux_list, all_gene_stat, gene_missing_list)

    temp_cycle_edge_obvs = get_cycle_edge_obvs(tumor_name, all_gene_stat, edge_expression, gene_missing_list, prm_1, prm_2)
    cycle_edge_flux_list[[tumor_name]][,"DE_cof"] = temp_cycle_edge_obvs$DE_cof

  }
  return(cycle_edge_flux_list)
}


#######################################################################################
# source("3_flux_subnet_add_gap.R")
subnet_edge_flux_list_main <- function(subnet_edge_expression, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3) {
  print("subnet_edge_flux")
  #init

  tumors_array = c(input_tumor_name)
  subnet_edge_flux_list = get_cycle_edge_flux_list(subnet_edge_expression, tumors_array, stat_all, deg_all, gene_missing_list, prm_1, prm_2, prm_3)

  return(subnet_edge_flux_list)
}


#######################################################################################
# source("4_flux_cycle_add_gap.R")
cycle_edge_flux_list_main <- function(cycle_edge_expression, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3) {
  print("cycle_edge_flux")
  #init


  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list = get_cycle_edge_flux_list(cycle_edge_expression, tumors_array, stat_all, deg_all, gene_missing_list, prm_1, prm_2, prm_3)

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
    #有gap现象的环的cycle_id
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
cycle_edge_flux_list_in_main <- function(cycle_edgesucs_expression_in, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3) {


  #init
  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list_in = get_cycle_edge_flux_list(cycle_edgesucs_expression_in, tumors_array, stat_all, deg_all, gene_missing_list, prm_1, prm_2, prm_3)

  return(cycle_edge_flux_list_in)

}


cycle_edge_flux_list_out_main <- function(cycle_edgesucs_expression_out, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3) {


  #init
  tumors_array = c(input_tumor_name)

  ##
  cycle_edge_flux_list_out = get_cycle_edge_flux_list(cycle_edgesucs_expression_out, tumors_array, stat_all, deg_all, gene_missing_list, prm_1, prm_2, prm_3)

  return(cycle_edge_flux_list_out)
}


#######################################################################################
# source("8_my_ug_degnode_1.R")

#e2<--up--c2<--up--c1
#outdegree
get_select_out_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  #tumor_name = "COAD"

  out_gapinfo = cycle_edge_flux_list_out[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycle_id %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$state == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  #  e1<----c2<--up--c1
  select_out_upinfo = data.frame()
  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycle_id = ug_cyc_pct_upnode[i, "cycle_id"]
    tmp_cout = ug_cyc_pct_upnode[i, "c_out"]

    tmp_row1 = out_gapinfo[which((out_gapinfo$cycle_id == tmp_cycle_id) & (out_gapinfo$node == tmp_cout)),]
    select_out_upinfo = rbind(select_out_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_out_upinfo = select_out_upinfo[which(select_out_upinfo$state == "up"),]
  return(select_out_upinfo)
}


get_select_in_upinfo_eachedge <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list) {
  in_gapinfo = cycle_edge_flux_list_in[[tumor_name]]
  cyc_pct = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_cyc_pct = cyc_pct[which(cyc_pct$cycle_id %in% ug_c),]

  #  c2<--up--c1
  #ug_cyc_pct_upnode = ug_cyc_pct[which(ug_cyc_pct$state == "up"),]
  ug_cyc_pct_upnode = ug_cyc_pct

  up_cin_node = ug_cyc_pct_upnode$c_in

  #  e1<----c2<--up--c1
  select_in_upinfo = data.frame()

  for (i in 1:length(ug_cyc_pct_upnode[,1])) {
    tmp_cycle_id = ug_cyc_pct_upnode[i, "cycle_id"]
    tmp_cin = ug_cyc_pct_upnode[i, "c_in"]

    tmp_row1 = in_gapinfo[which((in_gapinfo$cycle_id == tmp_cycle_id) & (in_gapinfo$node == tmp_cin)),]
    select_in_upinfo = rbind(select_in_upinfo, tmp_row1)
  }

  #  e1<--up--c2<--up--c1
  select_in_upinfo = select_in_upinfo[which(select_in_upinfo$state == "up"),]

  return(select_in_upinfo)
}

my_ug_degnode_1_main <- function(cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_out_upinfo_eachedge = list()
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ##up
    temp_select_out_upinfo_eachedge = get_select_out_upinfo_eachedge(tumor_name, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)

    select_out_upinfo_eachedge[[tumor_name]] = temp_select_out_upinfo_eachedge
  }

  return(select_out_upinfo_eachedge)
}



my_ug_degnode_2_main <- function(cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  tumors_array = c(input_tumor_name)
  select_in_upinfo_eachedge = list()

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    ## up
    temp_select_in_upinfo_eachedge = get_select_in_upinfo_eachedge(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list)

    select_in_upinfo_eachedge[[tumor_name]] = temp_select_in_upinfo_eachedge
  }

  return(select_in_upinfo_eachedge)
}


#######################################################################################
# source("9_my_ug_degnode_3.R")

#e1<--up--(c2<--up--c1)<--up--e1
# st2_cycle_id = select_st2_1_gapinfo_e1e2$cycle_id
# st3_cycle_id = select_st3_1_gapinfo_e1e2$cycle_id
# intersect_cycle_id = intersect(st2_cycle_id, st3_cycle_id)

plot_up_cycle<-function(tumor_name, plot_name, select_in_upinfo, select_out_upinfo, gapup_cycle_chain_list, never_considered_comp_names) {
  tmp_chain_str = gapup_cycle_chain_list[[tumor_name]]
  if (length(tmp_chain_str) > 0) {
    cyc_id_arr = names(tmp_chain_str)

    for (k in 1:length(cyc_id_arr)) {
      tmp_cyc_id = cyc_id_arr[k]
      tmp_cycle_chain_str = tmp_chain_str[[tmp_cyc_id]] #[1]C->U->c->G->c [2]c->G->c->U->c
      tmp_cycle_chain_arr = unlist(strsplit(tmp_cycle_chain_str[1],split = " -> ")) #'c' 'U' 'c' 'G' 'c'

      cyc_node_vector = c()
      for (i in 1:length(tmp_cycle_chain_arr)) {
        if((i%%2 == 1) & (i != length(tmp_cycle_chain_arr))) {
          cyc_node_vector = append(cyc_node_vector, tmp_cycle_chain_arr[i])
        }
      }


      ##
      up_out_node = select_out_upinfo[which(select_out_upinfo$cycle_id == tmp_cyc_id), "node"]
      up_out_od = select_out_upinfo[which(select_out_upinfo$cycle_id == tmp_cyc_id), "od"]
      up_in_node = select_in_upinfo[which(select_in_upinfo$cycle_id == tmp_cyc_id), "node"]
      up_in_ind = select_in_upinfo[which(select_in_upinfo$cycle_id == tmp_cyc_id), "ind"]

      #degree_node_arr = c(rbind(up_out_od, up_in_ind))
      degree_node_arr = c(up_out_od, up_in_ind)
      degree_node_arr = unique(degree_node_arr)

      library(igraph)
      all_node_vector = append(cyc_node_vector, degree_node_arr)
      all_node_vector = unique(all_node_vector)
      g <- make_empty_graph(n = length(all_node_vector))

      g <- set.vertex.attribute(g, "name", value=all_node_vector)

      for (i in 1:length(all_node_vector)) {
        if(all_node_vector[i] %in% never_considered_comp_names) {
          V(g)[name==all_node_vector[i]]$color <- "grey"
        }
      }

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

      # up_out_node
      # up_out_od
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

      tmp_graph_name = paste0(tumor_name, "_C", cyc_id_arr[k], ".png")
      tmp_plot_name = file.path(plot_name, tmp_graph_name)

      png(tmp_plot_name, 500, 500)
      plot(g)
      dev.off()
    }
  }
}

my_ug_degnode_3_main <- function(select_in_upinfo_eachedge, select_out_upinfo_eachedge, gapup_cycle_chain_list, never_considered_comp_names, input_tumor_name) {

  ##select_in_upinfo_eachedge select_out_upinfo_eachedge
  tumors_array = c(input_tumor_name)

  ##2.each_edge
  res_sub_path = "single_graphs"
  dir.create(file.path(res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_sub_path)

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    plot_up_cycle(tumor_name, res_file_path, select_in_upinfo_eachedge[[tumor_name]], select_out_upinfo_eachedge[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

  }

}


#######################################################################################
# source("10_cyc_classification.R")
#
# each edge
#
# 1.Up only head in and tail out
# 2.There is a degree of up in the middle of Up
# 3.Another way
#
# 1.Up only head in and tail out
# Only 2 nodes in the ring have the degree of up
# And when these two nodes are used as c_in or c_out in this ring, state==up
#
# 2.There is a degree of up in the middle of Up
# There are >3 nodes in the ring with up degree
# And when these nodes are used as c_in or c_out in this ring, the state at both ends are not all gaps
#
# 3.Another way
# There is a degree node (ind or od)
# This node is the degree of 2 nodes in this ring at the same time
# And when these two nodes are used as c_in or c_out in this ring, state==up
# And the edge between these two nodes, there is a gap edge
#
# 4. Burst from 1 node
#
# 5.Other classes (with up out-degree in the middle of Gap)
#


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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) == 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_state = tmp_node_row[,"state"]
        if(!("up" %in% tmp_node_row[,"state"])){
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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) > 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_state = tmp_node_row[,"state"]
        if(("gap" %in% tmp_node_row[,"state"]) & !("up" %in% tmp_node_row[,"state"])){
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

    up_indeg = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),]
    up_outdeg = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),]

    up_indeg_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),"ind"]
    up_outdeg_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),"od"]

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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) == 1) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_state = tmp_node_row[,"state"]
        if(!("up" %in% tmp_node_row[,"state"])){
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



get_ug_class_df <- function(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  cyc_df = as.data.frame(cycle_edge_flux_list[[tumor_name]])

  ug_class_df = as.data.frame(unique(cyc_df[which(cyc_df$cycle_id %in% ug_c), "cycle_id"]))
  ug_class_df[,"classid"] = 0
  colnames(ug_class_df) = c("cycle_id", "classid")

  class_1_c = get_class_1_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_2_c = get_class_2_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_3_c = get_class_3_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
  class_4_c = get_class_4_c(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)

  class_1_c = setdiff(class_1_c, class_3_c)
  class_2_c = setdiff(class_2_c, class_3_c)

  class_ist = c(class_1_c, class_2_c, class_3_c, class_4_c)
  class_ist = unique(class_ist)

  class_5_c = setdiff(ug_c, class_ist)

  ug_class_df[which(ug_class_df$cycle_id %in% class_1_c), "classid"] = 1
  ug_class_df[which(ug_class_df$cycle_id %in% class_2_c), "classid"] = 2
  ug_class_df[which(ug_class_df$cycle_id %in% class_3_c), "classid"] = 3
  ug_class_df[which(ug_class_df$cycle_id %in% class_4_c), "classid"] = 4
  ug_class_df[which(ug_class_df$cycle_id %in% class_5_c), "classid"] = 5

  return(ug_class_df)
}



cyc_classification_main <- function(cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  cycle_upgap_class_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    ug_class_df = get_ug_class_df(tumor_name, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    cycle_upgap_class_list[[tumor_name]] = ug_class_df
  }

  return(cycle_upgap_class_list)
}



#######################################################################################
# source("11_cyc_shift.R")
##
#Another way
get_cyc_shift_formula <- function(tumor_name, cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list) {
  indeg_df = cycle_edge_flux_list_in[[tumor_name]]
  outdeg_df = cycle_edge_flux_list_out[[tumor_name]]
  cyc_df = cycle_edge_flux_list[[tumor_name]]
  ug_c = names(gapup_cycle_chain_list[[tumor_name]])
  ug_chain_list = gapup_cycle_chain_list[[tumor_name]]
  new_old_path_df = data.frame()

  for (i in 1:length(ug_c)) {
    cycle_id = ug_c[i]
    print(cycle_id)
    tmp_cyc_indeg = indeg_df[which(indeg_df$cycle_id == ug_c[i]),]
    tmp_cyc_outdeg = outdeg_df[which(outdeg_df$cycle_id == ug_c[i]),]
    tmp_cyc = cyc_df[which(cyc_df$cycle_id == ug_c[i]),]

    up_indeg = tmp_cyc_indeg[which(tmp_cyc_indeg$state == "up"),]
    up_outdeg = tmp_cyc_outdeg[which(tmp_cyc_outdeg$state == "up"),]

    up_indeg_arr = up_indeg[,"ind"]
    up_outdeg_arr = up_outdeg[,"od"]

    ist_deg_arr = intersect(up_indeg_arr, up_outdeg_arr)
    if(length(ist_deg_arr) != 0) {
      #shared_deg_arr = c()
      for (j in 1:length(ist_deg_arr)) {
        indegnode = up_indeg[which(up_indeg$ind == ist_deg_arr[j]), "node"] #hnode
        outdegnode = up_outdeg[which(up_outdeg$od == ist_deg_arr[j]), "node"] #tnode

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

          new_path_str = paste0(permutation_node[1,"od"], " -> U -> ", ist_deg_arr[j], " -> U -> ", permutation_node[1,"ind"])

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
          temp_old_path_df = data.frame()
          temp_old_path_df[1, "cycle_id"] = cycle_id
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


cyc_shift_main <- function(cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name) {

  # init
  cycle_flux_path_df_list = list()
  tumors_array = c(input_tumor_name)
  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]
    new_old_path_df = get_cyc_shift_formula(tumor_name, cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list)
    cycle_flux_path_df_list[[tumor_name]] = new_old_path_df
  }

  return(cycle_flux_path_df_list)
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

get_shift_node_freq_list <- function(input_tumor_name, cycle_flux_path_df_list) {
  shift_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    shift_node_freq = as.data.frame(table(cycle_flux_path_df_list[[tumor_name]]$shift_node))
    shift_node_freq = shift_node_freq[order(shift_node_freq$Freq, decreasing = T),]

    shift_node_df = cycle_flux_path_df_list[[tumor_name]]
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

get_shift_edge_node_freq_list <- function(input_tumor_name, cycle_flux_path_df_list) {
  shift_edge_node_freq_list = list()
  for (i in 1:length(input_tumor_name)) {
    tumor_name = input_tumor_name[i]
    temp_cycle_flux_path_df = cycle_flux_path_df_list[[tumor_name]]
    shift_edgenode_freq = data.frame()
    if (length(rownames(temp_cycle_flux_path_df)) != 0) {
      shift_edgenode_freq = as.data.frame(table(temp_cycle_flux_path_df[,c("endpoint_node","shift_node")]))
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
    "state_value" = c("no_express", "up", "gap", "normal")
  )
  rownames(deg_dict) = c("00", "01", "10", "11")

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
#' prm_1&prm_2:
#' FOR DE_cof
#' (gene pvalue < prm_1) & (abs(gene fc) > prm_2)
#' -> DE num ++
#' -> DE_cof = DE num/gene num
#' @param
#' prm_3:
#' IF prm_3=1
#' For an edge, if any gene on the edge is up, the edge is up
#' IF prm_3=2
#' For an edge, if the gene with the largest mean of tumor value is up, then the edge is up
#' @param
#' Default:
#' prm_1=0.05
#' prm_2=1
#' prm_3=1
#' @keywords getCycleFlux
#' @export
#' @examples
#' getCycleFlux(input_net_file, input_stat_gene_file, input_deg_gene_file, 0.05, 1, 1)
#'
getCycleFlux <- function(input_net_file, input_stat_gene_file, input_deg_gene_file, prm_1=0.05, prm_2=1, prm_3=1) {
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
  gene_missing_list = get_gene_and_tofind_list_main(subnet_edge_expression, stat_all, deg_all)
  #source("3_flux_subnet_add_gap.R")
  subnet_edge_flux_list = subnet_edge_flux_list_main(subnet_edge_expression, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3)
  #source("4_flux_cycle_add_gap.R")
  cycle_edge_flux_list = cycle_edge_flux_list_main(cycle_edge_expression, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3)

  #source("5_flux_cycle_chain_list.R")
  all_chain_list_cid = get_cycle_chain_list_main(cycle_directed, cycle_edge_flux_list, input_tumor_name)
  compounds_dict = get_compounds_dict_main()

  #source("6_degnode_gap_analysis.R")
  gapup_cycle_chain_list = single_cycle_gap_analysis_main(cycle_edge_flux_list, all_chain_list_cid, input_tumor_name)
  never_considered_comp_names = never_considered_comp_names_main(compounds_dict, input_tumor_name)

  #source("7_degnode_add_gap.R")
  cycle_edge_flux_list_in = cycle_edge_flux_list_in_main(cycle_edgesucs_expression_in, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3)
  cycle_edge_flux_list_out = cycle_edge_flux_list_out_main(cycle_edgesucs_expression_out, gene_missing_list, stat_all, deg_all, input_tumor_name, prm_1, prm_2, prm_3)

  #source("8_my_ug_degnode_1.R")
  select_out_upinfo_eachedge = my_ug_degnode_1_main(cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)
  select_in_upinfo_eachedge = my_ug_degnode_2_main(cycle_edge_flux_list_in, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)

  #### result ####

  #source("9_my_ug_degnode_3.R")
  my_ug_degnode_3_main(select_in_upinfo_eachedge, select_out_upinfo_eachedge, gapup_cycle_chain_list, never_considered_comp_names, input_tumor_name)

  #source("10_cyc_classification.R")
  cycle_upgap_class_list = cyc_classification_main(cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)

  #source("11_cyc_shift.R")
  cycle_flux_path_df_list = cyc_shift_main(cycle_upgap_class_list, cycle_edge_flux_list_in, cycle_edge_flux_list_out, cycle_edge_flux_list, gapup_cycle_chain_list, input_tumor_name)

  #source("12_freq_stat.R")
  # gap_cycleid_freq_list = get_gap_cycleid_freq_list(input_tumor_name, cycle_edge_flux_list)
  # shift_node_freq_list = get_shift_node_freq_list(input_tumor_name, cycle_flux_path_df_list)
  # shift_edge_node_freq_list = get_shift_edge_node_freq_list(input_tumor_name, cycle_flux_path_df_list)

  result_list = list()
  result_list[["cycle_edge"]] = cycle_edge_flux_list
  result_list[["shift_edge"]] = cycle_flux_path_df_list

  return(result_list)
}







