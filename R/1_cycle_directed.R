

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




###
# edge_betweenness()
# edge.connectivity()
# edge_density()
# edge_disjoint_paths()

##################################################################################################################

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



##################################################################################################################

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

      ##
      cycle_node = directed_cycle_node[[i]]
      shared_cycles = c()
      for (j in 1:length(cycle_node)) {
        for (k in 1:length(directed_cycle_node)) {
          if(cycle_node[j] %in% directed_cycle_node[[k]]) {
            shared_cycles = c(shared_cycles, k)
          }
        }
      }

      shared_cycles = shared_cycles[shared_cycles != i]
      shared_cycles = unique(shared_cycles)
      cycle_directed[i, "shared_cycle"] = length(shared_cycles)
      cycle_directed[i, "shared_cycles"] = paste(shared_cycles, collapse = ";")

    }


    for(i in 1:length(directed_cycle)) {
      print(paste0("cyc_", i))
      cycle_directed[i, "cycle_id"] = i
      cycle_directed[i, "compound_chain"] = paste(directed_cycle[[i]], collapse = ";")

      ##
      cycle_node = directed_cycle_node[[i]]
      neighbor_cycles = c()
      for (j in 1:length(cycle_node)) {
        neighbor_nodes = neighbors(g, cycle_node[j], mode = "all")$name
        for (m in 1:length(neighbor_nodes)) {
          for (k in 1:length(directed_cycle_node)) {
            if(neighbor_nodes[m] %in% directed_cycle_node[[k]]) {
              neighbor_cycles = c(neighbor_cycles, k)
            }
          }
        }
      }
      neighbor_cycles = neighbor_cycles[neighbor_cycles != i]
      neighbor_cycles = unique(neighbor_cycles)
      cycle_directed[i, "neighcyc"] = length(neighbor_cycles)
      cycle_directed[i, "neighcycs"] = paste(neighbor_cycles, collapse = ";")
    }
  }

  return(cycle_directed)
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
get_cycle_edge_expression <- function(cycle_directed, g) {
  print("cycle_edge")
  cycle_edge_expression = data.frame()

  for(i in 1:length(rownames(cycle_directed))) {
    cycle_id = as.integer(cycle_directed[i, "cycle_id"])
    compound_chain = as.character(cycle_directed[i, "compound_chain"])
    compound_chain = unlist(strsplit(compound_chain, split = ";"))

    for (j in 1:length(compound_chain)) {
      cycle_node = unlist(strsplit(compound_chain[j], split = "->"))

      for (k in 1:(length(cycle_node)-1)) {
        temp_row = data.frame()
        c_in = cycle_node[k]
        c_out = cycle_node[k+1]

        temp_row[1, "cycle_id"] = cycle_id
        temp_row[1, "c_in"] = c_in
        temp_row[1, "c_out"] = c_out
        temp_row[1, "rid"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$rid
        temp_row[1, "enzyme"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$enzyme
        temp_row[1, "pathway"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$pathway
        temp_row[1, "gene_symbol"] = E(g)[get.edge.ids(g, c(c_in, c_out))]$gene_symbol

        cycle_edge_expression = rbind(cycle_edge_expression, temp_row)
      }
    }
  }

  return(cycle_edge_expression)
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
      node = cycle_node[k]
      ind_nodes = neighbors(g, node, mode = "in")$name
      ind_nodes = setdiff(ind_nodes, node)
      for (m in 1:length(ind_nodes)) {
        temp_row = data.frame()

        ind = ind_nodes[m]
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
      node = cycle_node[k]
      od_nodes = neighbors(g, node, mode = "out")$name
      od_nodes = setdiff(od_nodes, node)
      for (m in 1:length(od_nodes)) {
        temp_row = data.frame()

        od = od_nodes[m]
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


#####################################################################################
# subnet_edge_expression
# get_expression_main <- function(input_net_file, output_path) {
#   #init
#   load(file.path(input_net_file))
#
#   ##
#   g = build_net(hsa_net)
#   cycle_directed = get_cycle_directed(hsa_net, g)
#   subnet_edge_expression = get_subnet_edge_expression(g)
#   cycle_edge_expression = get_cycle_edge_expression(cycle_directed, g)
#   cycle_expression = get_cycle_expression(cycle_directed, g)
#   cycle_edgesucs_expression_in = get_cycle_edgesucs_expression_in(cycle_directed, g)
#   cycle_edgesucs_expression_out = get_cycle_edgesucs_expression_out(cycle_directed, g)
# }

#init
load(file.path(input_net_file))
##
g = build_net(hsa_net)
plot(g)
cycle_directed = get_cycle_directed(hsa_net, g)

subnet_edge_expression = data.frame()
cycle_edge_expression = data.frame()
cycle_expression = data.frame()
cycle_edgesucs_expression_in = data.frame()
cycle_edgesucs_expression_out = data.frame()

if(length(rownames(cycle_directed)) > 0) {
  subnet_edge_expression = get_subnet_edge_expression(g)
  cycle_edge_expression = get_cycle_edge_expression(cycle_directed, g)
  cycle_expression = get_cycle_expression(cycle_directed, g)
  cycle_edgesucs_expression_in = get_cycle_edgesucs_expression_in(cycle_directed, g)
  cycle_edgesucs_expression_out = get_cycle_edgesucs_expression_out(cycle_directed, g)
}


