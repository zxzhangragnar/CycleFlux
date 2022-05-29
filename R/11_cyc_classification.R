#
# each edge
#
# 1.Up only head in and tail out
# 2.There is a degree of up in the middle of Up
# 3.Another way
#
# 1.Up only head in and tail out
# Only 2 nodes in the ring have the degree of up
# And when these two nodes are used as c_in or c_out in this ring, ifup==up
#
# 2.There is a degree of up in the middle of Up
# There are >3 nodes in the ring with up degree
# And when these nodes are used as c_in or c_out in this ring, the ifgap at both ends are not all gaps
#
# 3.Another way
# There is a degree node (ind or od)
# This node is the degree of 2 nodes in this ring at the same time
# And when these two nodes are used as c_in or c_out in this ring, ifup==up
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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) == 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(!("up" %in% tmp_node_row[,"ifup"])){
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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) > 2) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(("gap" %in% tmp_node_row[,"ifgap"]) & !("up" %in% tmp_node_row[,"ifup"])){
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

    up_indeg = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),]
    up_outdeg = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),]

    up_indeg_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"ind"]
    up_outdeg_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"od"]

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
    up_indegnode_arr = tmp_cyc_indeg[which(tmp_cyc_indeg$ifup == "up"),"node"]
    up_outdegnode_arr = tmp_cyc_outdeg[which(tmp_cyc_outdeg$ifup == "up"),"node"]
    up_degnode_arr = append(up_indegnode_arr, up_outdegnode_arr)
    up_degnode_arr = unique(up_degnode_arr)

    if (length(up_degnode_arr) == 1) {
      for (j in 1:length(up_degnode_arr)) {
        tmpnode = up_degnode_arr[j]
        tmp_node_row = tmp_cyc[which((tmp_cyc$c_in == tmpnode) | (tmp_cyc$c_out == tmpnode)),]
        #tmp_node_row_ifup = tmp_node_row[,"ifup"]
        if(!("up" %in% tmp_node_row[,"ifup"])){
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



cyc_classification_main <- function(res_path, input_tumor_name) {

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




# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# input_tumor_name = "COAD"
# cyc_classification_main(res_path, input_tumor_name)




