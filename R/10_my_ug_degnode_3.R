
#e1<--up--(c2<--up--c1)<--up--e1
###################################################################

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
      #######################################################

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


##############################################################################
my_ug_degnode_3_main <- function(res_path, graph_path, input_tumor_name) {

  ##select_in_upinfo_eachedge select_out_upinfo_eachedge
  tumors_array = c(input_tumor_name)

  ##2.each_edge
  res_sub_path = "single_graphs"
  dir.create(file.path(graph_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(graph_path, res_sub_path)

  for (i in 1:length(tumors_array)) {
    tumor_name = tumors_array[i]

    plot_up_cycle(tumor_name, res_file_path, select_in_upinfo_eachedge[[tumor_name]], select_out_upinfo_eachedge[[tumor_name]], gapup_cycle_chain_list, never_considered_comp_names)

  }

}






# test
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
# input_tumor_name = "COAD"
# my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)


















##################################################################################
## disuse
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

