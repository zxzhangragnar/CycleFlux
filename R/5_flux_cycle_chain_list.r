
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


####################################################################################

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
        tmp_ifgap = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifgap"]
        tmp_ifup = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifup"]

        if(tmp_ifgap == "gap") {
          tmp_str = c(tmp_str, "G")
        }else if(tmp_ifup == "up") {
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

# get_all_chain_list<-function(tumor_name, cycle_edge_flux_list, cycle_directed, compounds_dict) {
#   #tumor_name = "COAD"
#   cyc_metric = cycle_edge_flux_list[[tumor_name]]
#   all_chain_cyc = list()
#   for (i in 1:length(cycle_directed[,1])) {
#     cycle_id = cycle_directed[i, "cycle_id"]
#     tmp_cyc_metric = cyc_metric[which(cyc_metric$cycle_id == cycle_id),]
#     tmp_cyc_chain = c()
#
#     compound_chain = cycle_directed[i, "compound_chain"]
#     compound_chain = unlist(strsplit(compound_chain, split = ";"))
#     for (j in 1:length(compound_chain)) {
#       cycle_node = unlist(strsplit(compound_chain[j], split = "->"))
#       tmp_str = c(cycle_node[1])
#       for (k in 1:(length(cycle_node)-1)) {
#         tmp_cin = cycle_node[k]
#         tmp_cout = cycle_node[k+1]
#         tmp_ifgap = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifgap"]
#         tmp_ifup = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifup"]
#
#         if(tmp_ifgap == "gap") {
#           tmp_str = c(tmp_str, "G")
#         }else if(tmp_ifup == "up") {
#           tmp_str = c(tmp_str, "U")
#         }else {
#           tmp_str = c(tmp_str, "")
#         }
#         tmp_str = c(tmp_str, cycle_node[k+1])
#       }
#
#       for (k in 1:length(tmp_str)) {
#         if (tmp_str[k] %in% compounds_dict$ENTRY) {
#           tmp_str[k] = compounds_dict[which(compounds_dict$ENTRY == tmp_str[k]), "...1"]
#         }
#       }
#       tmp_chain_str = paste(tmp_str, collapse = " -> ")
#       tmp_cyc_chain = append(tmp_cyc_chain, tmp_chain_str)
#     }
#
#     all_chain_cyc[[i]] = tmp_cyc_chain
#   }
#
#   return(all_chain_cyc)
# }
#

##############################################################################
get_cycle_chain_list_main <- function(output_path, res_path, package_path, input_tumor_name) {
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



get_compounds_dict_main <- function(package_path) {
  #init
  library(readr)
  compoundsVer = read_csv(file.path("tool_data/compoundsVer120921.csv"))

  compounds_dict = get_compoundsVer(compoundsVer)

  return(compounds_dict)
}


# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
#
# get_cycle_chain_list_main(output_path, res_path, package_path, input_tumor_name)






