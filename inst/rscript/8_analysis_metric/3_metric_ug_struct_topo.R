# web       1
# hub       2
# idv       3  
# nrl       4

#cycle_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/5_flux_cycle/result_final/cycle_flux_list.RData")

#cycle_edge_flux_list
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")

#gapup_cycle_chain_list.RData
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/6_graph_single_cycle/result_analysis/gapup_cycle_chain_list.RData")



################################################################################

###
tumor_name = "ESCA"
ug_c = names(gapup_cycle_chain_list[[tumor_name]])

######################
## updown_info
tmp_updown_info = cycle_flux_list[[tumor_name]]
ug_updown_info = tmp_updown_info[which(tmp_updown_info$cycle_id %in% ug_c),]
table(ug_updown_info$struct_id)

#exp_ug_c = c(37,55,63,123,225,227,278,344,345,347,357,405,426,469,500) #COAD
exp_ug_c = c(281,347,348,349,378,410,411,423,424,426) #ESCA

exp_ug_updown_info = tmp_updown_info[which(tmp_updown_info$cycle_id %in% exp_ug_c),]
table(exp_ug_updown_info$struct_id)
