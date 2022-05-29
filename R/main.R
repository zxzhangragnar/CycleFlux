

input_net_file = 'hsa_subnet.RData'
output_path = "/output_files"
res_path = "/res_files"
graph_path = "/cycle_flux/graph_files"
input_tumor_name = "COAD"
input_pathway_name = "all"

prm_1=0.05
prm_2=1
prm_3=1
prm_4=2
prm_5=0.5

###########################################################################
source("1_cycle_directed.R")
source("2_get_missing_gene.R")
gene_missing_list = get_gene_and_tofind_list_main(package_path)
source("3_flux_subnet_add_gap.R")
subnet_edge_flux_list = subnet_edge_flux_list_main(output_path, res_path, package_path, input_tumor_name, prm_1, prm_2, prm_3, prm_4)
source("4_flux_cycle_add_gap.R")
cycle_edge_flux_list = cycle_edge_flux_list_main(output_path, res_path, package_path, input_tumor_name, prm_1, prm_2, prm_3, prm_4)

source("5_flux_cycle_chain_list.R")
all_chain_list_cid = get_cycle_chain_list_main(output_path, res_path, package_path, input_tumor_name)
compounds_dict = get_compounds_dict_main(package_path)

source("6_degnode_gap_analysis.R")
gapup_cycle_chain_list = single_cycle_gap_analysis_main(output_path, res_path, package_path, input_tumor_name)
never_considered_comp_names = never_considered_comp_names_main(output_path, res_path, package_path, input_tumor_name)

source("8_degnode_add_gap.R")
cycle_edge_flux_list_in = cycle_edge_flux_list_in_main(output_path, res_path, package_path, input_tumor_name, prm_1, prm_2, prm_3, prm_4)
cycle_edge_flux_list_out = cycle_edge_flux_list_out_main(output_path, res_path, package_path, input_tumor_name, prm_1, prm_2, prm_3, prm_4)

source("9_my_ug_degnode_1.R")
select_out_upinfo_eachedge = my_ug_degnode_1_main(output_path, res_path, input_tumor_name)
select_in_upinfo_eachedge = my_ug_degnode_2_main(output_path, res_path, input_tumor_name)

#### result ####

source("10_my_ug_degnode_3.R")
my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)

source("11_cyc_classification.R")
cycle_upgap_class_list = cyc_classification_main(res_path, input_tumor_name)

source("12_cyc_shift.R")
cycle_shift_path_df_list = cyc_shift_main(res_path, input_tumor_name)

source("13_freq_stat.R")
gap_cycleid_freq_list = get_gap_cycleid_freq_list(input_tumor_name, cycle_edge_flux_list)
shift_node_freq_list = get_shift_node_freq_list(input_tumor_name, cycle_shift_path_df_list)
shift_edge_node_freq_list = get_shift_edge_node_freq_list(input_tumor_name, cycle_shift_path_df_list)

