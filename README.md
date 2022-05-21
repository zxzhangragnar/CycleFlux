# cycle_flux
cycle demo for test

# download

library(devtools)

install_github("zxzhangragnar/CycleFlux")

# load

library(CycleFlux)


# parameters
input_net_file = 'your_path/hsa_net.RData'

output_path = "your_path/output_files"

res_path = "your_path/res_files"

graph_path = "your_path/graph_files"

input_tumor_name = c("COAD", "ESCA")

input_pathway_name = "all"

input_tumor_data = c("your_path/TCGA-COAD.RData", "your_path/TCGA-ESCA.RData")

input_normal_data = c("your_path/TCGA-COAD_N.RData", "your_path/TCGA-ESCA_N.RData")


prm_1=0.05

prm_2=1

prm_3=1

prm_4=2

prm_5=0.5


# 1.

find_net_cycle(input_net_file, output_path)


# 2.

get_topology_info(input_net_file, output_path, res_path)

get_struct_sort_info(input_net_file, output_path, res_path)


# 3.

get_basic_gene_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)

get_subnet_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_1, prm_2, prm_3, prm_4)

get_cycle_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_1, prm_2, prm_3, prm_4)



# 4.

get_enzyme_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)

get_cycle_gapup_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_5)



# 5.

get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)

get_graph_single_cycle(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)

get_graph_single_cycle_degree_1(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)

get_cycle_shift(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)

#get_graph_single_cycle_degree_2(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)

get_graph_subnet(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)

#get_graph_topology_struct(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)


# 6.
get_freq_stat(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
