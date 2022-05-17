#' get_graph_basic function
#'
#' This function for cycle get_graph_basic
#'
#' @param get_graph_basic input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_graph_basic
#' @export
#' @examples
#' get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_graph_basic <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "4_flux_edge/result_final"))) {
    get_cycle_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/6_graph_single_cycle/1_funcs/1_gap_analysis.R", package = "CycleFlux"), local=ENV)
  ENV$single_cycle_gap_analysis_main(output_path, res_path, package_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}


# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# input_tumor_name = "COAD"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)







####################################################################################

#' get_graph_single_cycle function
#'
#' This function for cycle get_graph_single_cycle
#'
#' @param get_graph_single_cycle input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_graph_single_cycle
#' @export
#' @examples
#' get_graph_single_cycle(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_graph_single_cycle <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/result_analysis"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/6_graph_single_cycle/1_funcs/2_cycle_draw.R", package = "CycleFlux"), local=ENV)
  #ENV$cycle_draw_cname_main(res_path, graph_path, input_tumor_name)
  ENV$cycle_draw_cid_main(res_path, graph_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}

# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_single_cycle(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)



####################################################################################

#' get_graph_single_cycle_degree_1 function
#'
#' This function for cycle get_graph_single_cycle_degree_1
#'
#' @param get_graph_single_cycle_degree_1 input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_graph_single_cycle_degree_1
#' @export
#' @examples
#' get_graph_single_cycle_degree_1(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_graph_single_cycle_degree_1 <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/result_analysis"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/1_funcs/4_add_gap.R", package = "CycleFlux"), local=ENV)
  ENV$cycle_edge_flux_list_inout_main(output_path, res_path, package_path, input_tumor_name)

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/2_funcs_draw/6_my_ug_degnode_1.R", package = "CycleFlux"), local=ENV)
  ENV$my_ug_degnode_1_main(output_path, res_path, input_tumor_name)
  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/2_funcs_draw/7_my_ug_degnode_2.R", package = "CycleFlux"), local=ENV)
  ENV$my_ug_degnode_2_main(output_path, res_path, input_tumor_name)
  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/2_funcs_draw/8_my_ug_degnode_3.R", package = "CycleFlux"), local=ENV)
  ENV$my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}


# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_single_cycle_degree_1(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)





####################################################################################

#' get_cycle_shift function
#'
#' This function for cycle get_cycle_shift
#'
#' @param get_cycle_shift input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_cycle_shift
#' @export
#' @examples
#' get_cycle_shift(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_cycle_shift <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/3_funcs_cyc_classification/1_cyc_classification.R", package = "CycleFlux"), local=ENV)
  ENV$cyc_classification_main(res_path, input_tumor_name)

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc1/3_funcs_cyc_classification/2_cyc_shift.R", package = "CycleFlux"), local=ENV)
  ENV$cyc_shift_main(res_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}



# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_cycle_shift(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)








####################################################################################

#' get_graph_single_cycle_degree_2 function
#'
#' This function for cycle get_graph_single_cycle_degree_2
#'
#' @param get_graph_single_cycle_degree_2 input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_graph_single_cycle_degree_2
#' @export
#' @examples
#' get_graph_single_cycle_degree_2(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_graph_single_cycle_degree_2 <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/result_analysis"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc2/1_funcs/4_add_gap.R", package = "CycleFlux"), local=ENV)
  ENV$cycle_successors_expression_main(output_path, res_path, package_path, input_tumor_name)

  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc2/1_funcs_draw/6_my_ug_degnode_1.R", package = "CycleFlux"), local=ENV)
  ENV$suc2_my_ug_degnode_1_main(output_path, res_path, input_tumor_name)
  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc2/1_funcs_draw/7_my_ug_degnode_2.R", package = "CycleFlux"), local=ENV)
  ENV$suc2_my_ug_degnode_2_main(output_path, res_path, input_tumor_name)
  source(system.file("rscript/6_graph_single_cycle/2_funcs_suc2/1_funcs_draw/8_my_ug_degnode_3.R", package = "CycleFlux"), local=ENV)
  ENV$suc2_my_ug_degnode_3_main(res_path, graph_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}


# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_single_cycle_degree_2(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)







####################################################################################

#' get_graph_subnet function
#'
#' This function for cycle get_graph_subnet
#'
#' @param get_graph_subnet input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data
#' @keywords get_graph_subnet
#' @export
#' @examples
#' get_graph_subnet(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)
#'


get_graph_subnet <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/result_analysis"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/7_graph_subnet/1_funcs/2_graph_subnet.R", package = "CycleFlux"), local=ENV)
  ENV$graph_subnet_main(output_path, res_path, graph_path, input_tumor_name, input_pathway_name)

  setwd(ENV$old_wd)

  rm(ENV)
}



# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_pathway_name = "hsa00350"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_subnet(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)







####################################################################################

#' get_graph_topology_struct function
#'
#' This function for cycle get_graph_topology_struct
#'
#' @param get_graph_topology_struct input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data
#' @keywords get_graph_topology_struct
#' @export
#' @examples
#' get_graph_topology_struct(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)
#'

get_graph_topology_struct <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/result_analysis"))) {
    get_graph_basic(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/7_graph_topo_struct/1_funcs/1_graph_struct.R", package = "CycleFlux"), local=ENV)
  ENV$graph_struct_main(output_path, res_path, graph_path, input_tumor_name, input_pathway_name)

  setwd(ENV$old_wd)

  rm(ENV)
}


# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# graph_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/graph_files"
#
# input_tumor_name = "COAD"
# input_pathway_name = "hsa00350"
# input_tumor_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD.RData"
# input_normal_data = "E:/scFEA_universal/Data/TCGA_data/TCGA_convolution/TCGA_data/TCGA-COAD_N.RData"
# get_graph_topology_struct(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_pathway_name, input_tumor_data, input_normal_data)

