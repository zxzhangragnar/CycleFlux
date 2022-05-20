#' get_topology_info function
#'
#' This function for cycle get_topology_info
#'
#' @param get_topology_info input_net_file output_path, res_path
#' @keywords get_topology_info
#' @export
#' @examples
#' get_topology_info(input_net_file, output_path, res_path)
#'
#'
get_topology_info <- function(input_net_file, output_path, res_path) {

  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))

  if(!file.exists(file.path(output_path, "cycle_directed.RData"))) {
    find_net_cycle(input_net_file, output_path)
  }

  source(system.file("rscript/1_cycle_topology/1_funcs/0_rela_module_cyc.R", package = "CycleFlux"), local=ENV)
  ENV$rela_module_cyc_main(output_path, res_path)

  source(system.file("rscript/1_cycle_topology/1_funcs/1_cyc_stcid_topology.R", package = "CycleFlux"), local=ENV)
  ENV$cycle_stcid_topology_main(output_path, res_path)

  source(system.file("rscript/1_cycle_topology/1_funcs/2_cycle_plot_and_distance.R", package = "CycleFlux"), local=ENV)
  ENV$cycle_plot_and_distance_main(output_path, res_path)

  setwd(ENV$old_wd)

  rm(ENV)
}

# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# get_topology_info(input_net_file, output_path, res_path)





####################################################################################

#' get_struct_sort_info function
#'
#' This function for cycle get_struct_sort_info
#'
#' @param get_struct_sort_info input_net_file output_path res_path
#' hub_prm_1 = 5,
#' hub_prm_2 = 5,
#' hub_prm_3 = 3,
#' web_prm_1 = .18,
#' web_prm_2 = 3,
#' web_prm_3 = 7,
#' web_prm_4 = 2,
#' web_prm_5 = 10,
#' idv_prm_1 = .16,
#' idv_prm_2 = 1,
#' idv_prm_3 = 2,
#' idv_prm_4 = 20,
#' idv_prm_5 = .56,
#' idv_prm_6 = 20,
#'
#'
#' @keywords get_struct_sort_info
#' @export
#' @examples
#' get_struct_sort_info()
#'
#'

get_struct_sort_info <- function(
  input_net_file,
  output_path,
  res_path,
  hub_prm_1 = 5,
  hub_prm_2 = 5,
  hub_prm_3 = 3,
  web_prm_1 = .18,
  web_prm_2 = 3,
  web_prm_3 = 7,
  web_prm_4 = 2,
  web_prm_5 = 10,
  idv_prm_1 = .16,
  idv_prm_2 = 1,
  idv_prm_3 = 2,
  idv_prm_4 = 20,
  idv_prm_5 = .56,
  idv_prm_6 = 20
) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))

  if(!file.exists(file.path(res_path, "1_cycle_topology/result_topo"))) {
    get_topology_info(input_net_file, output_path, res_path)
  }

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/1_struct_hub.R", package = "CycleFlux"), local=ENV)
  ENV$struct_hub_main(output_path, res_path, hub_prm_1, hub_prm_2, hub_prm_3)

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/2_struct_web.R", package = "CycleFlux"), local=ENV)
  ENV$struct_web_main(output_path, res_path, web_prm_1, web_prm_2, web_prm_3, web_prm_4, web_prm_5)

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/3_struct_idv.R", package = "CycleFlux"), local=ENV)
  ENV$struct_idv_main(output_path, res_path, idv_prm_1, idv_prm_2, idv_prm_3, idv_prm_4, idv_prm_5, idv_prm_6)

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/4_struct_other.R", package = "CycleFlux"), local=ENV)
  ENV$struct_other_main(output_path, res_path)

  setwd(ENV$old_wd)

  rm(ENV)
}


# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# res_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/res_files"
# get_struct_sort_info(input_net_file, output_path, res_path)


















