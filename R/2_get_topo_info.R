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
#' hub_para1 = 5,
#' hub_para2 = 5,
#' hub_para3 = 3,
#' web_para1 = .18,
#' web_para2 = 3,
#' web_para3 = 7,
#' web_para4 = 2,
#' web_para5 = 10,
#' idv_para1 = .16,
#' idv_para2 = 1,
#' idv_para3 = 2,
#' idv_para4 = 20,
#' idv_para5 = .56,
#' idv_para6 = 20,
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
  hub_para1 = 5,
  hub_para2 = 5,
  hub_para3 = 3,
  web_para1 = .18,
  web_para2 = 3,
  web_para3 = 7,
  web_para4 = 2,
  web_para5 = 10,
  idv_para1 = .16,
  idv_para2 = 1,
  idv_para3 = 2,
  idv_para4 = 20,
  idv_para5 = .56,
  idv_para6 = 20
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
  ENV$struct_hub_main(output_path, res_path, hub_para1, hub_para2, hub_para3)

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/2_struct_web.R", package = "CycleFlux"), local=ENV)
  ENV$struct_web_main(output_path, res_path, web_para1, web_para2, web_para3, web_para4, web_para5)

  source(system.file("rscript/2_cycle_struct_sort/1_funcs/3_struct_idv.R", package = "CycleFlux"), local=ENV)
  ENV$struct_idv_main(output_path, res_path, idv_para1, idv_para2, idv_para3, idv_para4, idv_para5, idv_para6)

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


















