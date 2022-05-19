#' get_basic_gene_info function
#'
#' This function for cycle get_basic_gene_info
#'
#' @param get_basic_gene_info input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_basic_gene_info
#' @export
#' @examples
#' get_basic_gene_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'


get_basic_gene_info <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(output_path, "subnet_edge_expression.RData"))) {
    find_net_cycle(input_net_file, output_path)
  }

  source(system.file("rscript/3_flux_subnet/1_funcs_tool/1_get_missing_gene.R", package = "CycleFlux"), local=ENV)
  ENV$get_gene_and_tofind_list_main(output_path, res_path, package_path)

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
# get_basic_gene_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)






####################################################################################

#' get_subnet_edge_info function
#'
#' This function for cycle get_subnet_edge_info
#'
#' @param get_subnet_edge_info input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_subnet_edge_info
#' @export
#' @examples
#' get_subnet_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'


get_subnet_edge_info <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "3_flux_subnet/result_tool"))) {
    get_basic_gene_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/3_flux_subnet/2_funcs/4_add_gap.R", package = "CycleFlux"), local=ENV)
  ENV$subnet_edge_flux_list_main(output_path, res_path, package_path, input_tumor_name)

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
# get_subnet_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)















####################################################################################

#' get_cycle_edge_info function
#'
#' This function for cycle get_cycle_edge_info
#'
#' @param get_cycle_edge_info input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' par_up1 = 1
#' par_up2 = 10
#' par_gap1 = -2
#' par_gap2 = 1.71
#' @keywords get_cycle_edge_info
#' @export
#' @examples
#' get_cycle_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, par_up1, par_up2, par_gap1, par_gap2)
#'
#' (max gene foldchange>(par_up1))&(max gene meanval>par_up2) -> up
#' (max gene foldchange<(par_gap1))|(max gene meanval<par_gap2) -> gap
#'
#'
#'
#'
#'
#'



get_cycle_edge_info <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "3_flux_subnet/result_tool"))) {
    get_basic_gene_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/4_flux_edge/1_funcs/4_add_gap.R", package = "CycleFlux"), local=ENV)
  ENV$cycle_edge_flux_list_main(output_path, res_path, package_path, input_tumor_name)

  source(system.file("rscript/4_flux_edge/1_funcs/5_get_cycle_chain_list.r", package = "CycleFlux"), local=ENV)
  ENV$get_cycle_chain_list_main(output_path, res_path, package_path, input_tumor_name)

  # source(system.file("rscript/4_flux_edge/1_funcs/6_add_cycleother_metric.R", package = "CycleFlux"), local=ENV)
  # ENV$add_cycleother_metric_main(output_path, res_path, input_tumor_name)

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
#
# par_up1 = 1
# par_up2 = 10
# par_gap1 = -2
# par_gap2 = 1.71
#
# get_cycle_edge_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, par_up1, par_up2, par_gap1, par_gap2)



