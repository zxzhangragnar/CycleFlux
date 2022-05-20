#' get_enzyme_info function
#'
#' This function for cycle get_enzyme_info
#'
#' @param get_enzyme_info input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_enzyme_info
#' @export
#' @examples
#' get_enzyme_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'



get_enzyme_info <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data) {
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

  source(system.file("rscript/5_flux_cycle/1_funcs/4_get_enzyme_info_0.R", package = "CycleFlux"), local=ENV)
  ENV$get_enzyme_info_0_main(res_path, input_tumor_name, input_tumor_data, input_normal_data)

  source(system.file("rscript/5_flux_cycle/1_funcs/4_get_enzyme_info_1.R", package = "CycleFlux"), local=ENV)
  ENV$get_enzyme_info_1_main(input_net_file, output_path, res_path, package_path, input_tumor_name)

  setwd(ENV$old_wd)
  rm(ENV)
}


#' get_cycle_gapup_info function
#'
#' This function for cycle get_cycle_gapup_info
#'
#' @param get_cycle_gapup_info input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' prm_1=0.5
#' @keywords get_cycle_gapup_info
#' @export
#' @examples
#' get_cycle_gapup_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_1=0.5)
#'
#' if DE_cof > prm_1
#' -> DE = "obvious
#' else DE = "normal
#'

get_cycle_gapup_info <- function(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_1=0.5) {
  ENV = new.env()
  attach(ENV)

  library(reticulate)
  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "5_flux_cycle/result_enzyme"))) {
    get_enzyme_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/5_flux_cycle/1_funcs/5_cyc_info_1.R", package = "CycleFlux"), local=ENV)
  ENV$cyc_info_1_main(output_path, res_path, package_path, input_tumor_name, prm_1)

  source(system.file("rscript/5_flux_cycle/1_funcs/5_cyc_info_2.R", package = "CycleFlux"), local=ENV)
  ENV$cyc_info_2_main(res_path, input_tumor_name)

  source(system.file("rscript/5_flux_cycle/1_funcs/6_cyc_info_gap.R", package = "CycleFlux"), local=ENV)
  ENV$cyc_info_gap_main(res_path, input_tumor_name)

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
# prm_1=0.5
# get_enzyme_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
# get_cycle_gapup_info(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data, prm_1)
#

