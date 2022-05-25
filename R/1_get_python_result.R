#' find_net_cycle function
#'
#' This function for cycle find_net_cycle
#'
#' @param find_net_cycle input_net_file, output_path
#' @keywords find_net_cycle
#' @export
#' @examples
#' find_net_cycle(input_net_file, output_path)
#'
#'
find_net_cycle <- function(input_net_file, output_path) {
  ENV = new.env()
  attach(ENV)

  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))

  source(system.file("rscript/0_cycle_expression/1_funcs/1_cycle_directed.R", package = "CycleFlux"), local=ENV)
  ENV$get_expression_main(input_net_file, output_path)

  setwd(ENV$old_wd)

  rm(ENV)
}





# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# find_net_cycle(input_net_file, output_path)



