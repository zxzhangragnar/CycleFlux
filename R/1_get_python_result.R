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
  library(reticulate)

  tryCatch({
    use_condaenv("r-reticulate")
  },error=function() {
    conda_create("r-reticulate")
    # install python packages
    conda_install("r-reticulate", "networkx")
    conda_install("r-reticulate", "pyreadr")
    conda_install("r-reticulate", "pandas")
    conda_install("r-reticulate", "numpy")
  })

  old_wd = getwd()
  ##new wd
  setwd(system.file("python/", package = "CycleFlux"))

  source_python(system.file("python/control_centra.py", package = "CycleFlux"))
  write_directed(input_net_file, output_path)

  setwd(old_wd)
}





# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# find_net_cycle(input_net_file, output_path)



