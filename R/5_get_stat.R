#' get_freq_stat function
#'
#' This function for cycle get_freq_stat
#'
#' @param get_freq_stat input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data
#' @keywords get_freq_stat
#' @export
#' @examples
#' get_freq_stat(input_net_file, output_path, res_path, input_tumor_name, input_tumor_data, input_normal_data)
#'

get_freq_stat <- function(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data) {
  ENV = new.env()
  attach(ENV)


  ENV$old_wd = getwd()

  ##new wd
  setwd(system.file("rscript/", package = "CycleFlux"))
  package_path = getwd() #"E:/R/R-4.1.2/library/CycleFlux/rscript"

  if(!file.exists(file.path(res_path, "6_graph_single_cycle/2_funcs_suc1/cyc_successors_data"))) {
    get_cycle_shift(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
  }

  source(system.file("rscript/8_analysis_metric/1_funcs/1_freq_stat.R", package = "CycleFlux"), local=ENV)
  ENV$freq_stat_main(res_path, input_tumor_name)

  setwd(ENV$old_wd)

  rm(ENV)
}


# get_freq_stat(input_net_file, output_path, res_path, graph_path, input_tumor_name, input_tumor_data, input_normal_data)
