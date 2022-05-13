
get_compoundsVer <- function(compoundsVer) {
  compoundsVer120921<-as.matrix(compoundsVer)
  for(i in 1:nrow(compoundsVer120921))
  {
    compoundsVer120921[i,1]<-unlist(strsplit(compoundsVer120921[i,4],";",fixed=T))[1]
  }
  rownames(compoundsVer120921)<-compoundsVer120921[,1]
  compoundsVer120921 = as.data.frame(compoundsVer120921)

  return(compoundsVer120921)
}


####################################################################################

get_all_chain_list_cid<-function(tumor_name, cycle_edge_flux_list, cycle_directed) {
  #tumor_name = "COAD"
  cyc_metric = cycle_edge_flux_list[[tumor_name]]
  all_chain_cyc = list()
  for (i in 1:length(cycle_directed[,1])) {
    cycid = cycle_directed[i, "cycid"]
    tmp_cyc_metric = cyc_metric[which(cyc_metric$cycid == cycid),]
    tmp_cyc_chain = c()
    #1.生成反应链字符串
    chain_str = cycle_directed[i, "ord_cpds_str"]
    chain_str = substring(chain_str, 2,nchar(chain_str)-1)
    chain_str = unlist(strsplit(chain_str,split = ", "))
    for (j in 1:length(chain_str)) {
      tmp_chain_str = chain_str[j]
      tmp_chain_str = substring(tmp_chain_str, 2,nchar(tmp_chain_str)-1)

      tmp_chain_arr = unlist(strsplit(tmp_chain_str,split = "->"))
      #2.判断是否为up或gapp
      #开始处理
      for (k in 1:length(tmp_chain_arr)) {
        if(k%%2 == 0) {
          tmp_cin = tmp_chain_arr[k-1]
          tmp_cout = tmp_chain_arr[k+1]
          #print(tmp_chain_arr[k])
          # tmp_ifgap = tmp_cyc_metric[which(((tmp_cyc_metric$c_in == tmp_cin)|(tmp_cyc_metric$c_in == tmp_cout)) & ((tmp_cyc_metric$c_out == tmp_cin)|(tmp_cyc_metric$c_out == tmp_cout))), "ifgap"]
          # tmp_ifup = tmp_cyc_metric[which(((tmp_cyc_metric$c_in == tmp_cin)|(tmp_cyc_metric$c_in == tmp_cout)) & ((tmp_cyc_metric$c_out == tmp_cin)|(tmp_cyc_metric$c_out == tmp_cout))), "ifup"]
          tmp_ifgap = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifgap"]
          tmp_ifup = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifup"]

          if(tmp_ifgap == "gap") {
            tmp_chain_arr[k] = "G"
          }else if(tmp_ifup == "up") {
            tmp_chain_arr[k] = "U"
          }else {
            tmp_chain_arr[k] = ""
          }
        }
      }
      tmp_chain_result = paste(tmp_chain_arr, collapse = " -> ")
      #处理完毕
      tmp_cyc_chain = append(tmp_cyc_chain, tmp_chain_result)
    }
    all_chain_cyc[[i]] = tmp_cyc_chain
  }

  return(all_chain_cyc)
}

get_all_chain_list<-function(tumor_name, cycle_edge_flux_list, cycle_directed, compounds_dict) {
  #tumor_name = "COAD"
  cyc_metric = cycle_edge_flux_list[[tumor_name]]
  all_chain_cyc = list()
  for (i in 1:length(cycle_directed[,1])) {
    cycid = cycle_directed[i, "cycid"]
    tmp_cyc_metric = cyc_metric[which(cyc_metric$cycid == cycid),]
    tmp_cyc_chain = c()
    #1.生成反应链字符串
    chain_str = cycle_directed[i, "ord_cpds_str"]
    chain_str = substring(chain_str, 2,nchar(chain_str)-1)
    chain_str = unlist(strsplit(chain_str,split = ", "))
    for (j in 1:length(chain_str)) {
      tmp_chain_str = chain_str[j]
      tmp_chain_str = substring(tmp_chain_str, 2,nchar(tmp_chain_str)-1)

      tmp_chain_arr = unlist(strsplit(tmp_chain_str,split = "->"))
      #2.判断是否为up或gapp
      #开始处理
      #rid
      for (k in 1:length(tmp_chain_arr)) {
        if(k%%2 == 0) {
          tmp_cin = tmp_chain_arr[k-1]
          tmp_cout = tmp_chain_arr[k+1]
          #print(tmp_chain_arr[k])
          # tmp_ifgap = tmp_cyc_metric[which(((tmp_cyc_metric$c_in == tmp_cin)|(tmp_cyc_metric$c_in == tmp_cout)) & ((tmp_cyc_metric$c_out == tmp_cin)|(tmp_cyc_metric$c_out == tmp_cout))), "ifgap"]
          # tmp_ifup = tmp_cyc_metric[which(((tmp_cyc_metric$c_in == tmp_cin)|(tmp_cyc_metric$c_in == tmp_cout)) & ((tmp_cyc_metric$c_out == tmp_cin)|(tmp_cyc_metric$c_out == tmp_cout))), "ifup"]
          tmp_ifgap = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifgap"]
          tmp_ifup = tmp_cyc_metric[which((tmp_cyc_metric$c_in == tmp_cin) & (tmp_cyc_metric$c_out == tmp_cout)), "ifup"]

          if(tmp_ifgap == "gap") {
            tmp_chain_arr[k] = "G"
          }else if(tmp_ifup == "up") {
            tmp_chain_arr[k] = "U"
          }else {
            tmp_chain_arr[k] = ""
          }
        }
      }
      #cpd
      for (k in 1:length(tmp_chain_arr)) {
        if(k%%2 == 1) {
          tmp_chain_arr[k] = compounds_dict[which(compounds_dict$ENTRY == tmp_chain_arr[k]), "...1"]
        }
      }


      tmp_chain_result = paste(tmp_chain_arr, collapse = " -> ")
      #处理完毕
      tmp_cyc_chain = append(tmp_cyc_chain, tmp_chain_result)
    }
    all_chain_cyc[[i]] = tmp_cyc_chain
  }

  return(all_chain_cyc)
}


##############################################################################
get_cycle_chain_list_main <- function(output_path, res_path, package_path, input_tumor_name) {
  #init
  load(file.path(output_path, "cycle_directed.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))


  library(readr)
  compoundsVer = read_csv(file.path(package_path, "tool_data/compoundsVer120921.csv"))


  compounds_dict = get_compoundsVer(compoundsVer)
  tumors_array = c(input_tumor_name)

  ## cname
  all_chain_list = list()
  for (i in 1:length(tumors_array)) {
    temp_chain_list =  get_all_chain_list(tumors_array[i], cycle_edge_flux_list, cycle_directed, compounds_dict)
    all_chain_list[[tumors_array[i]]] = temp_chain_list
  }

  ## cid
  all_chain_list_cid = list()
  for (i in 1:length(tumors_array)) {
    temp_chain_list_cid =  get_all_chain_list_cid(tumors_array[i], cycle_edge_flux_list, cycle_directed)
    all_chain_list_cid[[tumors_array[i]]] = temp_chain_list_cid
  }


  #save
  res_sub_path = "4_flux_edge/result_final"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "TCGA_gap_all_chain_list.RData")
  res_file_path_cid = file.path(res_path, res_sub_path, "TCGA_gap_all_chain_list_cid.RData")
  res_file_path_compound_dict = file.path(res_path, res_sub_path, "compounds_dict.RData")

  # save(all_chain_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/TCGA_gap_all_chain_list.RData")
  # save(all_chain_list_cid, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/TCGA_gap_all_chain_list_cid.RData")
  save(all_chain_list, file=res_file_path)
  save(all_chain_list_cid, file=res_file_path_cid)
  save(compounds_dict, file=res_file_path_compound_dict)
}



# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# package_path = "E:/R/R-4.1.2/library/CycleFlux/rscript"
# input_tumor_name = "COAD"
#
# get_cycle_chain_list_main(output_path, res_path, package_path, input_tumor_name)






