
#######################################################################################################
# 目标：为cycle_edge_flux_list
# 增加8列
# "other_up_rct_num"
# "other_normal_rct_num"
# "other_up_rct_ratio"
# "other_up_highest_val"
# "other_normal_lowest_val"
# "other_up_highest_rct"
# "other_up_every_val"
# "other_all_every_val"

#######################################################################################################







#
# input_flux_inout的codes内执行顺序
# 4_flux_gap_reaction_addtable.R
# 3_flux_gap_reacion_pre.R
# 4_flux_gap_reaction.R 
# 4_flux_gap_reaction_after.R

##################################### part1 flux gap updown ###########################################
#######################################################################################################



# 
# 若这个环中有多个gap的断裂的点，
# 则以一个gap为支点，
# 两种方式
# (1)分析其它 nogap的反应
# (2)分析其它 的全部反应
#
#


#######################################################################################  
#方式一：(1)分析其它 nogap的反应
get_new_temp_part_cyc_df_situation1 <- function(cycle_directed, temp_cycle_edge_flux_list) {
  for (i in 1:length(cycle_directed[,1])) {
    cycid = cycle_directed[i,"cycid"]
    #cycid = 22
    
    temp_part_cyc_df = temp_cycle_edge_flux_list[temp_cycle_edge_flux_list$cycid == cycid,]
    #set new columns
    temp_part_cyc_df[, "other_up_rct_num"] = NA
    temp_part_cyc_df[, "other_normal_rct_num"] = NA
    temp_part_cyc_df[, "other_up_rct_ratio"] = NA
    temp_part_cyc_df[, "other_up_highest_val"] = NA
    temp_part_cyc_df[, "other_normal_lowest_val"] = NA
    temp_part_cyc_df[, "other_up_highest_rct"] = NA
    temp_part_cyc_df[, "other_up_every_val"] = NA
    temp_part_cyc_df[, "other_all_every_val"] = NA
    
    # gap part
    temp_gap_part_cyc_df = temp_part_cyc_df[temp_part_cyc_df$ifgap == "gap",]
    if (NA %in% temp_part_cyc_df$foldchangesign) {
      print(cycid)
    }else if (length(temp_gap_part_cyc_df[,1]) == 0) {
      print(cycid)
    }else {
      for (j in 1:length(temp_part_cyc_df[,1])) {
        if (temp_part_cyc_df[j,"ifgap"] == "gap") {
          # other part
          temp_rid = temp_part_cyc_df[j,"rid"]
          #temp_other_part_cyc_df = temp_part_cyc_df[temp_part_cyc_df$rid != temp_rid,]
          temp_other_part_cyc_df = temp_part_cyc_df[temp_part_cyc_df$ifgap != "gap",]
          
          # 老定义
          # other_down_rct_df = temp_other_part_cyc_df[temp_other_part_cyc_df$updown == "down",]
          # other_up_rct_df = temp_other_part_cyc_df[temp_other_part_cyc_df$updown == "up",]
          # 新定义
          other_down_rct_df = temp_other_part_cyc_df[temp_other_part_cyc_df$ifup == "normal",]
          other_up_rct_df = temp_other_part_cyc_df[temp_other_part_cyc_df$ifup == "up",]
          
          if(length(other_down_rct_df[,1]) == 0) {
            other_normal_rct_num = 0
            other_normal_lowest_val = 0
          }else {
            other_normal_rct_num = length(other_down_rct_df[,1]) 
            other_normal_lowest_val = min(other_down_rct_df$foldchangesign)
          }
          
          if(length(other_up_rct_df[,1]) == 0) {
            other_up_rct_num = 0
            other_up_rct_ratio = 0
            other_up_highest_val = 0
            other_up_highest_rct = 0
            other_up_every_val = 0
          }else {
            other_up_rct_num = length(other_up_rct_df[,1])
            other_up_rct_ratio = other_up_rct_num/(other_up_rct_num + other_normal_rct_num)
            other_up_highest_val = max(other_up_rct_df$foldchangesign)
            other_up_highest_rct = other_up_rct_df[other_up_rct_df$foldchangesign == other_up_highest_val,"rid"]
            if (length(other_up_highest_rct) > 1) {
              other_up_highest_rct = other_up_highest_rct[1]
            }
            other_up_every_val = mean(other_up_rct_df$foldchangesign)  
          }
          
          
          other_all_every_val = mean(temp_other_part_cyc_df$foldchangesign)
          
          ## put in dataframe    
          temp_part_cyc_df[j, "other_up_rct_num"] = other_up_rct_num
          temp_part_cyc_df[j, "other_normal_rct_num"] = other_normal_rct_num
          temp_part_cyc_df[j, "other_up_rct_ratio"] = other_up_rct_ratio
          temp_part_cyc_df[j, "other_up_highest_val"] = other_up_highest_val
          temp_part_cyc_df[j, "other_normal_lowest_val"] = other_normal_lowest_val
          temp_part_cyc_df[j, "other_up_highest_rct"] = other_up_highest_rct
          temp_part_cyc_df[j, "other_up_every_val"] = other_up_every_val
          temp_part_cyc_df[j, "other_all_every_val"] = other_all_every_val        
        }else {
          print("pass")
        }
        
      }
      
    }
    
    
    if(cycid == 0) {
      new_temp_part_cyc_df = as.data.frame(temp_part_cyc_df)
    }else {
      new_temp_part_cyc_df = rbind(new_temp_part_cyc_df,temp_part_cyc_df)
    }
    
  }
  
  return(new_temp_part_cyc_df)
}





#####################################################################################
add_cycleother_metric_main <- function(output_path, res_path, input_tumor_name) {
  # init
  load(file.path(output_path, "res_allpathway_cycle_union_directed.RData"))
  load(file.path(res_path, "4_flux_edge/result_final/cycle_edge_flux_list.RData"))
  
  ##
  tumors_array = c(input_tumor_name)
  
  new_cycle_edge_flux_list = list()
  for (i in 1:length(tumors_array)) {
    temp_cycle_edge_flux_list = cycle_edge_flux_list[[tumors_array[i]]]
    new_temp_part_cyc_df_s1 = get_new_temp_part_cyc_df_situation1(cycle_directed, temp_cycle_edge_flux_list)
    new_cycle_edge_flux_list[[tumors_array[i]]] = new_temp_part_cyc_df_s1
  }
  cycle_edge_flux_list = new_cycle_edge_flux_list
  
  
  #save
  res_sub_path = "4_flux_edge/result_final"
  dir.create(file.path(res_path, res_sub_path), recursive = TRUE, showWarnings = FALSE)
  res_file_path = file.path(res_path, res_sub_path, "cycle_edge_flux_list.RData")
  
  #save(cycle_edge_flux_list, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")
  save(cycle_edge_flux_list, file=res_file_path)
}




# test
# output_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output"
# res_path = "E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/"
# 
# input_tumor_name = "COAD"
# 
# add_cycleother_metric_main(output_path, res_path, input_tumor_name)
























