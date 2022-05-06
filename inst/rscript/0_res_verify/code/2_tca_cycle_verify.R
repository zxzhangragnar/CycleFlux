################################### debug ######################################################################
################################### &&& ########################################################################
################################################################################################################
#6.调用python程序，得到结果(directed)res_allpathway_cycle_union_directed.RData
#load("E:/OneDrive - jlu/scFEA_Universal_cloud/my_R/aimA/rdata_cycle_detect/res_allpathway_cycle_union_directed.RData")
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_cycle_union_directed.RData")


##############################################################
#判断是否有tca主环
#判断是否有tca主环
tca_compound_cycles <- cycle_directed
tca_compound_cycles <- subset(tca_compound_cycles, grepl("C00042", cpds))
tca_compound_cycles <- subset(tca_compound_cycles, grepl("C00122", cpds))
tca_compound_cycles <- subset(tca_compound_cycles, grepl("C00149", cpds))
tca_compound_cycles <- subset(tca_compound_cycles, grepl("C00036", cpds))
rm(tca_compound_cycles)
##############################################################
#10.(可选)调用list_create文件夹中merge_compound_list.R方法合并compound
