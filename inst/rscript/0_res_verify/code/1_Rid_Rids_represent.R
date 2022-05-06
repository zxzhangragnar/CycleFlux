###########################################################
## debug
# 同1个cin cout可能有多个Reaction, 在建图时用最后面的哪个Reaction覆盖前面所有的
# 即用最后面的那个 Rid代表全部的Rids
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/dp_part_net.RData")

test_dp = dp_part_net
for (i in 1:length(test_dp[,1])) {
  tmp_cpd = c(test_dp[i,"C_in"], test_dp[i,"C_out"])
  #tmp_cpd = sort(tmp_cpd) #注意 G[a][b]['val'] != G[b][a]['val']
  tmp_cpd_str = paste(tmp_cpd, collapse = ";")
  test_dp[i,"cpd"] = tmp_cpd_str
}

#part_dp = test_dp[,c("cpd","Rid")]
part_dp = test_dp[,c("cpd","Rid","Irreversible")]

#
table_part_dp =  aggregate(part_dp$Rid, list(part_dp$cpd), paste, collapse = ";")
colnames(table_part_dp) = c("cincout", "Rids")

#可逆反应
for (i in 1:length(part_dp[,1])) {
  temp_cincout = part_dp[i, "cpd"]
  temp_Rid = part_dp[i, "Rid"]
  temp_irreversible = part_dp[i, "Irreversible"]
  
  if((temp_irreversible == TRUE) | (is.na(temp_irreversible))) {
    temp_cincout_irreversible = unlist(strsplit(temp_cincout, split = ";"))
    temp_cincout_irreversible = rev(temp_cincout_irreversible)
    temp_cincout_irreversible = paste(temp_cincout_irreversible, collapse = ";")
    
    old_Rids_str = table_part_dp[which(table_part_dp$cincout == temp_cincout_irreversible), "Rids"]
    table_part_dp[which(table_part_dp$cincout == temp_cincout_irreversible), "Rids"] = paste0(old_Rids_str, ";", temp_Rid)
  }
}


#去重
for (i in 1:length(table_part_dp[,1])) {
  temp_Rids_str = table_part_dp[i, "Rids"]
  temp_Rids_arr = unlist(strsplit(temp_Rids_str,split = ";"))
  new_Rids_arr = unique(temp_Rids_arr)
  new_Rids_str = paste(new_Rids_arr, collapse = ";")
  table_part_dp[i, "Rids"] = new_Rids_str
}

#找出作为代表的Rid
for (i in 1:length(part_dp[,1])) {
  temp_cincout = part_dp[i, "cpd"]
  temp_rid = part_dp[i, "Rid"]
  temp_irreversible = temp_irreversible = part_dp[i, "Irreversible"]
  
  table_part_dp[which(table_part_dp$cincout == temp_cincout), "Rid"] = temp_rid
  #可逆反应的2个作为代表的Rid要相同
  if((temp_irreversible == TRUE) | (is.na(temp_irreversible))) {
    temp_cincout_irreversible = unlist(strsplit(temp_cincout, split = ";"))
    temp_cincout_irreversible = rev(temp_cincout_irreversible)
    temp_cincout_irreversible = paste(temp_cincout_irreversible, collapse = ";")
    #即始终和另一个可逆反应保持一致
    table_part_dp[which(table_part_dp$cincout == temp_cincout_irreversible), "Rid"] = temp_rid
  }
  
  
}

table_part_dp = table_part_dp[,c("cincout", "Rid", "Rids")]

save(table_part_dp, file="E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/0_res_verify/table_part_dp.RData")


#有多少个是Rids多于1个的
for (i in 1:length(table_part_dp[,1])) {
  temp_Rids_str = table_part_dp[i, "Rids"]
  temp_Rids_arr = unlist(strsplit(temp_Rids_str,split = ";"))
  table_part_dp[i, "len_Rids"] = length(temp_Rids_arr)
}

table(table_part_dp$len_Rids)
