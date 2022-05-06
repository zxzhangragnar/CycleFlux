
####### 1. 验证目前找到的环中的基因gene_symbol是否正确 ########
## cycle_directed
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_allpathway_cycle_union_directed.RData")

#cycle_edge_expression
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/cycle_edge_expression.RData")

#dp_part_net
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData")


# verify 验证
sel_rows = dp_part_net[which((dp_part_net$C_in == "C00021") & (dp_part_net$C_out == "C00019")),]
sel_rows = dp_part_net[which((dp_part_net$C_in == "C00019") & (dp_part_net$C_out == "C00021")),]


edge_express_genes = cycle_edge_expression[3, "gene_symbol"]
edge_express_genes = unlist(strsplit(edge_express_genes,split = ";"))
length1 = length(edge_express_genes)
edge_express_genes = unique(edge_express_genes)

dp_part_net_genes = c()
for (i in 1:length(sel_rows[,1])) {
  tmp_gene = as.character(sel_rows[i, "Gene_symbol"])
  tmp_gene = unlist(strsplit(tmp_gene,split = ";"))
  dp_part_net_genes = append(dp_part_net_genes, tmp_gene)
}
dp_part_net_genes = unique(dp_part_net_genes)



sort(dp_part_net_genes) == sort(edge_express_genes)
#全为TRUE






####### 2. simulate环数据, 并在input_flux_inout/codes中验证verify ########

## 输入数据为 cycle_edge_flux_list 的格式
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/4_flux_edge/result_final/cycle_edge_flux_list.RData")

sim_cyc = cycle_edge_flux_list[["COAD"]]

sim_cyc = sim_cyc[which(sim_cyc$cycid == 88),]

#fc 低的gene
# APOA4
# OTOP2
# APOC3
# SLC10A2

sim_cyc[1, "gene_symbol"] = "APOA4"

#fc 高的gene
# RNA5SP141
# MAGEA3
# RNA5SP204
# RNA5SP486


sim_cyc[2, "gene_symbol"] = "RNA5SP141"
sim_cyc[3, "gene_symbol"] = "RNA5SP486"
sim_cyc[4, "gene_symbol"] = "MAGEA3"



cycle_edge_flux_list = list()
cycle_edge_flux_list[["COAD"]] = sim_cyc

##跑input_flux_inout/codes中的1_ 2_ 3_即可
