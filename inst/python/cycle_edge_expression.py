
#得到每个环所包含的基因gene
def write_cycle_edge_expression(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    #将所有环中的化合物排序，只取一个方向的顺序
    import find_cycles_order
    if direct == "directed":
        cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)
    else:
        cycle_order_path_arr = find_cycles_order.order_undirected_cycle_cpds(G, cycles_arr)
    for i in range(len(cycle_order_path_arr)):
        circle = cycle_order_path_arr[i]
        #cyc = circle[0]
        for k in range(len(circle)):
            cyc = circle[k]
            SG = G.subgraph(cyc)
            if direct == "undirected":
                SG = SG.to_undirected()
            # col1
            cycid = i
            for j in range(1,len(cyc)):
                pre_cpd = str(cyc[j-1])
                now_cpd = str(cyc[j])
                # col2 col3 col4 col5
                rid = SG[pre_cpd][now_cpd]['Rid']
                ec_express = SG[pre_cpd][now_cpd]['EC']
                gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                #添加行
                cyc_df = pd.DataFrame([[cycid, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "c_in", "rid", "c_out", "enzyme", "gene_symbol", "pathway"])
                list_df = list_df.append(cyc_df, ignore_index=True)
    
    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_edge_expression.RData", list_df, df_name=str("cycle_edge_expression"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_edge_expression.RData", list_df, df_name=str("cycle_edge_expression"))
















## test 
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

# write_cycle_edge_expression(G, cycles_arr, "directed")


