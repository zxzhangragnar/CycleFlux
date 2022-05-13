
#得到每个环所包含的基因gene
def write_subnet_edge_expression(G, output_path, direct):
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    all_edges = G.edges()
    print(len(all_edges))

    SG = G
    if direct == "undirected":
        SG = G.to_undirected()
    # col1
    for e in all_edges:
        pre_cpd = e[0]
        now_cpd = e[1]
        # col2 col3 col4 col5
        rid = SG[pre_cpd][now_cpd]['Rid']
        ec_express = SG[pre_cpd][now_cpd]['EC']
        gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
        pathway_express = SG[pre_cpd][now_cpd]['Pathway']

        #添加行
        cyc_df = pd.DataFrame([[pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["c_in", "rid", "c_out", "enzyme", "gene_symbol", "pathway"])
        list_df = pd.concat((list_df, cyc_df), ignore_index=True)
        

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/subnet_edge_expression.RData", list_df, df_name=str("subnet_edge_expression"))
    else:
        pyreadr.write_rdata(output_path + "/subnet_edge_expression.RData", list_df, df_name=str("subnet_edge_expression"))
















## test 
