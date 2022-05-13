####
# 共应该有6个
# cycle_successors_1_expression_1 cycle_successors_1_expression_2
# cycle_successors_2_expression_1 cycle_successors_2_expression_2
# cycle_successors_3_expression_1 cycle_successors_3_expression_2
# 将这6个文件在cyc_successors_code的4个代码文件中跑一遍
# 得到6个输出结果
# 即为 这些环的successor的输入输出流的up和gap信息 
####




################################################################################################


###
def write_situation_comps_1(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd

    # 输入数据
    import cycle_degnode_successors
    self_loop_succ_nodes = cycle_degnode_successors.get_situation_1(G, cycles_arr)

    # 输出数据
    list_df_1 = pd.DataFrame()  #得到的列表
    for cyc in self_loop_succ_nodes:
        for node in self_loop_succ_nodes[cyc]:
            for od in self_loop_succ_nodes[cyc][node]:
                successor_incyc = self_loop_succ_nodes[cyc][node][od]
                if successor_incyc != []:
                    for suc in successor_incyc:     
                        pre_cpd = node
                        now_cpd = od

                        SG = G.subgraph([pre_cpd, now_cpd])
                        if direct == "undirected":
                            SG = SG.to_undirected()
                        
                        # col2 col3 col4 col5
                        rid = SG[pre_cpd][now_cpd]['Rid']
                        ec_express = SG[pre_cpd][now_cpd]['EC']
                        gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                        pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                        #添加行
                        tmp_row_df_1 = pd.DataFrame([[cyc, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "node", "rid", "od", "enzyme", "gene_symbol", "pathway"])
                        list_df_1 = pd.concat((list_df_1, tmp_row_df_1), ignore_index=True)

    # 输出数据
    list_df_2 = pd.DataFrame()  #得到的列表
    for cyc in self_loop_succ_nodes:
        for node in self_loop_succ_nodes[cyc]:
            for od in self_loop_succ_nodes[cyc][node]:
                successor_incyc = self_loop_succ_nodes[cyc][node][od]
                if successor_incyc != []:
                    for suc in successor_incyc:   
                        pre_cpd = od
                        now_cpd = suc

                        SG = G.subgraph([pre_cpd, now_cpd])
                        if direct == "undirected":
                            SG = SG.to_undirected()
                        
                        # col2 col3 col4 col5
                        rid = SG[pre_cpd][now_cpd]['Rid']
                        ec_express = SG[pre_cpd][now_cpd]['EC']
                        gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                        pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                        #添加行
                        tmp_row_df_2 = pd.DataFrame([[cyc, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "od", "rid", "suc", "enzyme", "gene_symbol", "pathway"])
                        list_df_2 = pd.concat((list_df_2, tmp_row_df_2), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_1_expression_1.RData", list_df_1, df_name=str("cycle_successors_1_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_1_expression_2.RData", list_df_2, df_name=str("cycle_successors_1_expression_2"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_1_expression_1.RData", list_df_1, df_name=str("cycle_successors_1_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_1_expression_2.RData", list_df_2, df_name=str("cycle_successors_1_expression_2"))





def write_situation_comps_2(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd

    # 输入数据
    import cycle_degnode_successors
    out_succ_nodes = cycle_degnode_successors.get_situation_2(G, cycles_arr)

    # 输出数据
    list_df_1 = pd.DataFrame()  #得到的列表
    for cyc in out_succ_nodes:
        for node in out_succ_nodes[cyc]:
            for od in out_succ_nodes[cyc][node]:
                for suc in out_succ_nodes[cyc][node][od]:
                    print(cyc)
                    pre_cpd = node
                    now_cpd = od

                    SG = G.subgraph([pre_cpd, now_cpd])
                    if direct == "undirected":
                        SG = SG.to_undirected()
                    
                    # col2 col3 col4 col5
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    ec_express = SG[pre_cpd][now_cpd]['EC']
                    gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                    pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                    #添加行
                    tmp_row_df_1 = pd.DataFrame([[cyc, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "node", "rid", "od", "enzyme", "gene_symbol", "pathway"])
                    list_df_1 = pd.concat((list_df_1, tmp_row_df_1), ignore_index=True)

    # 输出数据
    list_df_2 = pd.DataFrame()  #得到的列表
    for cyc in out_succ_nodes:
        for node in out_succ_nodes[cyc]:
            for od in out_succ_nodes[cyc][node]:
                for suc in out_succ_nodes[cyc][node][od]:
                    pre_cpd = od
                    now_cpd = suc

                    SG = G.subgraph([pre_cpd, now_cpd])
                    if direct == "undirected":
                        SG = SG.to_undirected()
                    
                    # col2 col3 col4 col5
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    ec_express = SG[pre_cpd][now_cpd]['EC']
                    gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                    pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                    #添加行
                    tmp_row_df_2 = pd.DataFrame([[cyc, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "od", "rid", "suc", "enzyme", "gene_symbol", "pathway"])
                    list_df_2 = pd.concat((list_df_2, tmp_row_df_2), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_2_expression_1.RData", list_df_1, df_name=str("cycle_successors_2_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_2_expression_2.RData", list_df_2, df_name=str("cycle_successors_2_expression_2"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_2_expression_1.RData", list_df_1, df_name=str("cycle_successors_2_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_2_expression_2.RData", list_df_2, df_name=str("cycle_successors_2_expression_2"))










def write_situation_comps_3(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd

    # 输入数据
    import cycle_degnode_successors
    in_succ_nodes = cycle_degnode_successors.get_situation_3(G, cycles_arr)

    # 输出数据
    list_df_1 = pd.DataFrame()  #得到的列表
    for cyc in in_succ_nodes:
        for node in in_succ_nodes[cyc]:
            for ind in in_succ_nodes[cyc][node]:
                for presuc in in_succ_nodes[cyc][node][ind]:
                    print(cyc)
                    pre_cpd = ind
                    now_cpd = node

                    SG = G.subgraph([pre_cpd, now_cpd])
                    if direct == "undirected":
                        SG = SG.to_undirected()
                    
                    # col2 col3 col4 col5
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    ec_express = SG[pre_cpd][now_cpd]['EC']
                    gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                    pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                    #添加行
                    tmp_row_df_1 = pd.DataFrame([[cyc, now_cpd, rid, pre_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "node", "rid", "ind", "enzyme", "gene_symbol", "pathway"])
                    list_df_1 = pd.concat((list_df_1, tmp_row_df_1), ignore_index=True)

    # 输出数据
    list_df_2 = pd.DataFrame()  #得到的列表
    for cyc in in_succ_nodes:
        for node in in_succ_nodes[cyc]:
            for ind in in_succ_nodes[cyc][node]:
                for presuc in in_succ_nodes[cyc][node][ind]:
                    print(cyc)
                    pre_cpd = presuc
                    now_cpd = ind

                    SG = G.subgraph([pre_cpd, now_cpd])
                    if direct == "undirected":
                        SG = SG.to_undirected()
                    
                    # col2 col3 col4 col5
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    ec_express = SG[pre_cpd][now_cpd]['EC']
                    gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                    pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                    #添加行
                    tmp_row_df_2 = pd.DataFrame([[cyc, now_cpd, rid, pre_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "ind", "rid", "presuc", "enzyme", "gene_symbol", "pathway"])
                    list_df_2 = pd.concat((list_df_2, tmp_row_df_2), ignore_index=True)
    
    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_3_expression_1.RData", list_df_1, df_name=str("cycle_successors_3_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_3_expression_2.RData", list_df_2, df_name=str("cycle_successors_3_expression_2"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_3_expression_1.RData", list_df_1, df_name=str("cycle_successors_3_expression_1"))
        pyreadr.write_rdata(output_path + "/cycle_successors_3_expression_2.RData", list_df_2, df_name=str("cycle_successors_3_expression_2"))




#########################################################################
# for cycle_degnode_edgesucs.py 和 suc1_codes
# 找出那些 长度为1的出度和入度


def write_edgesucs_out_comps(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd

    # 输入数据
    # comps_arr = read_edgesucs_out_comps(G, cycles_arr)
    # node_arr = comps_arr[0]
    # od_arr = comps_arr[1]

    import cycle_degree
    out_succ_nodes = cycle_degree.get_cycle_ngbnode_out_detail(G, cycles_arr)

    # 输出数据
    list_df_1 = pd.DataFrame()  #得到的列表
    for cyc in out_succ_nodes:
        for node in out_succ_nodes[cyc]:
            for od in out_succ_nodes[cyc][node]:
                print(cyc)
                pre_cpd = node
                now_cpd = od

                SG = G.subgraph([pre_cpd, now_cpd])
                if direct == "undirected":
                    SG = SG.to_undirected()
                
                # col2 col3 col4 col5
                rid = SG[pre_cpd][now_cpd]['Rid']
                ec_express = SG[pre_cpd][now_cpd]['EC']
                gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                pathway_express = SG[pre_cpd][now_cpd]['Pathway']
                #添加行
                tmp_row_df_1 = pd.DataFrame([[cyc, pre_cpd, rid, now_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "node", "rid", "od", "enzyme", "gene_symbol", "pathway"])
                list_df_1 = pd.concat((list_df_1, tmp_row_df_1), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_expression_out.RData", list_df_1, df_name=str("cycle_edgesucs_expression_out"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_expression_out.RData", list_df_1, df_name=str("cycle_edgesucs_expression_out"))


def write_edgesucs_in_comps(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd

    # 输入数据
    # comps_arr = read_edgesucs_in_comps(G, cycles_arr)
    # node_arr = comps_arr[0]
    # ind_arr = comps_arr[1]

    import cycle_degree
    in_succ_nodes = cycle_degree.get_cycle_ngbnode_in_detail(G, cycles_arr)

    # 输出数据
    list_df_1 = pd.DataFrame()  #得到的列表
    for cyc in in_succ_nodes:
        for node in in_succ_nodes[cyc]:
            for ind in in_succ_nodes[cyc][node]:
                print(cyc)
                pre_cpd = ind
                now_cpd = node

                SG = G.subgraph([pre_cpd, now_cpd])
                if direct == "undirected":
                    SG = SG.to_undirected()
                
                # col2 col3 col4 col5
                rid = SG[pre_cpd][now_cpd]['Rid']
                ec_express = SG[pre_cpd][now_cpd]['EC']
                gene_express = SG[pre_cpd][now_cpd]['Gene_symbol']
                pathway_express = SG[pre_cpd][now_cpd]['Pathway']

                #添加行
                tmp_row_df_1 = pd.DataFrame([[cyc, now_cpd, rid, pre_cpd, ec_express, gene_express, pathway_express]], columns=["cycid", "node", "rid", "ind", "enzyme", "gene_symbol", "pathway"])
                list_df_1 = pd.concat((list_df_1, tmp_row_df_1), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_expression_in.RData", list_df_1, df_name=str("cycle_edgesucs_expression_in"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_expression_in.RData", list_df_1, df_name=str("cycle_edgesucs_expression_in"))


#########################################################################





##
#test
###
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

# ## 2步度 表达
# write_situation_comps_1(G, cycles_arr, "directed")
# write_situation_comps_2(G, cycles_arr, "directed")
# write_situation_comps_3(G, cycles_arr, "directed")


# ## 1步度 表达
# write_edgesucs_out_comps(G, cycles_arr, "directed")
# write_edgesucs_in_comps(G, cycles_arr, "directed")

