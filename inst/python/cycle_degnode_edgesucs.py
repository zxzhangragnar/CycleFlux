# 
# Task:
#     找到：  1.环内 （环上某node --> 环外某node ）
# 			  2.出环 （环上某node --> 环外某node ）
# 			  3.进环 （环外某node --> 环上某node ）
#
#     找到环上每个节点的一步的出度
#############################################################################


#给二维字典添加元素
def addtwodimdict(thedict, key_a, key_b, val):
  if key_a in thedict:
    thedict[key_a].update({key_b: val})
  else:
    thedict.update({key_a:{key_b: val}})


#找出某两个元素C1,C2组成的边都包含于哪些环中
def get_edge_incycles(cycles_arr, C1, C2):
    incycles = list()
    for cyno in range(len(cycles_arr)):
        if (C1 in cycles_arr[cyno]) and (C2 in cycles_arr[cyno]):
            incycles.append(cyno)
    return incycles


######################################################################



#############################################################################################

#情况2
def write_edgesucs_out(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    #out_succ_nodes = get_situation_2(G, cycles_arr)
    import cycle_degree
    out_succ_nodes = cycle_degree.get_cycle_ngbnode_out_detail(G, cycles_arr)


    for cyc in out_succ_nodes:
        for node in out_succ_nodes[cyc]:
            for od in out_succ_nodes[cyc][node]:
                cycid = cyc
                node = node 
                od = od
                #reaction
                node_od_edge = G.get_edge_data(node, od)
                node_od_edge_Rid = node_od_edge['Rid']                 

                #添加行
                cyc_df = pd.DataFrame([[cycid, node, node_od_edge_Rid, od]], columns=["cycid", "node", "rid1", "od"])
                list_df = pd.concat((list_df, cyc_df), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_out.RData", list_df, df_name=str("cycle_edgesucs_out"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_out.RData", list_df, df_name=str("cycle_edgesucs_out"))





#情况3
def write_edgesucs_in(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    #in_succ_nodes = get_situation_3(G, cycles_arr)
    import cycle_degree
    in_succ_nodes = cycle_degree.get_cycle_ngbnode_in_detail(G, cycles_arr)

    for cyc in in_succ_nodes:
        for node in in_succ_nodes[cyc]:
            for ind in in_succ_nodes[cyc][node]:
                cycid = cyc
                node = node 
                ind = ind
                #reaction
                node_ind_edge = G.get_edge_data(ind, node)
                node_ind_edge_Rid = node_ind_edge['Rid']                 

                #添加行
                cyc_df = pd.DataFrame([[cycid, node, node_ind_edge_Rid, ind]], columns=["cycid", "node", "rid1", "ind"])
                list_df = pd.concat((list_df, cyc_df), ignore_index=True)


    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_in.RData", list_df, df_name=str("cycle_edgesucs_in"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_edgesucs_in.RData", list_df, df_name=str("cycle_edgesucs_in"))



##########################
# test

# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]



# write_edgesucs_out(G, cycles_arr, "directed")
# write_edgesucs_in(G, cycles_arr, "directed")



