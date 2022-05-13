# 
# Task:
#     找到：  1.环内 （环上某node --> 环外某node --> 环上某node）
# 			  2.出环 （环上某node --> 环外某node --> suc后继节点 ）
# 			  3.进环 （presuce前继节点 --> 环外某node --> 环上某node）
#
#
#     对某一个node 在表中要保证unique
#     若3中有这个节点，则2和1中就不要这个节点了(目前以下的代码取消了这一条件)
#     往后走一步的Reaction     #这些Reaction是否也在别的环里 1在 0不在
# 	  这些Reaction中的哪些也在别的环里，则日后这些不需要考虑
#############################################################################



#################



#给二维字典添加元素
def addtwodimdict(thedict, key_a, key_b, val):
  if key_a in thedict:
    thedict[key_a].update({key_b: val})
  else:
    thedict.update({key_a:{key_b: val}})


# #得到每个环所包含的反应Rid的字典
# def get_cycles_rids_dict(G, cycles_arr):
#     #将所有环中的化合物排序，只取一个方向的顺序
#     import find_cycles_order
#     cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)
#     cycles_rids_arr = list()
#     for circle in cycle_order_path_arr:
#         cyc = circle[0]
#         circle_str_arr = list()
#         SG = G.subgraph(cyc)
#         for j in range(1,len(cyc)):
#             pre_cpd = str(cyc[j-1])
#             now_cpd = str(cyc[j])
#             rid = SG[pre_cpd][now_cpd]['Rid']
#             circle_str_arr.append(rid)
#         cycles_rids_arr.append(circle_str_arr)

#     return cycles_rids_arr

# #得到每个环所包含的反应Rid的字典
# #cycles_rids_arr = get_cycles_rids_dict(G, cycles_arr)
# #print(cycles_rids_arr)



#找出某两个元素C1,C2组成的边都包含于哪些环中
def get_edge_incycles(cycles_arr, C1, C2):
    incycles = list()
    for cyno in range(len(cycles_arr)):
        if (C1 in cycles_arr[cyno]) and (C2 in cycles_arr[cyno]):
            incycles.append(cyno)
    return incycles


# print(get_edge_incycles(cycles_arr, "C01231", "C00236"))

######################################################################

#情况1 
def get_situation_1(G, cycles_arr):
    #出度 outdegree
    import cycle_degree
    cycle_ngbnode_dict = cycle_degree.get_cycle_ngbnode_out_detail(G, cycles_arr)

    ############## situation_1
    self_loop_succ_nodes = dict()
    #初始化
    for cyc in cycle_ngbnode_dict:
        for node in cycle_ngbnode_dict[cyc]:
            addtwodimdict(self_loop_succ_nodes, cyc, node, [])
    for cyc in cycle_ngbnode_dict:
        for node in cycle_ngbnode_dict[cyc]:
            od_arr = cycle_ngbnode_dict[cyc][node]
            temp_od_dict = dict()
            for od in od_arr:
                temp_od_list = list()
                od_successors = list(G.successors(od))
                for odsuc in od_successors:
                    if (odsuc in cycles_arr[cyc]) and (odsuc != node):
                        temp_od_list.append(odsuc)
                    #如环47也属于situation1的情况 https://www.kegg.jp/pathway/hsa00010
                    edge_incycles = get_edge_incycles(cycles_arr, od, odsuc)
                    if (odsuc == node) and (edge_incycles != []):
                    #if odsuc == node:
                        temp_od_list.append(odsuc)

                temp_od_dict[od] = temp_od_list

            addtwodimdict(self_loop_succ_nodes, cyc, node, temp_od_dict)
    
    #print(self_loop_succ_nodes)
    return self_loop_succ_nodes

###########################################################
#情况2

# 老方法, 已经废弃不用了
# def get_situation_2(G, cycles_arr):
#     #进环 （环上某node --> 环外某node ）

#     #出度 outdegree
#     import cycle_degree
#     cycle_ngbnode_dict = cycle_degree.get_cycle_ngbnode_out_detail(G, cycles_arr)

#     #选出环中不属于情况1的节点
#     out_succ_nodes = dict()
#     #初始化
#     for cyc in cycle_ngbnode_dict: 
#         for node in cycle_ngbnode_dict[cyc]:
#             od_arr = cycle_ngbnode_dict[cyc][node]
#             temp_od_dict = dict()
#             for od in od_arr:
#                 temp_od_list = list()
#                 od_successors = list(G.successors(od))
#                 cyc_judge = True
#                 for odsuc in od_successors:
#                     if odsuc in cycles_arr[cyc]:
#                         cyc_judge = False
#                     temp_od_list.append(odsuc)
#                 if cyc_judge:
#                     temp_od_dict[od] = temp_od_list
#                 else:
#                     temp_od_dict[od] = []
#             addtwodimdict(out_succ_nodes, cyc, node, temp_od_dict)
#     return out_succ_nodes
#     #看这些odsuc的Rid #这些Reaction是否也在别的环里 1在 0不在
#     #这些Reaction中的哪些也在别的环里，则日后这些不需要考虑



def get_situation_2(G, cycles_arr):
    #进环 （环上某node --> 环外某node ）

    #出度 outdegree
    import cycle_degree
    cycle_ngbnode_dict = cycle_degree.get_cycle_ngbnode_out_detail(G, cycles_arr)

    #选出环中不属于情况1的节点
    out_succ_nodes = dict()
    #初始化
    for cyc in cycle_ngbnode_dict: 
        for node in cycle_ngbnode_dict[cyc]:
            od_arr = cycle_ngbnode_dict[cyc][node]
            temp_od_dict = dict()
            for od in od_arr:
                temp_od_list = list()
                od_successors = list(G.successors(od))
                #cyc_judge = True
                for odsuc in od_successors:
                    if odsuc not in cycles_arr[cyc]:
                        #cyc_judge = False
                        temp_od_list.append(odsuc)
                temp_od_dict[od] = temp_od_list
            addtwodimdict(out_succ_nodes, cyc, node, temp_od_dict)
    return out_succ_nodes
    #看这些odsuc的Rid #这些Reaction是否也在别的环里 1在 0不在
    #这些Reaction中的哪些也在别的环里，则日后这些不需要考虑



###########################################################
#情况3

#老方法 不用了
# def get_situation_3(G, cycles_arr):
#     #进环 （环上某node --> 环外某node ）
#     #出度 outdegree
#     import cycle_degree
#     cycle_ngbnode_dict = cycle_degree.get_cycle_ngbnode_in_detail(G, cycles_arr)
#     #选出环中不属于情况1的节点
#     in_presucc_nodes = dict()
#     #初始化
#     for cyc in cycle_ngbnode_dict: 
#         for node in cycle_ngbnode_dict[cyc]:
#             ind_arr = cycle_ngbnode_dict[cyc][node]
#             temp_ind_dict = dict()
#             for ind in ind_arr:
#                 temp_od_list = list()
#                 ind_predecessors = list(G.predecessors(ind))
#                 cyc_judge = True
#                 for indsuc in ind_predecessors:
#                     if indsuc in cycles_arr[cyc]:   #情况4 与 情况1完全相同
#                         cyc_judge = False
#                     temp_od_list.append(indsuc)
#                 if cyc_judge:
#                     temp_ind_dict[ind] = temp_od_list
#                 else:
#                     temp_ind_dict[ind] = []
#             addtwodimdict(in_presucc_nodes, cyc, node, temp_ind_dict)
#     return in_presucc_nodes


def get_situation_3(G, cycles_arr):
    #进环 （环上某node --> 环外某node ）

    #出度 outdegree
    import cycle_degree
    cycle_ngbnode_dict = cycle_degree.get_cycle_ngbnode_in_detail(G, cycles_arr)

    #选出环中不属于情况1的节点
    in_presucc_nodes = dict()
    #初始化
    for cyc in cycle_ngbnode_dict: 
        for node in cycle_ngbnode_dict[cyc]:
            ind_arr = cycle_ngbnode_dict[cyc][node]
            temp_ind_dict = dict()
            for ind in ind_arr:
                temp_od_list = list()
                ind_predecessors = list(G.predecessors(ind))
                for indsuc in ind_predecessors:
                    if indsuc not in cycles_arr[cyc]:   #情况4 与 情况1完全相同
                        temp_od_list.append(indsuc)
                temp_ind_dict[ind] = temp_od_list
            addtwodimdict(in_presucc_nodes, cyc, node, temp_ind_dict)
    return in_presucc_nodes


###########################################################
#情况4
#从 环上的节点a  ---> 环上的节点b 这种情况  其实与 情况1 完全相同 所以不考虑

#############################################################################################

#情况1
def write_situation_1(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表
    self_loop_succ_nodes = get_situation_1(G, cycles_arr)

    for cyc in self_loop_succ_nodes:
        for node in self_loop_succ_nodes[cyc]:
            for od in self_loop_succ_nodes[cyc][node]:
                od = od
                cycid = cyc
                node = node
                successor_incyc = self_loop_succ_nodes[cyc][node][od]
                
                #reaction
                node_od_edge = G.get_edge_data(node, od)
                node_od_edge_Rid = node_od_edge['Rid']

                if successor_incyc != []:
                    for suc in successor_incyc:     
                        od_suc_edge = G.get_edge_data(od, suc)
                        od_suc_edge_Rid = od_suc_edge['Rid']
                        #针对环47的情况，判断为可逆但在别的环中的这种情况
                        od_reverse_judge = 0
                        reverse_edge_incycle = list()
                        if suc == node:
                            od_reverse_judge = 1
                            reverse_edge_incycle = get_edge_incycles(cycles_arr, od, suc)

                        cyc_df = pd.DataFrame([[cycid, node, node_od_edge_Rid, od, od_suc_edge_Rid, suc, od_reverse_judge, reverse_edge_incycle]], columns=["cycid", "node", "rid1", "od", "rid2", "suc_incyc", "is_reverse", "reverse_edge_incycle"])
                        list_df = pd.concat((list_df, cyc_df), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_1.RData", list_df, df_name=str("cycle_degnode_situation_1"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_1.RData", list_df, df_name=str("cycle_degnode_situation_1"))

#情况2
def write_situation_2(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    out_succ_nodes = get_situation_2(G, cycles_arr)

    for cyc in out_succ_nodes:
        for node in out_succ_nodes[cyc]:
            for od in out_succ_nodes[cyc][node]:
                for suc in out_succ_nodes[cyc][node][od]:
                    cycid = cyc
                    node = node 
                    od = od
                    suc = suc
                    #reaction
                    node_od_edge = G.get_edge_data(node, od)
                    node_od_edge_Rid = node_od_edge['Rid']
                    od_suc_edge = G.get_edge_data(od, suc)
                    od_suc_edge_Rid = od_suc_edge['Rid']                    

                    ifrid2_inothcyc_judge = True
                    rid2_inothcyc_no = get_edge_incycles(cycles_arr, od, suc)
                    if rid2_inothcyc_no != []:
                        ifrid2_inothcyc_judge = False

                    if ifrid2_inothcyc_judge:
                        ifrid2_inothcyc = 0
                    else:
                        ifrid2_inothcyc = 1
                    #添加行
                    cyc_df = pd.DataFrame([[cycid, node, node_od_edge_Rid, od, od_suc_edge_Rid, suc, ifrid2_inothcyc, rid2_inothcyc_no]], columns=["cycid", "node", "rid1", "od", "rid2", "suc", "ifrid2_inothcyc", "rid2_inothcyc_no"])
                    list_df = pd.concat((list_df, cyc_df), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_2.RData", list_df, df_name=str("cycle_degnode_situation_2"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_2.RData", list_df, df_name=str("cycle_degnode_situation_2"))





#情况3

def write_situation_3(G, cycles_arr, output_path, direct):
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表

    in_succ_nodes = get_situation_3(G, cycles_arr)

    for cyc in in_succ_nodes:
        for node in in_succ_nodes[cyc]:
            for ind in in_succ_nodes[cyc][node]:
                for suc in in_succ_nodes[cyc][node][ind]:
                    cycid = cyc
                    node = node 
                    ind = ind
                    suc = suc
                    #reaction
                    node_ind_edge = G.get_edge_data(ind, node)
                    node_ind_edge_Rid = node_ind_edge['Rid']
                    ind_suc_edge = G.get_edge_data(suc, ind)
                    ind_suc_edge_Rid = ind_suc_edge['Rid']                    

                    ifrid2_inothcyc_judge = True
                    rid2_inothcyc_no = get_edge_incycles(cycles_arr, suc, ind)
                    if rid2_inothcyc_no != []:
                        ifrid2_inothcyc_judge = False

                    if ifrid2_inothcyc_judge:
                        ifrid2_inothcyc = 0
                    else:
                        ifrid2_inothcyc = 1
                    #添加行
                    cyc_df = pd.DataFrame([[cycid, node, node_ind_edge_Rid, ind, ind_suc_edge_Rid, suc, ifrid2_inothcyc, rid2_inothcyc_no]], columns=["cycid", "node", "rid1", "ind", "rid2", "presuc", "ifrid2_inothcyc", "rid2_inothcyc_no"])
                    list_df = pd.concat((list_df, cyc_df), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_successors_3.RData", list_df, df_name=str("cycle_degnode_situation_3"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_successors_3.RData", list_df, df_name=str("cycle_degnode_situation_3"))



##########################
# test

# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]



# write_situation_1(G, cycles_arr, "directed")
# write_situation_2(G, cycles_arr, "directed")
# write_situation_3(G, cycles_arr, "directed")



