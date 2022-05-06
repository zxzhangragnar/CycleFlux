###########################################################################
# (第一行为环1，第二行为环2)
# 第1列：环1中全部代谢物
# 第2列：参与者(AD),非参与者(BC)  （例如：环上的节点中，AD参与反应R1,FH参与反应R2..E什么都没参与）
# 第3列：RID
# 第4列：环ID
# 第5列：属于第几类结构ID

# unprts数据格式： 
# 二级字典：R1:C1,C2...   R2:C3,C4....
# {'C00124': ['R10619':['C00001','C00002'], 'R01095':['C00003','C00004']]}

# {"参与者1"：['反应1':['此反应的非参与者1','非参与者2'],'反应2':['非参与者3','非参与者4']], "参与者2":[] }


# 非参与者 4种情况：           
# （这四种情况其实可以归纳为一种情况，即 包含这个compound的 Rid反应，除了当前compound 其他化合物均为 非参与者）
# 1.compound指向的节点  ok
# 2.指向compound的节点  ok
# 3.与compound共同作为输入的化合物的节点 
# 4.与compound共同作为输出的化合物的节点
# 	（怎么找？ 先找到compound在当前rid输出能得到的几种化合物,
# 	然后分别查找他们在当前rid的输入化合物有哪些？除去compound，即为非参与者）

###########################################################################
""" 

import pyreadr as pyreadr
import pandas as pd

# 对dataframe增加一行
# df = pd.DataFrame([["a",1],["b",2]], columns=["A", "B"])
# df2 = pd.DataFrame([["c",1],["d",2]], columns=["A", "B"])
# df = df.append(df2, ignore_index=True)
# print(df)

# prepare a pandas dataframe
# df = pd.DataFrame([["a",1],["b",2]], columns=["A", "B"])
#    A    B
# 1  a    1
# 2  b    2

#读取RData数据，并建图
import read_n_build_graph
graph_res = read_n_build_graph.read_n_build_graph_func('partnet.RData','part_net_after_delet')
G = graph_res[0]
cin_arr = graph_res[1]
cout_arr = graph_res[2]

# 在网络中找环
print("find cycles")
import find_cycles
cycles_arr = find_cycles.find_cyc_func("undirected", G) #无向图环
#1.环中全部代谢物：[['C00376', 'C20332', 'C20331', 'C02094'], ['C00154', 'C00249', 'C02588']]
#4.环ID 就是 cycle_arr[环ID]
#print(cycles_arr)

#找节点对应的RID G.edges.data()
edges_list = G.edges.data()
#print(edges_list) 
#二维数组[('C12270', 'C00025', {'Rid': 'R10687'}), ('C12270', 'C00025', {'Rid': 'R10687'})]
two_edges = [('C12270', 'C00025', {'Rid': 'R10687'}), ('C12270', 'C00025', {'Rid': 'R10687'})]
print(two_edges[0][0]) 
#C12270
print(two_edges[0][2]['Rid'])
#R10687

"""
#################################################################################


def addtwodimdict(thedict, key_a, key_b, val):
  if key_a in thedict:
    thedict[key_a].update({key_b: val})
  else:
    thedict.update({key_a:{key_b: val}})
#addtwodimdict(all_path_dist, now_cycle_node, temp_node, 0)

# 通过边的Rid属性找边
def get_edge_by_attribute(G, attr_name, attr_value):
    #attr_name = 'Rid', attr_value = 'R07003'
    selected_edges = [(u,v) for u,v,e in G.edges(data=True) if e[attr_name] == attr_value]
    #[('C00051', 'C14793'), ('C14787', 'C14793')]
    return selected_edges

#获取compound的Rid
def get_cpd_Rid(G, compound): #compound = 'C12270'
    print(str(compound))
    node_edges = list(G.edges(compound)) #得到这个化合物参与的各个edge
    #[('C12270', 'C01042'), ('C12270', 'C00025')]
    edge_data_arr = dict()
    for edge in node_edges:
        edge_data = G.get_edge_data(edge[0], edge[1]) #得到这条edge的Rid
        edge_data_arr[edge_data['Rid']] = list() #此字典用于存放每个Rid的非参与者
        Rid_edge_nodes = get_edge_by_attribute(G, 'Rid', edge_data['Rid']) #得到参与这个RID的所有化合物
        for formula in Rid_edge_nodes: #2维数组
            for node in formula:
                if node != compound:   #unprt 为除了当前环上compound 的其他所有参与反应的化合物
                    edge_data_arr[edge_data['Rid']].append(node)    #加入 unprt
        edge_data_arr[edge_data['Rid']] = list(set(edge_data_arr[edge_data['Rid']]))
    return edge_data_arr
#获取'C14787'的Rid
# cpd_Rid = get_cpd_Rid(G, 'C14787')
# print(cpd_Rid)
# {'R07003': ['C00051', 'C14793'], 'R07004': ['C14792', 'C00051'], 'R07014': ['C06205']}





#获取第cyno个环的 prts 和 unprts
def get_cyc_prts_unprts(G, cyno, cycles_arr):
    #print("get_cyc_prts_unprts")
    cyc_cpds = cycles_arr[cyno]
    cyc_prts_unprts_dict = dict()
    for cpd in cyc_cpds:
        cpd_Rid_dict = get_cpd_Rid(G, cpd)
        for key in cpd_Rid_dict:
            addtwodimdict(cyc_prts_unprts_dict, cpd, key, cpd_Rid_dict[key])
    return cyc_prts_unprts_dict
#cyc_prts_unprts_dict = get_cyc_prts_unprts(1, cycles_arr)
#print(cyc_prts_unprts_dict)
#{'C00124': ['R10619':['C00001','C00002'], 'R01095':['C00003','C00004']]}


#获取第cyno个环的Rids
def get_cyc_Rids(G, cyno, cycles_arr):
    #print("get_cyc_Rids")
    cyc_prts = get_cyc_prts_unprts(G, cyno, cycles_arr)
    cyc_Rids = list()
    for key1 in cyc_prts:
        for key2 in cyc_prts[key1]:
            cyc_Rids.append(key2) #key2 = R10619,R01095
    cyc_Rids = list(set(cyc_Rids))
    return cyc_Rids
#cyc_Rids = get_cyc_Rids(1, cycles_arr)
#print(cyc_Rids)






########################################################################################
#for get_pathway_subnet

########################################################################################
# 完全表信息 complete for get_pathway_subnet

#directed
def write_complete_pathway_subnet_directed_cycle_df(G, cycles_arr, ord_cpds_arr, ridord_cpds_arr, direct, output_file_name):

    import pyreadr as pyreadr
    import pandas as pd

    #得到incyc rids (环内的反应 按顺序排列)
    import find_cycles_order
    incyc_rid_arr = find_cycles_order.get_directed_incyc_rids(G, cycles_arr)

    import cycle_shared_cyc
    comp_cycle_dict = cycle_shared_cyc.get_comp_cycle_dict_directed(G, cycles_arr)
    cyc_near_cyc_dict = cycle_shared_cyc.cycle_shared_cyc_main_detail(cycles_arr, comp_cycle_dict)
    cyc_sign_dict = cycle_shared_cyc.cycle_shared_cyc_main(cyc_near_cyc_dict)
    import cycle_degree
    deg_arr = cycle_degree.get_cycle_degree(G, cycles_arr, "alldegree")
    

    list_df = pd.DataFrame()  #得到的列表

    #degree degree_detail
    import cycle_degree
    deg_arr = deg_arr
    indeg_arr = cycle_degree.get_cycle_degree(G, cycles_arr, "indegree")
    outdeg_arr = cycle_degree.get_cycle_degree(G, cycles_arr, "outdegree")

    deg_det_arr = cycle_degree.get_cycle_ngbnode(G, cycles_arr, "alldegree")
    indeg_det_arr = cycle_degree.get_cycle_ngbnode(G, cycles_arr, "indegree")
    outdeg_det_arr = cycle_degree.get_cycle_ngbnode(G, cycles_arr, "outdegree")
    
    #significant significant_detail
    sig_arr = cyc_sign_dict
    sigdet_arr = cyc_near_cyc_dict

    import cycle_degree_cyc
    neighcyc_arr = cycle_degree_cyc.cycle_neighcyc_func(G, cycles_arr)
    neighcyc_detail_arr = cycle_degree_cyc.cycle_neighcyc_func_detail(G, cycles_arr)

    #node_num
    # import cycle_node_nums
    # num_arr = cycle_node_nums.get_cycle_node_nums(cycles_arr, len(cycles_arr))

    #distance
    import cycle_distance
    distance_arr = cycle_distance.get_all_cycle_dist(G, cycles_arr)

    for cyno in range(len(cycles_arr)):
        cycid = cyno
        cpds = cycles_arr[cyno]
        num = len(cycles_arr[cyno])

        #order compound
        ord_cpds = ord_cpds_arr[cyno]
        ord_cpds_str = ridord_cpds_arr[cyno]
        
        prts_unprts = get_cyc_prts_unprts(G, cyno, cycles_arr)
        Rids = get_cyc_Rids(G, cyno, cycles_arr)
        incyc_rids = incyc_rid_arr[cyno]

        #degree degree_detail
        neighcyc = neighcyc_arr[cyno]
        neighcyc_detail = neighcyc_detail_arr[cyno]

        deg = deg_arr[cyno]
        indeg = indeg_arr[cyno]
        outdeg = outdeg_arr[cyno]

        deg_det = deg_det_arr[cyno]
        indeg_det = indeg_det_arr[cyno]
        outdeg_det = outdeg_det_arr[cyno]

        #significant significant_detail
        sig = sig_arr[cyno]
        sigdet = sigdet_arr[cyno]
        
        #distance
        distance = distance_arr[cyno]

        cyc_df = pd.DataFrame([[cycid, cpds, num, ord_cpds, ord_cpds_str, prts_unprts, Rids, incyc_rids, neighcyc, neighcyc_detail, deg, deg_det, indeg, indeg_det, outdeg, outdeg_det, sig, sigdet, distance]], columns=["cycid", "cpds",  "len", "ord_cpds", "ord_cpds_str", "prts_unprts", "Rids", "incyc_rids", "neighcyc", "neighcycs", "deg", "degs", "indeg", "indegs", "outdeg", "outdegs", "sharedcyc", "sharedcycs", "distance"])
        list_df = list_df.append(cyc_df, ignore_index=True)
    #写入文件.RData
    import pyreadr as pyreadr
    pyreadr.write_rdata(output_file_name, list_df, df_name=str("cycle_"+direct))



#undirected
def write_complete_pathway_subnet_undirected_cycle_df(G, cycles_arr, ord_cpds_arr, direct, output_file_name):

    import pyreadr as pyreadr
    import pandas as pd

    #得到incyc rids
    import find_cycles_order
    incyc_rid_arr = find_cycles_order.get_undirected_incyc_rids(G, cycles_arr)

    import cycle_shared_cyc
    comp_cycle_dict = cycle_shared_cyc.get_comp_cycle_dict_undirected(G, cycles_arr)
    cyc_near_cyc_dict = cycle_shared_cyc.cycle_shared_cyc_main_detail(cycles_arr, comp_cycle_dict)
    cyc_sign_dict = cycle_shared_cyc.cycle_shared_cyc_main(cyc_near_cyc_dict)

    list_df = pd.DataFrame()  #得到的列表
    
    #significant significant_detail
    import cycle_degree_cyc
    neighcyc_arr = cycle_degree_cyc.cycle_neighcyc_func(G, cycles_arr)
    neighcyc_detail_arr = cycle_degree_cyc.cycle_neighcyc_func_detail(G, cycles_arr)
    import cycle_degree
    undegree_arr = cycle_degree.get_undirected_cycle_degree(G, cycles_arr)


    sig_arr = cyc_sign_dict
    sigdet_arr = cyc_near_cyc_dict
    
    #node_num
    # import cycle_node_nums
    # num_arr = cycle_node_nums.get_cycle_node_nums(cycles_arr, len(cycles_arr))

    #distance
    import cycle_distance
    distance_arr = cycle_distance.get_all_cycle_dist(G, cycles_arr)

    for cyno in range(len(cycles_arr)):
        cycid = cyno
        cpds = cycles_arr[cyno]
        num = len(cycles_arr[cyno])
        #order compound
        ord_cpds = ord_cpds_arr[cyno]
        
        prts_unprts = get_cyc_prts_unprts(G, cyno, cycles_arr)
        Rids = get_cyc_Rids(G, cyno, cycles_arr)
        incyc_rids = incyc_rid_arr[cyno]
        
        #significant significant_detail
        neighcyc = neighcyc_arr[cyno]
        neighcyc_detail = neighcyc_detail_arr[cyno]
        undegree = undegree_arr[cyno]

        sig = sig_arr[cyno]
        sigdet = sigdet_arr[cyno]
 
        #distance
        distance = distance_arr[cyno]

        cyc_df = pd.DataFrame([[cycid, cpds, num, ord_cpds, prts_unprts, Rids, incyc_rids, neighcyc, neighcyc_detail, undegree, sig, sigdet, distance]], columns=["cycid", "cpds",  "len", "ord_cpds", "prts_unprts", "Rids", "incyc_rids", "neighcyc", "neighcycs", "degree", "sharedcyc", "sharedcycs", "distance"])
        list_df = list_df.append(cyc_df, ignore_index=True)
    #写入文件.RData
    import pyreadr as pyreadr
    pyreadr.write_rdata(output_file_name, list_df, df_name=str("cycle_"+direct))


#########################################################################################################
#test
#见 get_pathway_subnet.py












