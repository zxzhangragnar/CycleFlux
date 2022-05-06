###########################################################################

#1. 拆：每1行是一个化合物

# 看化合物和网的关系都是什么？
# (第一行为化合物1，第二行为化合物2)
# 第1列：化合物ID
# 第2列：化合物属于哪个环
# 第3列：化合物的反应RID
# 第4列：参与者和非参与者（2种情况：1.除本身外的cout 2.除本身外的cin）
# 第5列：化合物属于的环属于第几类结构（共3类结构）

# unprts数据格式： 
# 二级字典：R1:C1,C2...   R2:C3,C4....
# {'C00124': ['R10619':['C00001','C00002'], 'R01095':['C00003','C00004']]}

# {"参与者1"：['反应1':['此反应的非参与者1','非参与者2'],'反应2':['非参与者3','非参与者4']], "参与者2":[] }

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

# 给2维字典添加元素的函数
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


#获取compound的Rid的key键
def get_cpd_Rid_key(G, compound):
    cpd_rid_dict = get_cpd_Rid(G, compound)
    cpd_rid_keys = list()
    for key in cpd_rid_dict:
        cpd_rid_keys.append(key)
    return cpd_rid_keys



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





#####################################################################################
#获取 有向图 节点的 alldeg 具体节点 和 度的数目
def get_cpd_alldeg_directed(G, compound):
    import cycle_degree
    nodename = compound
    neighbor_alldeg_nodes = cycle_degree.get_node_alldegree(G, nodename)[0]
    neighbor_alldeg = len(neighbor_alldeg_nodes)

    return neighbor_alldeg, neighbor_alldeg_nodes

#获取 无向图 节点的 alldeg 具体节点 和 度的数目
def get_cpd_alldeg_undirected(G, compound):
    import cycle_degree
    nodename = compound
    neighbor_alldeg_nodes = cycle_degree.get_undirected_cpd_degree(G, nodename)
    neighbor_alldeg = len(neighbor_alldeg_nodes)

    return neighbor_alldeg, neighbor_alldeg_nodes



##################################################################################################
# for get_pathway_subnet.py
def write_complete_compound_df(G, cycles_arr, direct, output_file_name):

    list_df = get_complete_compound_df(G, cycles_arr, direct)
    #写入文件.RData
    import pyreadr as pyreadr
    pyreadr.write_rdata(output_file_name, list_df, df_name=str("compound_"+direct))    






## 3.21 for cycle_shared_cyc.py
def get_complete_compound_df(G, cycles_arr, direct):
    import pyreadr as pyreadr
    import pandas as pd
    
    import cycle_degree
    #compound node distance
    import cycle_distance
    all_node_distance = cycle_distance.get_all_pairs_dist(G)

    list_df = pd.DataFrame()  #得到的列表
    for cycno in range(len(cycles_arr)):
        print("cycle:" + str(cycno))
        for cpd in cycles_arr[cycno]:
            cpdname = cpd
            cycle = cycno
            Rid = get_cpd_Rid_key(G, cpd)
            #cpd 的度信息 deg
            if direct == "directed":
                deg = get_cpd_alldeg_directed(G, cpd)[0]
                degnodes = get_cpd_alldeg_directed(G, cpd)[1]
            else:
                deg = get_cpd_alldeg_undirected(G, cpd)[0]
                degnodes = get_cpd_alldeg_undirected(G, cpd)[1]

            prts_unprts = get_cpd_Rid(G, cpd)
            # distance
            distance = all_node_distance[cpdname]

            #indeg nodes
            indeg_nodes = cycle_degree.get_node_indegree(G, cpdname)[0]
            indeg_num = len(indeg_nodes)
            #outdeg nodes
            outdeg_nodes = cycle_degree.get_node_outdegree(G, cpdname)[0]
            outdeg_num = len(outdeg_nodes)

            cyc_df = pd.DataFrame([[cpdname, cycle, Rid, deg, degnodes, prts_unprts, distance, indeg_nodes, indeg_num, outdeg_nodes, outdeg_num]], columns=["cpdname", "cycle", "Rid", "deg", "degnodes", "prts_unprts", "distance", "indeg_nodes", "indeg_num", "outdeg_nodes", "outdeg_num"])

            list_df = list_df.append(cyc_df, ignore_index=True)
    
    return list_df



#########################################################################################################
#test
#见 get_pathway_subnet.py


