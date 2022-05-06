#算法一：
#1.如何解决每次找到的环不一样的问题
# 最终解决方案：
# 无向图找环时，取 最大连通子图G1，并在调用nx.cycle_basis(G1, root)时设置根节点root
# 做到上述2点，每次找到的513个环cycles_arr的内容完全相同,上述各问题均可解决

#2.如何解决 simple_cycles MemoryError问题
# 算法二：
# 1.先按照无向图找环
# 2.定理：有向图的环一定为无向图的环的子集(有向图环在无向图中一定也是环，无向图环在有向图中不一定为环)
# 3.对找到的522个有向环取子图逐一判断看是不是有向环，对子图进行判断，看是否能找到一个有向环，且这个有向环包含子图全部节点

# 对于无向图每个环所包含的节点nodelist，构造一个新的有向图 SG = G.subgraph(nodelist)
# 在这个新的有向图中找环
# 若能找到包含这个有向图全部节点的环，则代表这个无向环同时也为有向图的环
# 将其加入到有向环结果list中即可	


#############################################################################################


import networkx as nx

#找无向图环的函数 底层调用的是nx.cycle_basis
def find_undirected_cycle(G):
    # G1为最大连通子图(有向)
    import find_connected_subgraph
    G1 = find_connected_subgraph.get_max_connected_subgraph_func(G)
    G1 = G1.to_undirected()
    # root 为找环时的根节点
    try: #用try 防止有的pathway子图找不到root
        root = list(G1.nodes)[0]
        cycles_arr = list(nx.cycle_basis(G1, root))
    except:
        #print("no root")
        cycles_arr = list(nx.cycle_basis(G1))
    #root = list(G1.nodes)[0]
    #cycles_arr = list(nx.cycle_basis(G1, root))
    #cycles_arr = list(nx.cycle_basis(G1))

    # 过滤环的长度 (2,16)
    jud_cycles_arr = [x for x in cycles_arr if(len(x) > 2 and len(x) < 16)]
    cycles_arr = jud_cycles_arr

    return cycles_arr


#找有向图环的函数相当于对nx.simple_cycles的封装
def find_directed_cycle(G):
    #存放有向图环的数组
    directed_cycle_arr = list()
    cycles_arr = find_undirected_cycle(G) #无向图环
    for cyc in cycles_arr:
        SG = G.subgraph(cyc) #拿出无向环的几个节点构建一个新的子图
        dir_simple_cycles = list(nx.simple_cycles(SG))  #对子图找有向环
        for dir_cyc in dir_simple_cycles: #若找到的环中存在一个与原来无向环中包含的节点完全相同的环
            if sorted(dir_cyc) == sorted(cyc): #判断两个数组是否相同
                directed_cycle_arr.append(cyc)

    return directed_cycle_arr


#选择主函数
def find_cyc_func(dir, G): #原有4个参数(dir, G, cin_arr, cout_arr)
    # 有向图找环
    if dir == "directed":
        #import tarjan_find_cycle
        #cycles_arr = list(tarjan_find_cycle.tarjan_find_cycle(cin_arr,cout_arr))
        cycles_arr = find_directed_cycle(G)
    # 无向图找环
    else:
        cycles_arr = find_undirected_cycle(G)
    #cycles_arr[第几个环][环中的第几个节点]
    return cycles_arr


###########################################################################
#test

# import pyreadr as pyreadr
# import networkx as nx

# #读取RData数据，并建图
# import read_n_build_graph
# graph_res = read_n_build_graph.read_n_build_graph_func('partnet.RData','part_net_after_delet')
# G = graph_res[0]
# cin_arr = graph_res[1]
# cout_arr = graph_res[2]

###########################################
#1.如何解决每次找到的环不一样的问题
#正解:取最大连通子图G1找环,并在找环时设置root节点

# 在网络中找环
"""
print("find cycles")
cycles_arr = find_cyc_func("undirected", G) #无向图环
"""
# #print(cycles_arr[514:522])
# print(cycles_arr)
# print(len(cycles_arr)) 


###########################################
#2.如何解决 simple_cycles MemoryError问题
#通过调用cycles_arr = list(nx.simple_cycles(G1))
# 正解：算法二
# 1.先按照无向图找环
# 2.定理：有向图的环一定为无向图的环的子集(有向图环在无向图中一定也是环，无向图环在有向图中不一定为环)
# 3.对找到的522个有向环取子图逐一判断看是不是有向环，对
"""
cycles_arr = find_cyc_func("directed", G) #有向图环

print(cycles_arr)
print(len(cycles_arr))
"""
############################################################





#####################################################################
# test 9.26
# import get_pathway_subnet
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]