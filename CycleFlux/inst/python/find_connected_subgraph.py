import networkx as nx

def find_connected_subgraph_func(G):

    #看图的连通性   找最大联通子图
    #先将有向图转成无向图（因为connected_components只对无向图适用）
    undirected_G = G.to_undirected()
    #再找最大联通子图
    conn_subgraph_arr = nx.connected_components(undirected_G)
    #排序
    conn_subgraph_arr = sorted(conn_subgraph_arr,key=len,reverse=True)
    return conn_subgraph_arr

#取最大连通子图(有向图)
def get_max_connected_subgraph_func(G):
    conn_subgraph_arr = find_connected_subgraph_func(G)
    G1 = G
    for i in range(1,len(conn_subgraph_arr)):  #删去全部独立小子图的节点，只保留最大的那个子图
        G1.remove_nodes_from(conn_subgraph_arr[i])
    return G1  #最大连通子图

###############################################################################
# test1
# 取全部的连通子图
""" import read_n_build_graph
graph_res = read_n_build_graph.read_n_build_graph_func('partnet.RData','part_net_after_delet')
G = graph_res[0]

conn_subgraph_arr = find_connected_subgraph_func(G)
print(conn_subgraph_arr)
#打印各子图的元素及其节点个数
for c in conn_subgraph_arr:
    print(c)     
    print(len(c))  

#子图数量 共115个子图
print(len(conn_subgraph_arr))
 """
###############################################################################
# test2
#G1为最大连通子图(有向)
""" G1 = get_max_connected_subgraph_func(G)

conn_subgraph_arr1 = find_connected_subgraph_func(G1)
print(len(conn_subgraph_arr1)) 
print(conn_subgraph_arr1[0]) 

for c in conn_subgraph_arr1:
    print(c)     
    print(len(c))  
 """





