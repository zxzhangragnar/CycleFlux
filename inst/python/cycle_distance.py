#计算环之间的距离， 每个环与其它所有环之间的距离
#遍历两个环上所有的节点，分别求他们之间的dij距离，取其中的最小值，即为两个环之间的距离


import networkx as nx
import pyreadr as pyreadr

##########################################################################
# 给2维字典添加元素的函数
def addtwodimdict(thedict, key_a, key_b, val):
  if key_a in thedict:
    thedict[key_a].update({key_b: val})
  else:
    thedict.update({key_a:{key_b: val}})


#now_cycle_node = all_node[0]
# 得到所有节点之间的距离
def get_all_pairs_dist(G):
    all_path_dist = dict(nx.all_pairs_shortest_path_length(G))
    return all_path_dist


#计算两个环之间的距离
def two_cycle_dist(this_cycle, that_cycle, all_path_dist):
    temp_dist = float("inf")   #inf正无穷，代表无路径
    for i in range(len(this_cycle)):
        this_node = this_cycle[i]
        for j in range(len(that_cycle)):
            that_node = that_cycle[j]
            #注意，若没有路径则会报错，要用try..catch
            try:
                length = all_path_dist[this_node][that_node]
            except:
                #print("no path")
                temp_dist = temp_dist
            else:
                if (length < temp_dist and length != 0):
                    temp_dist = length

    return temp_dist


#计算所有环之间的距离
def get_all_cycle_dist(G, cycles_arr):
    #转为无向图再求距离
    G = G.to_undirected()
    #所有节点之间路径的长度 
    all_path_dist = get_all_pairs_dist(G)
    #所有环之间路径的长度
    all_cycle_dist = dict()
    #cycles_arr = [[1,2,3],[4,5,6]]
    for i in range(len(cycles_arr)):
        for j in range(len(cycles_arr)):
            this_cycle = cycles_arr[i]
            that_cycle = cycles_arr[j]
            if i == j:
                addtwodimdict(all_cycle_dist, i, j, 0)  #距离为0，代表环到本身的距离
            else:
                dist = two_cycle_dist(this_cycle, that_cycle, all_path_dist)
                addtwodimdict(all_cycle_dist, i, j, dist) #两个环之间的距离 dist

    return all_cycle_dist





######################################################################################
#test
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

# res = get_all_cycle_dist(G, cycles_arr)
# print(res[3])
#

## 得到所有节点之间的距离
# all_node_distance = get_all_pairs_dist(G)
# print(all_node_distance['C01083'])

