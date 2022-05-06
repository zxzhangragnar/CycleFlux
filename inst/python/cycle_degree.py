#环的度节点
"""
包含两部分 directed 和 undirected
本函数内的degree指的是环周围的 节点(compound)的数目
(与之对应的是 cycle_degree_cyc 中的degree_cyc指的是 环周围 环(cyc)的数目)

本函数的all in out degree 只适用于 有向环 directed 
本函数的undegree 只适用于 无向环 undirected

找的各个环的度degree
算法：对环上所有的节点进行遍历，找到每个节点与自己距离为1(邻居节点)的点有哪些，刨去和自己同属一个环中的点，就是这个环的邻接点，
环的邻接点的个数就是环的度


有向图邻居节点
neighbors
有向图前向节点
predecessors
有向图后继节点
successors


[('C16643', 0), ('C11134', 0), ('C11136', 0), ('C00029', 1), 
('C05503', 0), ('C11131', 0), ('C00191', 0), ('C00190', 1), ('C05504', 0),
('C11376', 1), ('C19605', 0), ('C11132', 0), ('C14869', 0), ('C16578', 0),
('C00167', 20), ('C05787', 1), ('C11061', 0), ('C16577', 0), ('C11133', 0),
('C19606', 0), ('C11135', 0)]

"""

#######################################################################################################
######################################## directed ###################################################

import networkx as nx

#######################################################################################################
# 有向图 入度

# 有向图：统计各个节点的indegree度的详细信息，包括度的具体节点，和每个节点分别有多少个度 indegree用predecessors()
def get_node_indegree(G, nodename):
    import networkx as nx
    #nodename = 'C00029'
    #direct == "in_degree"
    #按有向环处理，找每个环的 入度
    neighbors = list(G.predecessors(nodename))
    # ['C00167', 'C00029', 'C00052', 'C00103']
    neighbor_degree_dict = dict() 
    # [('C00167', 1), ('C00029', 1), ('C00052', 1), ('C00103', 1)]
    for ngb in neighbors:
        neighbor_degree_dict[ngb] = 1
    
    return neighbors,neighbor_degree_dict

def get_cycle_degree_in(G, cycles_arr):
    cycle_degree_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        # 当前环的度
        temp_cycle_degree = 0
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_indegree(G, nodename)[0]
            neighbor_degree_dict = dict(get_node_indegree(G, nodename)[1])
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for interkey in same_cycle_nodes:
                neighbor_degree_dict.pop(interkey)
            #计算删除后的环的度
            for key in neighbor_degree_dict:
                temp_cycle_degree = temp_cycle_degree + neighbor_degree_dict[key]

        cycle_degree_dict[i] = temp_cycle_degree


    return cycle_degree_dict


def get_cycle_ngbnode_in(G, cycles_arr):
    cycle_ngbnode_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        # 当前环的度对应的具体节点
        temp_cycle_ngbnode_arr = list()
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_indegree(G, nodename)[0]
            neighbor_degree_dict = dict(get_node_indegree(G, nodename)[1])
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for interkey in same_cycle_nodes:
                neighbor_degree_dict.pop(interkey)
            #计算删除后的环的度 和 这些度对应的具体节点
            for key in neighbor_degree_dict:
                temp_cycle_ngbnode_arr.append(key)

        cycle_ngbnode_dict[i] = temp_cycle_ngbnode_arr

    return cycle_ngbnode_dict

#######################################################################################################
# 有向图 出度

# 有向图：统计各个节点的outdegree度的详细信息，包括度的具体节点，和每个节点分别有多少个度 outdegree用successors()
def get_node_outdegree(G, nodename):
    import networkx as nx
    #nodename = 'C00029'
    #direct == "out_degree"
    #按有向环处理，找每个环的 出度
    neighbors = list(G.successors(nodename))
    # ['C00167', 'C00029', 'C00052', 'C00103']
    neighbor_degree_dict = dict() 
    # [('C00167', 1), ('C00029', 1), ('C00052', 1), ('C00103', 1)]
    for ngb in neighbors:
        neighbor_degree_dict[ngb] = 1

    return neighbors,neighbor_degree_dict
 
def get_cycle_degree_out(G, cycles_arr):
    cycle_degree_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        # 当前环的度
        temp_cycle_degree = 0
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_outdegree(G, nodename)[0]
            neighbor_degree_dict = dict(get_node_outdegree(G, nodename)[1])
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for interkey in same_cycle_nodes:
                neighbor_degree_dict.pop(interkey)
            #计算删除后的环的度
            for key in neighbor_degree_dict:
                temp_cycle_degree = temp_cycle_degree + neighbor_degree_dict[key]

        cycle_degree_dict[i] = temp_cycle_degree


    return cycle_degree_dict

def get_cycle_ngbnode_out(G, cycles_arr):
    cycle_ngbnode_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        # 当前环的度对应的具体节点
        temp_cycle_ngbnode_arr = list()
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_outdegree(G, nodename)[0]
            neighbor_degree_dict = dict(get_node_outdegree(G, nodename)[1])
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for interkey in same_cycle_nodes:
                neighbor_degree_dict.pop(interkey)
            #计算删除后的环的度 和 这些度对应的具体节点
            for key in neighbor_degree_dict:
                temp_cycle_ngbnode_arr.append(key)

        cycle_ngbnode_dict[i] = temp_cycle_ngbnode_arr

    return cycle_ngbnode_dict


# for cycle_degnode_successors situation_1 situation_2
def get_cycle_ngbnode_out_detail(G, cycles_arr):
    cycle_ngbnode_detail_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        temp_cycle_ngbnode_detail_dict = dict()
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_outdegree(G, nodename)[0]
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for node in same_cycle_nodes:
                neighbor_nodes.remove(node)
            temp_cycle_ngbnode_detail_dict[nodename] = neighbor_nodes
            
        cycle_ngbnode_detail_dict[i] = temp_cycle_ngbnode_detail_dict
    return cycle_ngbnode_detail_dict



# for cycle_degnode_successors situation_3
def get_cycle_ngbnode_in_detail(G, cycles_arr):
    cycle_ngbnode_detail_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        temp_cycle_ngbnode_detail_dict = dict()
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_indegree(G, nodename)[0]
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for node in same_cycle_nodes:
                neighbor_nodes.remove(node)
            temp_cycle_ngbnode_detail_dict[nodename] = neighbor_nodes
            
        cycle_ngbnode_detail_dict[i] = temp_cycle_ngbnode_detail_dict
    return cycle_ngbnode_detail_dict



#######################################################################################################
# 有向图 所有度=入度+出度

# networkx_example/drawing/plot_degree.py
# 有向图：统计各个节点的degree度的详细信息，包括度的具体节点，和每个节点分别有多少个度
def get_node_alldegree(G, nodename):
    #direct == "all_degree" all_degree原名叫undirected
    #nodename = 'C00029'
    #找每个环的度
    import networkx as nx
    #neighbors = list(G.neighbors(nodename))
    neighbors = get_node_indegree(G, nodename)[0] + get_node_outdegree(G, nodename)[0] #所有度=入度+出度
    # ['C00167', 'C00029', 'C00052', 'C00103']
    neighbor_degree_dict = dict() 
    # [('C00167', 1), ('C00029', 1), ('C00052', 1), ('C00103', 1)]
    for ngb in neighbors:
        neighbor_degree_dict[ngb] = 1
    
    return neighbors,neighbor_degree_dict



#统计环的度
#算法：对环上所有的节点进行遍历，找到每个节点与自己距离为1的点有哪些，刨去和自己同属一个环中的点，就是这个环的邻接点，环的邻接点的个数就是环的度
def get_cycle_degree_all(G, cycles_arr):
    #所有度=入度+出度
    cycle_indegree_dict = get_cycle_degree_in(G, cycles_arr)
    cycle_outdegree_dict = get_cycle_degree_out(G, cycles_arr)
    cycle_degree_dict = dict()
    for i in range(len(cycles_arr)):
        cycle_degree_dict[i] = cycle_indegree_dict[i] + cycle_outdegree_dict[i]

    return cycle_degree_dict


# 统计环的度对应的具体节点  （环的邻接节点）
def get_cycle_ngbnode_all(G, cycles_arr):
    cycle_ngbnode_dict = dict()
    for i in range(len(cycles_arr)):
        temp_cycle = cycles_arr[i]
        # 当前环的度对应的具体节点
        temp_cycle_ngbnode_arr = list()
        for j in range(len(cycles_arr[i])):
            nodename = cycles_arr[i][j]
            neighbor_nodes = get_node_alldegree(G, nodename)[0]
            neighbor_degree_dict = dict(get_node_alldegree(G, nodename)[1])
            #取交集,找到与当前节点属于相同环的节点
            same_cycle_nodes = list(set(neighbor_nodes).intersection(set(temp_cycle)))   
            #删除字典中上述交集的键值对
            for interkey in same_cycle_nodes:
                neighbor_degree_dict.pop(interkey)
            #计算删除后的环的度 和 这些度对应的具体节点
            for key in neighbor_degree_dict:
                temp_cycle_ngbnode_arr.append(key)

        cycle_ngbnode_dict[i] = temp_cycle_ngbnode_arr

    return cycle_ngbnode_dict





##################
#总函数
def get_cycle_degree(G, cycles_arr, direct):
    if direct == 'indegree':      #入度
        return get_cycle_degree_in(G, cycles_arr)
    elif direct == 'outdegree':   #出度
        return get_cycle_degree_out(G, cycles_arr)
    else:                         #所有度 
        return get_cycle_degree_all(G, cycles_arr)

def get_cycle_ngbnode(G, cycles_arr, direct):
    if direct == 'indegree':      #入度
        return get_cycle_ngbnode_in(G, cycles_arr)
    elif direct == 'outdegree':   #出度
        return get_cycle_ngbnode_out(G, cycles_arr)
    else:                         #所有度 
        return get_cycle_ngbnode_all(G, cycles_arr)




#######################################################################################################
######################################## undirected ###################################################

#无向环的度 degree 此degree 指的是 环周围的节点
#函数写法来自 Functions_old cycle_significant.py 和 cycle_near_cyle.py

#此结果为 环 周围 有多少个节点
def get_undirected_cycle_degree(G, cycles_arr):
    
    cycle_significant_dict = dict()
    for i in range(len(cycles_arr)):
        now_sign_arr = list()
        now_cycle = i
        for j in range(len(cycles_arr[i])):
            compound = cycles_arr[i][j]
            temp_sign_arr = [n for n in G.neighbors(compound)]
            now_sign_arr= list(set(now_sign_arr).union(set(temp_sign_arr)))
        #排除掉当前这个环本身的化合物
        for comp in now_sign_arr:
            if comp in cycles_arr[i]:
                now_sign_arr.remove(comp)
        cycle_significant_dict[now_cycle] = len(now_sign_arr)
        
    return cycle_significant_dict

#此结果为 环 周围的 节点， 这些节点都是什么
def get_undirected_cycle_ngbnode(G, cycles_arr):
    
    cycle_significant_dict = dict()
    for i in range(len(cycles_arr)):
        now_sign_arr = list()
        now_cycle = i
        for j in range(len(cycles_arr[i])):
            compound = cycles_arr[i][j]
            temp_sign_arr = [n for n in G.neighbors(compound)]
            now_sign_arr= list(set(now_sign_arr).union(set(temp_sign_arr)))
        #排除掉当前这个环本身的化合物
        for comp in now_sign_arr:
            if comp in cycles_arr[i]:
                now_sign_arr.remove(comp)
        cycle_significant_dict[now_cycle] = now_sign_arr
        
    return cycle_significant_dict



#此结果为 节点 周围 有多少个节点
def get_undirected_cpd_degree(G, compound):
    
    neighbor_alldeg_nodes = [n for n in G.neighbors(compound)]

    return neighbor_alldeg_nodes
    












#######################################################
#test

# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

########################################################################################
# 主函数测试

################# directed #################
#所以环的度 dict   
#{0: 11, 1: 6, 2: 32, 3: 116, 4: 117, 5: 89, 6: 39, 7: 30, 8: 11, 9: 11, 10: 149}
"""
cycle_degree_dict = get_cycle_degree(G, cycles_arr, "alldegree")
print(cycle_degree_dict)
#所有环的度对应的具体节点  （环的邻接节点） dict
cycle_ngbnode_dict = get_cycle_ngbnode(G, cycles_arr, "alldegree")
print(cycle_ngbnode_dict)
"""
#201: ['C00243', 'C00611', 'C00984', 'C00197', 'C01231', 'C00394', 'C00668', 'C00167']}

#################
#入度 indegree
""" cycle_degree_dict = get_cycle_degree(G, cycles_arr, "indegree")
print(cycle_degree_dict)
cycle_ngbnode_dict = get_cycle_ngbnode(G, cycles_arr, "indegree")
print(cycle_ngbnode_dict)
#201: ['C00243', 'C00611', 'C00984', 'C00197', 'C01231', 'C00394', 'C00167']}
 """
#################
#出度 outdegree
""" cycle_degree_dict = get_cycle_degree(G, cycles_arr, "outdegree")
print(cycle_degree_dict)
cycle_ngbnode_dict = get_cycle_ngbnode(G, cycles_arr, "outdegree")
print(cycle_ngbnode_dict)
#201: ['C00243', 'C00611', 'C00984', 'C00668', 'C00197', 'C01231', 'C00167']}
 """


################# undirected #################
# undirected degree
"""
undirected_cycle_degree_dict = get_undirected_cycle_degree(G, cycles_arr)
print(undirected_cycle_degree_dict)
#所有环的度对应的具体节点  （环的邻接节点） dict
undirected_cycle_ngbnode_dict = get_undirected_cycle_ngbnode(G, cycles_arr)
print(undirected_cycle_ngbnode_dict)
"""



#######################################################################################
# 子函数测试
# nodename = 'C00167'
# direct = "directed"
# # neighbor_nodes = get_node_indegree_outdegree(G, nodename, direct)[0]
# # neighbor_degree_dict = get_node_indegree_outdegree(G, nodename, direct)[1]
# print("neighbor")
# neighbor_nodes = get_node_alldegree(G, nodename)[0]
# neighbor_degree_dict = get_node_alldegree(G, nodename)[1]
# nerighbor_nodes_len = len(neighbor_nodes)
# print(nerighbor_nodes_len)
# print(neighbor_nodes)
# print(neighbor_degree_dict)

# print("indegree")
# node_indegree = get_node_indegree(G, nodename)[0]
# in_neighbor_nodes = get_node_indegree(G, nodename)[1]
# print(node_indegree)
# print(in_neighbor_nodes)


# print("outdegree")
# node_outdegree = get_node_outdegree(G, nodename)[0]
# out_neighbor_nodes = get_node_outdegree(G, nodename)[1]
# print(node_outdegree)
# print(out_neighbor_nodes)




###################################################################################

# new test
# 找到某个节点compound的 deg 和 deg

### directed
# nodename = "C00118"
# neighbor_alldeg_nodes = get_node_alldegree(G, nodename)[0]
# neighbor_alldeg = len(neighbor_alldeg_nodes)

# print(neighbor_alldeg_nodes)
# print(neighbor_alldeg)

### undirected
# nodename = "C00118"
# neighbor_alldeg_nodes = get_undirected_cpd_degree(G, nodename)
# neighbor_alldeg = len(neighbor_alldeg_nodes)

# print(neighbor_alldeg_nodes)
# print(neighbor_alldeg)











#####################################################
## 2022.4.2
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

# #cycle_ngbnode_dict = get_cycle_ngbnode_out_detail(G, cycles_arr)
# cycle_ngbnode_dict = get_cycle_ngbnode_in_detail(G, cycles_arr)

# print("------------------------")
# print(cycle_ngbnode_dict[201])
# print("------------------------")
# print(cycle_ngbnode_dict[201]["C00319"])
# print("------------------------")
# print(cycle_ngbnode_dict[201]["C06124"])
# print("------------------------")


# print("------------------------")
# print(cycle_ngbnode_dict[374])
# print("------------------------")
# print(cycle_ngbnode_dict[374]["C00984"])
# print("------------------------")
# #print(cycle_ngbnode_dict[374]["C06124"])
# print("------------------------")


