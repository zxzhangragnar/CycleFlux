########################################################################################
# 7.23 joint meeting

# 画出某个图
# def draw_this_cycle(G, cycles_arr, no):
#     #cyc = cycles_arr[117] #117环长度为9
#     cyc = cycles_arr[no]
#     print(cyc)
#     SG = G.subgraph(cyc)
#     import draw_graph
#     draw_graph.draw_graph(SG)


#######################################################################################

# tool function
# 将环中节点按前后链表连接顺序输出


import networkx as nx

def get_path_directed(SG, head, cycle_len):
    SG = SG.to_undirected()
    import networkx as nx
    paths = []
    for target in SG.neighbors(head):
        all_paths = list(nx.all_simple_paths(SG, source=head, target=target))
        for l in all_paths:
            if len(l) == cycle_len:
                paths = paths + [l + [head]]

    #3.24
    #必须与环上节点的个数相同
    final_path_res = list()
    path_res = paths
    for tp in path_res:
        if len(tp) == (cycle_len+1):
            final_path_res.append(tp)
    #保证数组中相邻节点 在环中也是相邻的
    tra_path_res = judge_traverse_tool(SG, final_path_res)

    #环387的情况, 有向环保留4个，再筛选
    # if (len(tra_path_res) > 2):
    #     tra_path_res = get_2_reverse_path(tra_path_res)
    
    return tra_path_res


#视作无向图找全部path
def get_path_undirected(SG, head, cycle_len):
    SG = SG.to_undirected()
    import networkx as nx
    paths = []
    for target in SG.neighbors(head):
        all_paths = list(nx.all_simple_paths(SG, source=head, target=target))
        for l in all_paths:
            if len(l) == cycle_len:
                paths = paths + [l + [head]]

    #3.24
    #必须与环上节点的个数相同
    final_path_res = list()
    path_res = paths
    for tp in path_res:
        if len(tp) == (cycle_len+1):
            final_path_res.append(tp)
    #保证数组中相邻节点 在环中也是相邻的
    tra_path_res = judge_traverse_tool(SG, final_path_res)

    #无向环保证最多只有2个, 若有两个, 则这两个互为rev的 即环387的情况
    if (len(tra_path_res) > 2):
        tra_path_res = get_2_reverse_path(tra_path_res)
    
    return tra_path_res



def get_2_reverse_path(test_2d_list):
    result_list = []
    for list1 in test_2d_list:
        if len(result_list) == 2:
            break
        for list2 in test_2d_list:
            if list1 == list(reversed(list2)):
                result_list.append(list1)
                result_list.append(list2)
                break
    return result_list

# tool function
# 判断找到的两种顺序是否都构成有向通路
#[['C00062', 'C00086', 'C00179', 'C00062'], ['C00062', 'C00179', 'C00086', 'C00062']]
# 是否正反顺序都能走到底(首尾相接)

# 若输入的SG是无向图 无向环 1->2<-3->1 此时不会被排除掉, 但会排除掉387环 多于2个path的情况
# 若输入的SG是有向图 无向环 1->2<-3->1 此时会被排除掉
def judge_traverse_tool(SG, final_path_res):
    result_list = list()
    for new_cyc in final_path_res:
        if_ord = True
        for i in range(len(new_cyc)-1):
            ssG = SG.subgraph([new_cyc[i],new_cyc[i+1]])
            if not nx.has_path(ssG, new_cyc[i], new_cyc[i+1]):
                #print(i)
                if_ord = False
        if if_ord:
            result_list.append(new_cyc)
    return result_list







def jud_traverse(G, cycle_order_path_arr):

    import networkx as nx
    # 判断 是否正反顺序都能走到底(首尾相接)
    juded_cycle_order_path_arr = list()
    for i in range(len(cycle_order_path_arr)):
        tmp_cyc_arr = cycle_order_path_arr[i]
        SG = G.subgraph(tmp_cyc_arr[0])
        this_cyc_arr = judge_traverse_tool(SG, tmp_cyc_arr)
        juded_cycle_order_path_arr.append(this_cyc_arr)

    return juded_cycle_order_path_arr


# 有向图 排序主函数
def order_directed_cycle_cpds(G, cycles_arr):
    cycle_order_path_arr = list()
    # 先按无向图处理，找到节点前后连接的顺序，得到两种正反顺序的数组
    for i in range(len(cycles_arr)):
        cyc = cycles_arr[i]
        SG = G.subgraph(cyc)  # 拿出无向环的几个节点构建一个新的子图
        #SG = SG.to_undirected()
        temp_path = get_path_directed(SG, cyc[0], len(cyc))
        
        cycle_order_path_arr.append(temp_path)

    # 判断 是否正反顺序都能走到底(首尾相接)
    juded_cycle_order_path_arr = jud_traverse(G, cycle_order_path_arr)

    # for c in range(len(juded_cycle_order_path_arr)):
    #     print("cycle:" + str(c))
    #     print(juded_cycle_order_path_arr[c])
    return juded_cycle_order_path_arr


# 无向图 排序主函数
def order_undirected_cycle_cpds(G, cycles_arr):
    cycle_order_path_arr = list()
    # 先按无向图处理，找到节点前后连接的顺序，得到两种正反顺序的数组
    for i in range(len(cycles_arr)):
        cyc = cycles_arr[i]
        SG = G.subgraph(cyc)  # 拿出无向环的几个节点构建一个新的子图
        #SG = SG.to_undirected()
        temp_path = get_path_undirected(SG, cyc[0], len(cyc))
        cycle_order_path_arr.append(temp_path)
    return cycle_order_path_arr

# for get_pathway_subnet.py 只取一个方向的顺序，['C00062', 'C00086', 'C00179', 'C00062']


def pthsub_order_undirected_cycle_cpds(G, cycles_arr):
    cycle_order_path_arr = list()
    # 先按无向图处理，找到节点前后连接的顺序，得到两种正反顺序的数组
    for i in range(len(cycles_arr)):
        cyc = cycles_arr[i]
        SG = G.subgraph(cyc)  # 拿出无向环的几个节点构建一个新的子图
        #SG = SG.to_undirected()
        temp_path = get_path_undirected(SG, cyc[0], len(cyc))
        cycle_order_path_arr.append(temp_path[0])
    return cycle_order_path_arr


# 9.25 incyc rids
def get_undirected_incyc_rids(G, cycles_arr):
    cycle_order_path_arr = pthsub_order_undirected_cycle_cpds(G, cycles_arr)
    circle = cycle_order_path_arr
    incyc_arr = list()
    for cyc in circle:
        circle_incyc_rid = list()
        SG = G.subgraph(cyc)
        SG = SG.to_undirected() #转成无向图
        for j in range(len(cyc)):
            if j == 0:
                rid = ""
            else:
                pre_cpd = str(cyc[j-1])
                now_cpd = str(cyc[j])
                rid = SG[pre_cpd][now_cpd]['Rid']
                circle_incyc_rid.append(rid)

        incyc_arr.append(circle_incyc_rid)

    return incyc_arr


#######################################################################################

# 有向环中 cpds 按首尾相接顺序排列并用rid链接
def ridorder_directed_cycle_cpds(G, cycles_arr):
    cycle_order_path_arr = order_directed_cycle_cpds(G, cycles_arr)
    ridorder_path_arr = list()
    for circle in cycle_order_path_arr:
        circle_str_arr = list()
        for cyc in circle:
            SG = G.subgraph(cyc)
            temp_cyc_str = ""
            for j in range(len(cyc)):
                if j == 0:
                    #temp_cyc_str = temp_cyc_str + str(cyc[j]) + '->'
                    temp_cyc_str = ""
                elif j == len(cyc)-1:
                    pre_cpd = str(cyc[j-1])
                    now_cpd = str(cyc[j])
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    temp_cyc_str = temp_cyc_str + pre_cpd + \
                        '->' + str(rid) + '->' + now_cpd
                else:
                    pre_cpd = str(cyc[j-1])
                    now_cpd = str(cyc[j])
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    temp_cyc_str = temp_cyc_str + \
                        pre_cpd + '->' + str(rid) + '->'

            circle_str_arr.append(temp_cyc_str)

        ridorder_path_arr.append(circle_str_arr)

    return ridorder_path_arr


# 9.25 incyc rids
def get_directed_incyc_rids(G, cycles_arr):
    cycle_order_path_arr = order_directed_cycle_cpds(G, cycles_arr)
    incyc_arr = list()
    for circle in cycle_order_path_arr:
        circle_incyc_arr = list()
        for cyc in circle:
            SG = G.subgraph(cyc)
            circle_incyc_rid = list()
            for j in range(len(cyc)):
                if j == 0:
                    rid = ""
                else:
                    pre_cpd = str(cyc[j-1])
                    now_cpd = str(cyc[j])
                    rid = SG[pre_cpd][now_cpd]['Rid']
                    circle_incyc_rid.append(rid)

            circle_incyc_arr.append(circle_incyc_rid)

        incyc_arr.append(circle_incyc_arr)

    return incyc_arr




####################################################################################
# import read_n_build_graph
# graph_res = read_n_build_graph.read_n_build_graph_func('partnet.RData','part_net_after_delet')
# G = graph_res[0]
# cin_arr = graph_res[1]
# cout_arr = graph_res[2]

#######################################
# order

# import find_cycles
# cycles_arr = find_cycles.find_cyc_func("directed", G)


# order_directed_cycle_cpds_arr = order_directed_cycle_cpds(G, cycles_arr)
# print(order_directed_cycle_cpds_arr)

# ridorder_directed_cycle_cpds_arr = ridorder_directed_cycle_cpds(G, cycles_arr)
# print(ridorder_directed_cycle_cpds_arr)


#######################################
# undirected
# import find_cycles
# cycles_arr = find_cycles.find_cyc_func("undirected", G)

# order_undirected_cycle_cpds_arr = pthsub_order_undirected_cycle_cpds(G, cycles_arr)
# print(order_undirected_cycle_cpds_arr)


# 画图
# draw_this_cycle(G, cycles_arr,4)
# draw_this_cycle(G, cycles_arr,23)



############################################################################################################
#########
#9.25 test
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]
# incyc_rid = get_undirected_incyc_rids(G, cycles_arr)
# print(len(incyc_rid))
# incyc_rid = get_directed_incyc_rids(G, cycles_arr)
# print(len(incyc_rid))

############################################################################################################
#########
##debug 3.24
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

## 无向
# cycle_order_path_arr = order_undirected_cycle_cpds(G, cycles_arr)
# for i in range(len(cycle_order_path_arr)):
#     print(str(i))
#     print(cycles_arr[i])
#     print(cycle_order_path_arr[i])
#     # if len(cycle_order_path_arr[i]) == 0:
#     #     print(str(i))
#     #     print(cycles_arr[i])
#     #     print(cycle_order_path_arr[i])

## 有向
# cycle_order_path_arr = order_directed_cycle_cpds(G, cycles_arr)
# for i in range(len(cycle_order_path_arr)):
#     # print(str(i))
#     # print(cycles_arr[i])
#     # print(cycle_order_path_arr[i])
#     if len(cycle_order_path_arr[i]) == 0:
#         print(str(i))
#         print(cycles_arr[i])
#         print(cycle_order_path_arr[i])

## bug环:cyc387
# cyc = cycles_arr[387]
# #print(cyc)
# SG = G.subgraph(cyc)  # 拿出无向环的几个节点构建一个新的子图

## 画出来这个环
# import matplotlib.pyplot as plt
# nx.draw(SG, with_labels=True)
# plt.show()
