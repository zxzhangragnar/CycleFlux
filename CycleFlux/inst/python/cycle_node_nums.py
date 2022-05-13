#每个环所包含的节点个数
def get_cycle_node_nums(cycles_arr, cycles_num):
    cycles_node_nums_arr = dict()   
    #cycles_arr[第几个环][环中的第几个节点]
    for i in range(cycles_num):
        cycles_node_nums_arr[i] = len(cycles_arr[i])
    #按环的长度排序
    #cycles_node_nums_arr = sorted(cycles_node_nums_arr.items(), key=lambda item:item[1])
    return cycles_node_nums_arr


#############################
#test
"""
import read_n_build_graph
graph_res = read_n_build_graph.read_n_build_graph_func('partnet.RData','part_net_after_delet')
G = graph_res[0]
cin_arr = graph_res[1]
cout_arr = graph_res[2]

import find_cycles
cycles_arr = find_cycles.find_cyc_func("undirected", G)
#cycles_arr = find_cycles.find_cyc_func("directed", G)
cycles_num = len(cycles_arr)

res = get_cycle_node_nums(cycles_arr, cycles_num)

print(res)
"""