

import networkx as nx
from numpy import append

#now_cycle_node = all_node[0]
# 得到所有节点之间的距离
def get_all_pairs_dist(G):
    all_path_dist = dict(nx.all_pairs_shortest_path_length(G))
    return all_path_dist















##################################################################################################
# for get_pathway_subnet.py

def write_complete_compound_df(G, output_file_name):

    import pyreadr as pyreadr
    import pandas as pd

    #所有节点之间的距离
    all_node_distance = get_all_pairs_dist(G)
    list_df = pd.DataFrame()  #得到的列表

    for key in all_node_distance:
        # print(key)
        # print(all_node_distance[key])

        # distance
        cpdname = key
        distance = all_node_distance[key]

        cyc_df = pd.DataFrame([[cpdname, distance]], columns=["cpdname", "distance"])

        list_df = pd.concat((list_df, cyc_df), ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    pyreadr.write_rdata(output_file_name, list_df, df_name=str("cycasnode_net_distance"))    












#######################################################################################################
#######################################################################################################
#######################################################################################################
#
#    将这三种结构都用 具有相应度的节点替代
#    如何实现：删去环，并用相应的节点进行替换？
#
#    算法： 将图视为完全无向图，不考虑入度出度
#    stp1:在图G中删掉这个环上所有节点
#    stp2:然后增加一个新节点newnode
#    stp3:然后在图G中增加几条边，这些边为原来这个"环的邻接节点" 到新节点newnode 的边




#找到环的入度点和出度点
def get_in_neighcyc_detail_arr(neighcyc_detail_arr, cycles_arr, all_path_dist):
    in_neighcyc_detail_arr = dict()
    for cycno in neighcyc_detail_arr:  #环a
        temp_nodes = cycles_arr[cycno]
        neigh_cycs = neighcyc_detail_arr[cycno]
        now_sign_arr = list()

        for ncycno in neigh_cycs:      #环[b,c,d,e,f]
            neigh_nodes = cycles_arr[ncycno]
            for tn in temp_nodes:
                for nn in neigh_nodes:
                    try:
                        length = all_path_dist[tn][nn]
                        if length == 1:
                            now_sign_arr.append(ncycno)
                            break
                    except:
                        pass
        
        in_neighcyc_detail_arr[cycno] = now_sign_arr

    return in_neighcyc_detail_arr


def get_out_neighcyc_detail_arr(neighcyc_detail_arr, cycles_arr, all_path_dist):
    out_neighcyc_detail_arr = dict()
    for cycno in neighcyc_detail_arr:  #环a
        temp_nodes = cycles_arr[cycno]
        neigh_cycs = neighcyc_detail_arr[cycno]
        now_sign_arr = list()

        for ncycno in neigh_cycs:      #环[b,c,d,e,f]
            neigh_nodes = cycles_arr[ncycno]
            for tn in temp_nodes:
                for nn in neigh_nodes:
                    try:
                        length = all_path_dist[tn][nn]
                        if length == 1:
                            now_sign_arr.append(ncycno)
                            break
                    except:
                        pass
        
        out_neighcyc_detail_arr[cycno] = now_sign_arr

    return out_neighcyc_detail_arr




#########################################################################################################
#test
def del_cyc_add_node_main(G, cycles_arr, indegree_arr, outdegree_arr, in_neighcyc_detail_arr, out_neighcyc_detail_arr):

    ##################################################################################
    #原始状态
    # print("old node number")
    # print(len(G.nodes))

    newnode_list = list() #存放删掉 cycles 后 替代它的新节点

    #1.新建 老节点到新节点的边
    for cycno in range(len(cycles_arr)):        
        new_node = 'cyc' + str(cycno)
        newnode_list.append(new_node)
        #入度
        for old_ind_node in indegree_arr[cycno]:
            #print(len(G.nodes))
            if old_ind_node in G:
                G.add_edges_from([(old_ind_node, new_node)]) #新建老邻居到新节点的边
        #出度
        for old_outd_node in outdegree_arr[cycno]:
            #print(len(G.nodes))
            if old_outd_node in G:
                G.add_edges_from([(new_node, old_outd_node)]) #新建老邻居到新节点的边      

    #2.新建 老环到新环节点的边
    for cycno in range(len(newnode_list)):        
        new_cycnode = newnode_list[cycno]
        for neighcyc in in_neighcyc_detail_arr[cycno]:
            old_in_neigh = 'cyc' + str(neighcyc)
            G.add_edges_from([(old_in_neigh, new_cycnode)]) #新建老邻居环到新环的边

        for neighcyc in out_neighcyc_detail_arr[cycno]:
            old_out_neigh = 'cyc' + str(neighcyc)
            G.add_edges_from([(new_cycnode, old_out_neigh)]) #新建老邻居环到新环的边
    

    #3.删去原来的环所包含的节点
    for cycno in range(len(cycles_arr)): 
        temp_cyc_node = cycles_arr[cycno]
        for node in temp_cyc_node:
            #print(len(G.nodes))
            if node in G:
                G.remove_node(node)

    # print("new node number")
    # print(len(G.nodes))


    return G,newnode_list
##################################################################################



###################################################################################
### test
def get_new_G(G, merge_cycle_arr):
    cycles_arr = merge_cycle_arr
    print(cycles_arr)

    #总函数
    #环的邻居节点
    import cycle_degree
    indegree_arr = cycle_degree.get_cycle_ngbnode(G, cycles_arr, 'indegree')
    outdegree_arr = cycle_degree.get_cycle_ngbnode(G, cycles_arr, 'outdegree')

    #环的邻居环
    import cycle_degree_cyc
    neighcyc_arr = cycle_degree_cyc.cycle_neighcyc_func(G, cycles_arr)
    neighcyc_detail_arr = cycle_degree_cyc.cycle_neighcyc_func_detail(G, cycles_arr)

    # print("neighcyc_detail_arr")
    # print(neighcyc_detail_arr)

    # 增加：
    # in_neighcycs和out_neighcycs
    all_path_dist = dict(nx.all_pairs_shortest_path_length(G))



    ###
    in_neighcyc_detail_arr = get_in_neighcyc_detail_arr(neighcyc_detail_arr, cycles_arr, all_path_dist)
    out_neighcyc_detail_arr = get_out_neighcyc_detail_arr(neighcyc_detail_arr, cycles_arr, all_path_dist)


    #############################################
    cyc_to_node_res = del_cyc_add_node_main(G, cycles_arr, indegree_arr, outdegree_arr, in_neighcyc_detail_arr, out_neighcyc_detail_arr)
    new_G = cyc_to_node_res[0]

    return new_G

# new_G = get_new_G(G, merge_cycle_arr)
# # test write to RData
# write_complete_compound_df(new_G, 'cycasnode_net_distance.RData')

# print("ok")






###

# import get_pathway_subnet
# undirected_info = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = undirected_info[0]

# G_node_num = nx.number_of_nodes(G)
# G_edge_num = nx.number_of_edges(G)

# print(G_node_num)
# print(G_edge_num)







