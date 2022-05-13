#环的公用环(互相之间有share公用节点)
######################################################################################
# 2021.7.30
# joint meeting 
# 310.
#
# 7.
# 更改 环的显著性significant 划分标准  
# 更改 环周围的环cycle near cycle 划分标准 
# 环周围的环cycle near cycle 和 significant 决定网状结构的划分

# 如何找到环的cycle near cycle：
# 相邻的环 互为 彼此的一部分 (拥有shared的节点)
# 环的significant表示，符合cyclenearcycle互为彼此一部分的环 有多少个

# 8.
# 环的degree：环和环之间以一条 edge(Rid) 相连 
# 环的cycle near cycle：相邻的环 互为 彼此的一部分
# 						          标准---至少有一个点为公用点
# 环的significant：环的significant表示，符合cyclenearcycle互为彼此一部分的环 有多少个
# 								  标准---符合cyclenearcycle标准的点有多少个 
# 环的 网状结构：一群 符合cyclenearcycle互为彼此一部分的环 连接在一起的结构
# （网上的每个环：环上的节点，都至少为其它2个环上的节点  则它们共同构成网状结构）
#


# 算法：
# 对每个环遍历,当前环里每个化合物都在哪些环中，这些环为这个环cyncy的环
# 每个节点，都在哪些环里


######################################################################################


def get_comp_cycle_dict_directed(G, cycles_arr):

    import pyreadr as pyreadr
    #result = pyreadr.read_r('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_compound_directed_merged.RData')
    
    #3.22
    #替代方案
    import write_compound_list
    complete_compound_df = write_compound_list.get_complete_compound_df(G, cycles_arr, "directed")

    import pandas as pd
    complete_compound_df = pd.DataFrame(complete_compound_df)
    complete_compound_df['cycle'] = complete_compound_df['cycle'].astype(str)
    compound_df_merged = complete_compound_df.groupby('cpdname').agg({'cycle':';'.join}).reset_index()

    cycle_list = list(compound_df_merged["cycle"].array)
    compound_list = list(compound_df_merged["cpdname"].array)
    comp_cycle_dict = dict()
    for i in range(len(cycle_list)):
        this_cyc_arr = cycle_list[i].split(';')
        for j in range(len(this_cyc_arr)):
            this_cyc_arr[j] = int(this_cyc_arr[j])

        this_comp = compound_list[i]
        comp_cycle_dict[this_comp] = this_cyc_arr

    return comp_cycle_dict


def get_comp_cycle_dict_undirected(G, cycles_arr):

    import pyreadr as pyreadr
    #result = pyreadr.read_r('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_compound_undirected_merged.RData')

    #3.22
    #替代方案
    import write_compound_list
    complete_compound_df = write_compound_list.get_complete_compound_df(G, cycles_arr, "undirected")

    import pandas as pd
    complete_compound_df = pd.DataFrame(complete_compound_df)
    complete_compound_df['cycle'] = complete_compound_df['cycle'].astype(str)
    compound_df_merged = complete_compound_df.groupby('cpdname').agg({'cycle':';'.join}).reset_index()

    cycle_list = list(compound_df_merged["cycle"].array)
    compound_list = list(compound_df_merged["cpdname"].array)

    comp_cycle_dict = dict()
    for i in range(len(cycle_list)):
        this_cyc_arr = cycle_list[i].split(';')
        for j in range(len(this_cyc_arr)):
            this_cyc_arr[j] = int(this_cyc_arr[j])

        this_comp = compound_list[i]
        comp_cycle_dict[this_comp] = this_cyc_arr

    return comp_cycle_dict


##########################################################################################






######################################################################
################# cycle_near_cycle 主函数 ################
# 计算环周围有多少个环,这些环具体都是哪些环 detail
def cycle_shared_cyc_main_detail(cycles_arr, comp_cycle_dict):
    # 所有节点共用的环 存放在字典dict中
    cyc_near_cyc_dict = dict()

    for i in range(len(cycles_arr)):
        shared_cycles = list() #存放当前环的共环
        for j in range(len(cycles_arr[i])):
            this_comp = cycles_arr[i][j]
            this_comp_shared_cycles = comp_cycle_dict[this_comp]
            #shared_cycles.append(this_comp_shared_cycles)
            shared_cycles= list(set(shared_cycles).union(set(this_comp_shared_cycles)))

        cyc_near_cyc_dict[i] = shared_cycles

    return cyc_near_cyc_dict



#此结果为 环 周围 有多少个环
def cycle_shared_cyc_main(cyc_near_cyc_dict):
    cycle_significant_dict = dict()
    for key in cyc_near_cyc_dict: 
        cycle_significant_dict[key] = len(cyc_near_cyc_dict[key])
    return cycle_significant_dict

######################################################################################
#test

# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]

# comp_cycle_dict = get_comp_cycle_dict_undirected()

# cyc_near_cyc_dict = cycle_shared_cyc_main_detail(cycles_arr, comp_cycle_dict)
# print(cyc_near_cyc_dict)

# cyc_sign_dict = cycle_shared_cyc_main(cyc_near_cyc_dict)
#print(cyc_sign_dict)





####
# new test for not use /main_output/res_compound_directed_merged.RData
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]


# import write_compound_list
# result = write_compound_list.get_complete_compound_df(G, cycles_arr, "directed")

# import pandas as pd

# import pyreadr as pyreadr
# #result = pyreadr.read_r('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_output/res_compound_directed_merged.RData')
# result = pd.DataFrame(result)
# result['cycle'] = result['cycle'].astype(str)
# agg_data = result.groupby('cpdname').agg({'cycle':';'.join}).reset_index()

# #new_agg_data = agg_data[['cpdname','cycle']]
# #print(new_agg_data)

# cycle_list = list(agg_data["cycle"].array)
# compound_list = list(agg_data["cpdname"].array)

# print(cycle_list)
# print(cycle_list[1])
# print(compound_list)
# print(compound_list[1])