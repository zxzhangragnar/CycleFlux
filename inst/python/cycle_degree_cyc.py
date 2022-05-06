#环的邻居环(中间有边相连)
"""
本文件 包含的两个函数同时 适用于 directed 和 undirected 只要传入相应的cycles_arr即可
本函数内的degree_cyc指的是环周围的 环(cyc)的数目

本函数是在原来的cycle_significant基础上进行改进

"""

#此结果为 环 周围 有多少个环
def cycle_neighcyc_func(G, cycles_arr):
    cycle_neighcyc_dict = dict()
    for i in range(len(cycles_arr)):
        now_sign_arr = list()
        now_cycle = i
        for j in range(len(cycles_arr[i])):
            compound = cycles_arr[i][j]
            compound_arr = [n for n in G.neighbors(compound)]
            comp_cycles_arr = compound_belongto_cycles(compound_arr, cycles_arr)
            now_sign_arr= list(set(now_sign_arr).union(set(comp_cycles_arr)))
            #排除掉当前这个环本身
            if now_cycle in now_sign_arr:
                now_sign_arr.remove(now_cycle)
        cycle_neighcyc_dict[now_cycle] = len(now_sign_arr)
        
    return cycle_neighcyc_dict

#此结果为 环 周围的 环， 这些环都是什么
def cycle_neighcyc_func_detail(G, cycles_arr):
    cycle_neighcyc_dict = dict()
    for i in range(len(cycles_arr)):
        now_sign_arr = list()
        now_cycle = i
        for j in range(len(cycles_arr[i])):
            compound = cycles_arr[i][j]
            compound_arr = [n for n in G.neighbors(compound)]
            comp_cycles_arr = compound_belongto_cycles(compound_arr, cycles_arr)
            now_sign_arr= list(set(now_sign_arr).union(set(comp_cycles_arr)))
            #排除掉当前这个环本身
            if now_cycle in now_sign_arr:
                now_sign_arr.remove(now_cycle)
        cycle_neighcyc_dict[now_cycle] = now_sign_arr
        
    return cycle_neighcyc_dict

#求一个compound数组都属于哪些环
def compound_belongto_cycles(compound_arr, cycles_arr):
    comp_cycles_arr = list()
    for comp in compound_arr:
        for i in range(len(cycles_arr)):
            if comp in cycles_arr[i]:
                comp_cycles_arr.append(i)

    comp_cycles_arr = list(set(comp_cycles_arr)) # 用set去重（去重后顺序是乱的）
    #comp_cycles_arr.sort(key=comp_cycles_arr.index) #加上本句顺序就不会乱了
    return comp_cycles_arr

######################################################################################
#test

# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]


#此结果为 环 周围 有多少个环
# cycle_significant_dict = cycle_neighcyc_func(G, cycles_arr)
# cycle_significant_dict = sorted(cycle_significant_dict.items(), key=lambda item:item[1], reverse=True)
# print(cycle_significant_dict)

#此结果为 环 周围的 环， 这些环都是什么
# cycle_significant_dict_detail = cycle_neighcyc_func_detail(G, cycles_arr)
#print(cycle_significant_dict_detail) 

