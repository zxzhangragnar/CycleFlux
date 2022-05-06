

#####################################################################################
#正解二：
# 1. 先整体找一次环
# 2. 再对每个pathway找一次环
# 3. 最后对每个pathway的comp子图找一次环
# 4. 3次的结果做并集
# 5. 得到全部环


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

################# 读数据 ###########################################################
#数据 from get_pathway_subnet.R

## 用于在python中生成 all_pathway_partnet_rowidx_df
def get_all_pathway_partnet_rowidx(dp_part_net):
  #1.pathway_arr
  pathway_arr = list(dp_part_net["Pathway"].array)
  pathway_arr = list(set(pathway_arr))

  if '' in pathway_arr:
    pathway_arr.remove('')

  #2.标记每个pathway对应dp_part_net行号
  all_pathway_partnet_rowidx = list()
  for p in pathway_arr:
    temp_pathway_index = list(dp_part_net.query("Pathway == '" + p + "'").index)
    all_pathway_partnet_rowidx.append(temp_pathway_index)

  #3.pathway_subnet_arr
  pathway_subnet_arr = list()
  for i in range(len(all_pathway_partnet_rowidx)):
      this_pth_rowinx = all_pathway_partnet_rowidx[i]
      for j in range(len(this_pth_rowinx)):
          this_pth_rowinx[j] = int(this_pth_rowinx[j])
      sdf = dp_part_net.loc[this_pth_rowinx,] 
      pathway_subnet_arr.append(sdf)

  return pathway_arr,pathway_subnet_arr






#二维数组排序函数（[3,5]!=[5,3]所以要先sorted再比较）
def sort_2d_array(the_array):
    for i in range(len(the_array)):
        the_array[i] = sorted(the_array[i])
    return the_array




# tool function
def directed_tool_part2(pathway_subnet_arr, pathway_arr):
    #2
    #存放各个pathway子图找到的环
    pathway_subnet_cycles_dict = dict()
    pathway_ord_cpds_dict = dict()
    pathway_ridord_cpds_dict = dict()
    import read_n_build_graph
    for i in range(len(pathway_subnet_arr)):
        this_pathway_subnet = pathway_subnet_arr[i]
        
        sub_graph_res = read_n_build_graph.read_n_build_subnet_graph_func(this_pathway_subnet)
        sub_G = sub_graph_res[0]

        import find_cycles
        sub_cycles_arr = find_cycles.find_cyc_func("directed", sub_G) 
        
        #环的顺序
        import find_cycles_order
        sub_ord_cpds_arr = find_cycles_order.pthsub_order_undirected_cycle_cpds(sub_G, sub_cycles_arr)
        sub_ridord_cpds_arr = find_cycles_order.ridorder_directed_cycle_cpds(sub_G, sub_cycles_arr)
        #对找到的每个环排序，为了后续取并集时比较两个数组 （[3,5]!=[5,3]所以要先sorted再比较）
        sub_cycles_arr = sort_2d_array(sub_cycles_arr)
        #sub_ord_cpds_arr = sort_2d_array(sub_ord_cpds_arr)

        pathway_subnet_cycles_dict[pathway_arr[i]] = sub_cycles_arr
        pathway_ord_cpds_dict[pathway_arr[i]] = sub_ord_cpds_arr
        pathway_ridord_cpds_dict[pathway_arr[i]] = sub_ridord_cpds_arr


    print("pathway_subnet_cycles_dict")
    return pathway_subnet_cycles_dict, pathway_ord_cpds_dict, pathway_ridord_cpds_dict

# tool function
def directed_tool_merge3(cycles_arr, ord_cpds_arr, ridord_cpds_arr, pathway_subnet_cycles_dict, pathway_ord_cpds_dict, pathway_ridord_cpds_dict):
    #3
    #先对每个pathway找到的所有环进行sorted，然后对二维数组取并集

    ##########################################################
    #在python先直接相加，在R中去掉重复的行即可


    #此数组存放环的并集
    merge_cycle_arr = cycles_arr
    print("merge_cycle_arr")
    print(len(merge_cycle_arr))

    #取#1和#2的环的并集
    for key in pathway_subnet_cycles_dict:
        ths_pth_cyc = pathway_subnet_cycles_dict[key]
        merge_cycle_arr = merge_cycle_arr + ths_pth_cyc  #单纯相加,后面再去重

    print("merge_cycle_arr")
    print(len(merge_cycle_arr))


    ##########################################################
    #单纯取化合物顺序的并集会产生bug
    #求差集序号

    #要将删掉的cycle环序号记录下来，在ord_cpds_arr中删掉这些序号即可
    #union_ord_cpds_arr
    #此数组存放环化合物顺序的并集
    merge_ord_cpds_arr = ord_cpds_arr
    merge_ridord_cpds_arr = ridord_cpds_arr
    print("merge_ord_cpds_arr")
    print(len(merge_ord_cpds_arr))
    print("merge_ridord_cpds_arr")
    print(len(merge_ridord_cpds_arr))
    #取#1和#2的环化合物顺序的并集
    for key in pathway_ord_cpds_dict:
        ths_pth_cyc = pathway_ord_cpds_dict[key]
        merge_ord_cpds_arr = merge_ord_cpds_arr + ths_pth_cyc

    for key in pathway_ridord_cpds_dict:
        ths_pth_cyc = pathway_ridord_cpds_dict[key]
        merge_ridord_cpds_arr = merge_ridord_cpds_arr + ths_pth_cyc


    print("merge_ord_cpds_arr")
    print(len(merge_ord_cpds_arr))

    print("merge_ridord_cpds_arr")
    print(len(merge_ridord_cpds_arr))

    return merge_cycle_arr, merge_ord_cpds_arr, merge_ridord_cpds_arr


###########################################################################################
#输入数据 并得到关于环的全部信息
#directed
def get_pathway_subnet_info_directed(input_rdata):

    import pyreadr as pyreadr
    result = pyreadr.read_r(input_rdata)

    dp_part_net = result["dp_part_net"]

    res_all_pathway_partnet_rowidx = get_all_pathway_partnet_rowidx(dp_part_net)
    pathway_arr = res_all_pathway_partnet_rowidx[0]
    pathway_subnet_arr = res_all_pathway_partnet_rowidx[1]


    #####################################################################################
    # 1. 先整体找一次环
    # 2. 再对每个pathway找一次环
    # 3. 做并集得到全部环

    #1
    import read_n_build_graph
    graph_res = read_n_build_graph.read_n_build_subnet_graph_func(dp_part_net)
    G = graph_res[0]


    import networkx as nx
    import find_cycles
    cycles_arr = find_cycles.find_cyc_func("directed", G) 
    #[['C00026', 'C00311', 'C05379'], ['C00026', 'C05381', 'C16254', 'C00091']]
    cycles_arr = sort_2d_array(cycles_arr)

    #环的顺序
    import find_cycles_order
    ord_cpds_arr = find_cycles_order.pthsub_order_undirected_cycle_cpds(G, cycles_arr)
    ridord_cpds_arr = find_cycles_order.ridorder_directed_cycle_cpds(G, cycles_arr)

    #2
    #存放各个pathway子图找到的环
    res_directed_tool_part2 = directed_tool_part2(pathway_subnet_arr, pathway_arr)
    pathway_subnet_cycles_dict = res_directed_tool_part2[0]
    pathway_ord_cpds_dict = res_directed_tool_part2[1]
    pathway_ridord_cpds_dict = res_directed_tool_part2[2]

    #3
    #先对每个pathway找到的所有环进行sorted，然后对二维数组取并集
    res_directed_tool_merge3 = directed_tool_merge3(cycles_arr, ord_cpds_arr, ridord_cpds_arr, pathway_subnet_cycles_dict, pathway_ord_cpds_dict, pathway_ridord_cpds_dict)
    merge_cycle_arr = res_directed_tool_merge3[0]
    merge_ord_cpds_arr = res_directed_tool_merge3[1]
    merge_ridord_cpds_arr = res_directed_tool_merge3[2]


    ##### 去重 #####
    import get_pathway_subnet_dedup
    #去重 只保留第一次出现的元素的序号
    dedup_res = get_pathway_subnet_dedup.dedup_and_get_diff_index(merge_cycle_arr)
    dedup_merge_cycle_arr = dedup_res[0]
    reserve_index = dedup_res[1]

    #在ord和ridord中也只保留reserve_index对应的元素
    merge_ord_cpds_arr = get_pathway_subnet_dedup.reserve_2d_array_by_index(reserve_index, merge_ord_cpds_arr)
    merge_ridord_cpds_arr = get_pathway_subnet_dedup.reserve_2d_array_by_index(reserve_index, merge_ridord_cpds_arr)
    merge_cycle_arr = dedup_merge_cycle_arr

    print("###################################")
    print(len(merge_cycle_arr))
    print(len(merge_ord_cpds_arr))
    print(len(merge_ridord_cpds_arr))
    print("###################################")


    # 此后3行因为bug 不加这三行 G 会找不到C00018
    import read_n_build_graph
    graph_res = read_n_build_graph.read_n_build_subnet_graph_func(dp_part_net)
    G = graph_res[0]

    return G, merge_cycle_arr, merge_ord_cpds_arr, merge_ridord_cpds_arr




###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################


###########################################################################################
#undirected 此函数与上面的完全相同 只是undirected 且去掉了所有的ridoreder部分 但有order部分
########################################################################################
########################################################################################


#二维数组排序函数（[3,5]!=[5,3]所以要先sorted再比较）
def sort_2d_array(the_array):
    for i in range(len(the_array)):
        the_array[i] = sorted(the_array[i])
    return the_array


# tool function
def undirected_tool_part2(pathway_subnet_arr, pathway_arr):
    #2
    #存放各个pathway子图找到的环
    pathway_subnet_cycles_dict = dict()
    pathway_ord_cpds_dict = dict()
    import read_n_build_graph
    for i in range(len(pathway_subnet_arr)):
        this_pathway_subnet = pathway_subnet_arr[i]
        
        sub_graph_res = read_n_build_graph.read_n_build_subnet_graph_func(this_pathway_subnet)
        sub_G = sub_graph_res[0]

        import find_cycles
        sub_cycles_arr = find_cycles.find_cyc_func("undirected", sub_G) 

        
        #环的顺序
        import find_cycles_order
        sub_ord_cpds_arr = find_cycles_order.pthsub_order_undirected_cycle_cpds(sub_G, sub_cycles_arr)
        #对找到的每个环排序，为了后续取并集时比较两个数组 （[3,5]!=[5,3]所以要先sorted再比较）
        sub_cycles_arr = sort_2d_array(sub_cycles_arr)
        #sub_ord_cpds_arr = sort_2d_array(sub_ord_cpds_arr)

        pathway_subnet_cycles_dict[pathway_arr[i]] = sub_cycles_arr
        pathway_ord_cpds_dict[pathway_arr[i]] = sub_ord_cpds_arr

    return pathway_subnet_cycles_dict, pathway_ord_cpds_dict


# tool function
def undirected_tool_merge3(cycles_arr, ord_cpds_arr, pathway_subnet_cycles_dict, pathway_ord_cpds_dict):
    #4
    #先对每个pathway找到的所有环进行sorted，然后对二维数组取并集
    #对#2和#1的结果取并集

    #二维数组取并集函数，(用不到了) 去重放到R中处理，去掉重复的行即可
    # def union_2d_list(a,b) :
    #     #a = [[1,2],[3,5],[7,6]]
    #     #b = [[0,9],[7,8],[5,3]]
    #     union = set(tuple(t) for t in a+b)
    #     union = list(list(t) for t in union)
    #     return union



    ##########################################################
    #在python先直接相加，在R中去掉重复的行即可


    #此数组存放环的并集
    merge_cycle_arr = cycles_arr
    print("merge_cycle_arr")
    print(len(merge_cycle_arr))

    #取#1和#2的环的并集
    for key in pathway_subnet_cycles_dict:
        ths_pth_cyc = pathway_subnet_cycles_dict[key]
        merge_cycle_arr = merge_cycle_arr + ths_pth_cyc  #单纯相加,到R中再去重


    print("merge_cycle_arr")
    print(len(merge_cycle_arr))


    ##########################################################
    #单纯取化合物顺序的并集会产生bug
    #求差集序号

    #要将删掉的cycle环序号记录下来，在ord_cpds_arr中删掉这些序号即可
    #union_ord_cpds_arr
    #此数组存放环化合物顺序的并集
    merge_ord_cpds_arr = ord_cpds_arr
    #merge_ridord_cpds_arr = ridord_cpds_arr
    print("merge_ord_cpds_arr")
    print(len(merge_ord_cpds_arr))
    #print("merge_ridord_cpds_arr")
    #print(len(merge_ridord_cpds_arr))
    #取#1和#2的环化合物顺序的并集
    for key in pathway_ord_cpds_dict:
        ths_pth_cyc = pathway_ord_cpds_dict[key]
        merge_ord_cpds_arr = merge_ord_cpds_arr + ths_pth_cyc

    # for key in pathway_ridord_cpds_dict:
    #     ths_pth_cyc = pathway_ridord_cpds_dict[key]
    #     merge_ridord_cpds_arr = merge_ridord_cpds_arr + ths_pth_cyc


    # for key in pathwaycomp_ridord_cpds_dict:
    #     ths_pth_cyc = pathwaycomp_ridord_cpds_dict[key]
    #     merge_ridord_cpds_arr = merge_ridord_cpds_arr + ths_pth_cyc

    print("merge_ord_cpds_arr")
    print(len(merge_ord_cpds_arr))

    # print("merge_ridord_cpds_arr")
    # print(len(merge_ridord_cpds_arr))

    return merge_cycle_arr, merge_ord_cpds_arr


def get_pathway_subnet_info_undirected(input_rdata):
    import pyreadr as pyreadr
    result = pyreadr.read_r(input_rdata)

    dp_part_net = result["dp_part_net"]

    res_all_pathway_partnet_rowidx = get_all_pathway_partnet_rowidx(dp_part_net)
    pathway_arr = res_all_pathway_partnet_rowidx[0]
    pathway_subnet_arr = res_all_pathway_partnet_rowidx[1]


    #####################################################################################
    # 1. 先整体找一次环
    # 2. 再对每个pathway找一次环
    # 3. 做并集得到全部环

    #1
    import read_n_build_graph
    graph_res = read_n_build_graph.read_n_build_subnet_graph_func(dp_part_net)
    G = graph_res[0]
 
    import find_cycles
    cycles_arr = find_cycles.find_cyc_func("undirected", G) 
    #[['C00026', 'C00311', 'C05379'], ['C00026', 'C05381', 'C16254', 'C00091']]
    #对找到的每个环排序，为了后续取并集时比较两个数组 （[3,5]!=[5,3]所以要先sorted再比较）
    cycles_arr = sort_2d_array(cycles_arr)


    #环的顺序
    import find_cycles_order
    ord_cpds_arr = find_cycles_order.pthsub_order_undirected_cycle_cpds(G, cycles_arr)

    #2
    #存放各个pathway子图找到的环
    res_undirected_tool_part2 = undirected_tool_part2(pathway_subnet_arr, pathway_arr)
    pathway_subnet_cycles_dict = res_undirected_tool_part2[0]
    pathway_ord_cpds_dict = res_undirected_tool_part2[1]


    #3
    res_undirected_tool_merge3 = undirected_tool_merge3(cycles_arr, ord_cpds_arr, pathway_subnet_cycles_dict, pathway_ord_cpds_dict)
    merge_cycle_arr = res_undirected_tool_merge3[0]
    merge_ord_cpds_arr = res_undirected_tool_merge3[1]

    ##### 去重 #####
    import get_pathway_subnet_dedup
    #去重 只保留第一次出现的元素的序号
    dedup_res = get_pathway_subnet_dedup.dedup_and_get_diff_index(merge_cycle_arr)
    dedup_merge_cycle_arr = dedup_res[0]
    reserve_index = dedup_res[1]

    #在ord和ridord中也只保留reserve_index对应的元素
    merge_ord_cpds_arr = get_pathway_subnet_dedup.reserve_2d_array_by_index(reserve_index, merge_ord_cpds_arr)
    merge_cycle_arr = dedup_merge_cycle_arr

    print("###################################")
    print(len(merge_cycle_arr))
    print(len(merge_ord_cpds_arr))
    print("###################################")



    # 此后3行因为bug 不加这三行 G 会找不到C00018
    import read_n_build_graph
    
    graph_res = read_n_build_graph.read_n_build_subnet_graph_func(dp_part_net)
    G = graph_res[0]

    return G, merge_cycle_arr, merge_ord_cpds_arr












###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################


#write directed
def write_directed_list_main(input_rdata, output_path, result_type):
    directed_info = get_pathway_subnet_info_directed(input_rdata)
    G = directed_info[0]
    merge_cycle_arr = directed_info[1]
    merge_ord_cpds_arr = directed_info[2]
    merge_ridord_cpds_arr = directed_info[3]
    #### 完整版信息
    if result_type == "cycle":
        import write_cycle_list
        write_cycle_list.write_complete_pathway_subnet_directed_cycle_df(G, merge_cycle_arr, merge_ord_cpds_arr, merge_ridord_cpds_arr, "directed", output_path+"/res_allpathway_cycle_union_directed.RData")
    else: # "compound"
        import write_compound_list
        write_compound_list.write_complete_compound_df(G, merge_cycle_arr, "directed", output_path+"/res_allpathway_compound_union_directed.RData")

#write undirected
def write_undirected_list_main(input_rdata, output_path, result_type):
    undirected_info = get_pathway_subnet_info_undirected(input_rdata)
    G = undirected_info[0]
    merge_cycle_arr = undirected_info[1]
    merge_ord_cpds_arr = undirected_info[2]

    #### 完整版信息
    if result_type == "cycle":
        import write_cycle_list
        write_cycle_list.write_complete_pathway_subnet_undirected_cycle_df(G, merge_cycle_arr, merge_ord_cpds_arr, "undirected", output_path+"/res_allpathway_cycle_union_undirected.RData")
    else: # "compound"
        import write_compound_list
        write_compound_list.write_complete_compound_df(G, merge_cycle_arr, "undirected", output_path+"/res_allpathway_compound_union_undirected.RData")





#######################
#test 9.25
# get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')


#############################################################################################

#### 完整版信息
#import write_cycle_list
#write_cycle_list.write_complete_pathway_subnet_directed_cycle_df(G, merge_cycle_arr, merge_ord_cpds_arr, merge_ridord_cpds_arr, "directed", "res_allpathway_cycle_union_directed.RData")
#write_cycle_list.write_complete_pathway_subnet_undirected_cycle_df(G, merge_cycle_arr, merge_ord_cpds_arr, "undirected", "res_allpathway_cycle_union_undirected.RData")


####################################################################################
#### 生成完整版信息
# write_directed_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData', "compound")
# write_undirected_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData', "compound")

# write_directed_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData', "cycle")
# write_undirected_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData', "cycle")





