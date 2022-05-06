#####################################################################################

# 对环cyc和ord、ridord 的去重还是应该在 python中 完成

# 3者 先用单纯相加，得到总的结果，3者的长度应相同 且一一对应
# 再 先对环去重 得到去重后的结果 记录下来不同的元素的序号idxs

# 再在 ord、ridord中删去同样的序号idxs
# 完成去重



##########################
# tools 

#去重 只保留第一次出现的元素的序号
def dedup_and_get_diff_index(arr):
    # 记录2者不同的元素的序号
    set_collection = list()
    index_collection = list()
    for i in range(len(arr)):
        if arr[i] not in set_collection:
            set_collection.append(arr[i])
            index_collection.append(i)
    
    return set_collection, index_collection



# 按序号保留2维数组中的元素
def reserve_2d_array_by_index(reserve_index, a):
    # diff_index = [0,1]
    reserve_arr = [a[i] for i in range(len(a)) if(i in reserve_index)]
    return reserve_arr




###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

#test
# a = [[3,4],[5,6],[8,9],[10,11]]   
# b = [[3,4],[10,11]]   
# diff_index = get_diff_index(a, b)
# print(diff_index)
# diff_arr = del_2d_array_by_index(diff_index, a)
# print(diff_arr)



#test dedup_get_diff_index
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# merge_cycle_arr = cyc_res[1]

# set_collection = list()
# index_collection = list()
# for i in range(len(merge_cycle_arr)):
#     if merge_cycle_arr[i] not in set_collection:
#         set_collection.append(merge_cycle_arr[i])
#         index_collection.append(i)
# print(len(set_collection))
# print(len(index_collection))

