
################################# pathway #####################################
###############################################################################
#原cycle_expression_pathway内容

#得到每个环所包含的反应pathway的字典
def get_cycles_pathway_arr(G, cycles_arr, direct):
    #将所有环中的化合物排序，只取一个方向的顺序
    import find_cycles_order
    if direct == "directed":
        cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)
    else:
        cycle_order_path_arr = find_cycles_order.order_undirected_cycle_cpds(G, cycles_arr)
    cycles_pathway_arr = list()
    for circle in cycle_order_path_arr:
        #cyc = circle[0]
        circle_str_arr = list()
        for k in range(len(circle)):
            cyc = circle[k]
            SG = G.subgraph(cyc)
            if direct == "undirected":
                SG = SG.to_undirected()
            for j in range(1,len(cyc)):
                pre_cpd = str(cyc[j-1])
                now_cpd = str(cyc[j])
                Pathway = SG[pre_cpd][now_cpd]['Pathway']
                #circle_str_arr.append(Pathway)
                # append for array
                Pathway_arr = Pathway.split(';')
                for p in Pathway_arr:
                    circle_str_arr.append(p)
        cycles_pathway_arr.append(circle_str_arr)

    return cycles_pathway_arr

#得到每个环所包含的反应Rid的字典


########################分析 cycle的pathway表达量#####################################

#统计每个环中各个基因gene的表达量
def get_cycles_pathway_dict(G, cycles_arr, direct):
    cycles_pathway_arr = get_cycles_pathway_arr(G, cycles_arr, direct)
    #统计列表中不同元素的个数
    from collections import Counter
    #651: {'SULT2B1': 2, 'STS': 1, 'CHST15': 1}
    cycles_pathway_dict = dict()
    for i in range(len(cycles_pathway_arr)):
        pathwaydict = Counter(cycles_pathway_arr[i])
        cycles_pathway_dict[i] = dict(pathwaydict)

    return cycles_pathway_dict




################################# ec ##########################################
###############################################################################
#原cycle_expression_ec内容

#得到每个环所包含的反应EC的字典
def get_cycles_ec_arr(G, cycles_arr, direct):
    #将所有环中的化合物排序，只取一个方向的顺序
    import find_cycles_order
    if direct == "directed":
        cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)
    else:
        cycle_order_path_arr = find_cycles_order.order_undirected_cycle_cpds(G, cycles_arr)
    cycles_ec_arr = list()
    for circle in cycle_order_path_arr:
        #cyc = circle[0]
        circle_str_arr = list()
        for k in range(len(circle)):
            cyc = circle[k]
            SG = G.subgraph(cyc)
            if direct == "undirected":
                SG = SG.to_undirected()
            for j in range(1,len(cyc)):
                pre_cpd = str(cyc[j-1])
                now_cpd = str(cyc[j])
                EC = SG[pre_cpd][now_cpd]['EC']
                #circle_str_arr.append(EC)
                # append for array
                EC_arr = EC.split(';')
                for e in EC_arr:
                    circle_str_arr.append(e)
        cycles_ec_arr.append(circle_str_arr)

    return cycles_ec_arr

#得到每个环所包含的反应Rid的字典



########################分析 cycle的ec表达量#####################################

#统计每个环中各个基因gene的表达量
def get_cycles_ec_dict(G, cycles_arr, direct):
    cycles_ec_arr = get_cycles_ec_arr(G, cycles_arr, direct)
    #统计列表中不同元素的个数
    from collections import Counter
    #651: {'SULT2B1': 2, 'STS': 1, 'CHST15': 1}
    cycles_ec_dict = dict()
    for i in range(len(cycles_ec_arr)):
        ecdict = Counter(cycles_ec_arr[i])
        cycles_ec_dict[i] = dict(ecdict)

    return cycles_ec_dict







################################# gene ########################################
###############################################################################

#得到每个环所包含的基因gene
def get_cycles_gene_arr(G, cycles_arr, direct):
    #将所有环中的化合物排序，只取一个方向的顺序
    import find_cycles_order
    if direct == "directed":
        cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)
    else:
        cycle_order_path_arr = find_cycles_order.order_undirected_cycle_cpds(G, cycles_arr)
    cycles_gene_arr = list()
    for circle in cycle_order_path_arr:
        #print(circle)
        #cyc = circle[0]
        circle_str_arr = list()
        for k in range(len(circle)):
            cyc = circle[k]
            SG = G.subgraph(cyc)
            if direct == "undirected":
                SG = SG.to_undirected()
            for j in range(1,len(cyc)):
                pre_cpd = str(cyc[j-1])
                now_cpd = str(cyc[j])
                Gene_symbol = SG[pre_cpd][now_cpd]['Gene_symbol']
                Gene_symbol_arr = Gene_symbol.split(';')
                for gene in Gene_symbol_arr:
                    circle_str_arr.append(gene)

        cycles_gene_arr.append(circle_str_arr)

    return cycles_gene_arr

#得到每个环所包含的基因gene




########################分析 cycle的gene表达量#####################################

#统计每个环中各个基因gene的表达量
def get_cycles_gene_dict(G, cycles_arr, direct):
    cycles_gene_arr = get_cycles_gene_arr(G, cycles_arr, direct)
    #统计列表中不同元素的个数
    from collections import Counter
    #651: {'SULT2B1': 2, 'STS': 1, 'CHST15': 1}
    cycles_gene_dict = dict()
    for i in range(len(cycles_gene_arr)):
        genedict = Counter(cycles_gene_arr[i])
        cycles_gene_dict[i] = dict(genedict)

    return cycles_gene_dict




############################# write RData ############################################
#写入RData
def write_cycle_expression(G, cycles_arr, output_path, direct):
    
    import pyreadr as pyreadr
    import pandas as pd
    list_df = pd.DataFrame()  #得到的列表
    #gene
    cycles_gene_dict = get_cycles_gene_dict(G, cycles_arr, direct)
    #pathway
    cycles_pathway_dict = get_cycles_pathway_dict(G, cycles_arr, direct)
    #ec
    cycles_ec_dict = get_cycles_ec_dict(G, cycles_arr, direct)

    for cyc in cycles_gene_dict:
        cycid = cyc
        gene_express = cycles_gene_dict[cyc]
        pathway_express = cycles_pathway_dict[cyc]
        ec_express = cycles_ec_dict[cyc]
        #添加行
        cyc_df = pd.DataFrame([[cycid, gene_express, ec_express, pathway_express]], columns=["cycid", "gene_express", "ec_express", "pathway_express"])
        list_df = list_df.append(cyc_df, ignore_index=True)

    #写入文件.RData
    import pyreadr as pyreadr
    if direct == "directed":
        pyreadr.write_rdata(output_path + "/cycle_expression.RData", list_df, df_name=str("cycle_expression"))
    else:
        pyreadr.write_rdata(output_path + "/cycle_expression.RData", list_df, df_name=str("cycle_expression"))

######################################################################################################
#test
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]
# write_cycle_expression(G, cycles_arr, "undirected")




### 4.16
# import get_pathway_subnet
# cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData')
# G = cyc_res[0]
# cycles_arr = cyc_res[1]
# import find_cycles_order
# cycle_order_path_arr = find_cycles_order.order_directed_cycle_cpds(G, cycles_arr)

# print(cycle_order_path_arr)

