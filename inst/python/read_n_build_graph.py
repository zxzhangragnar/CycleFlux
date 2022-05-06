# input： RData文件名 和 数据集名  output：图 G 
# 读数据并建图
def read_n_build_graph_func(RData_name,dataset_name):
    print("create entire net")
    #读RData数据
    import pyreadr as pyreadr
    result = pyreadr.read_r(RData_name)
    input_dataframe = result[dataset_name]

    res_graph = build_graph_core(input_dataframe)
    G = res_graph[0]
    cin_arr = res_graph[1]
    cout_arr = res_graph[2]

    return G,cin_arr,cout_arr


#针对某个pathway建立子图，参考get_pathway_subnet.py 和 get_pathway_subnet.R
def read_n_build_subnet_graph_func(pathway_subnet):
    print("create subnet")
    input_dataframe = pathway_subnet

    res_graph = build_graph_core(input_dataframe)
    G = res_graph[0]
    cin_arr = res_graph[1]
    cout_arr = res_graph[2]

    return G,cin_arr,cout_arr


# 通过读取df，完成重新建图
#def read_n_build_bydf_func(df)


def add_edge_attribute(G, hnode, tnode, edge_ec_temp, edge_reaction_id_temp, edge_pathway_temp, edge_reaction_formula_temp, edge_gene_symbol_temp):
    
    # EC
    # Reaction_id
    # Pathway
    # Reaction_formula
    # Gene_symbol
    
    #老的edge的值
    edge_ec_old = G[hnode][tnode]['EC']
    edge_reaction_id_old = G[hnode][tnode]['Reaction_id']
    edge_pathway_old = G[hnode][tnode]['Pathway']
    edge_reaction_formula_old = G[hnode][tnode]['Reaction_formula']
    edge_gene_symbol_old = G[hnode][tnode]['Gene_symbol']

    #新的edge的值
    edge_ec_new = edge_ec_old + ";" + edge_ec_temp
    edge_reaction_id_new = edge_reaction_id_old + ";" + edge_reaction_id_temp
    edge_pathway_new = edge_pathway_old + ";" + edge_pathway_temp
    edge_reaction_formula_new = edge_reaction_formula_old + ";" + edge_reaction_formula_temp
    edge_gene_symbol_new = edge_gene_symbol_old + ";" + edge_gene_symbol_temp

    #去重
    edge_ec_new_arr = edge_ec_new.split(";")
    edge_ec_new_arr = list(set(edge_ec_new_arr))
    edge_ec_new_str = ";".join(edge_ec_new_arr)
    
    edge_reaction_id_new_arr = edge_reaction_id_new.split(";")
    edge_reaction_id_new_arr = list(set(edge_reaction_id_new_arr))
    edge_reaction_id_new_str = ";".join(edge_reaction_id_new_arr) #存放这条边edge包含的全部反应Rid

    edge_pathway_new_arr = edge_pathway_new.split(";")
    edge_pathway_new_arr = list(set(edge_pathway_new_arr))
    edge_pathway_new_str = ";".join(edge_pathway_new_arr)

    edge_reaction_formula_new_arr = edge_reaction_formula_new.split(";")
    edge_reaction_formula_new_arr = list(set(edge_reaction_formula_new_arr))
    edge_reaction_formula_new_str = ";".join(edge_reaction_formula_new_arr)

    edge_gene_symbol_new_arr = edge_gene_symbol_new.split(";")
    edge_gene_symbol_new_arr = list(set(edge_gene_symbol_new_arr))
    edge_gene_symbol_new_str = ";".join(edge_gene_symbol_new_arr)

    return edge_ec_new_str, edge_reaction_id_new_str, edge_pathway_new_str, edge_reaction_formula_new_str, edge_gene_symbol_new_str


## 工具函数
def build_graph_core(input_df):
    df = input_df
    #从读取入的数据中，拆分出不同列名的数据   
    Rid = df["Rid"]
    EC = df["EC"]
    Irreversible = df["Irreversible"]
    Pathway = df["Pathway"]
    Reaction_formula = df["Reaction_formula"]
    cin = df["C_in"]
    cout = df["C_out"]
    Gene_symbol = df["Gene_symbol"]
    
    #for expression
    Rid_arr = Rid.array
    EC_arr = EC.array
    Irreversible_arr = Irreversible.values

    Pathway_arr = Pathway.array
    Reaction_formula_arr = Reaction_formula.array
    cin_arr = list(cin.array)
    cout_arr = list(cout.array)
    Gene_symbol_arr = Gene_symbol.array

    #print(str(Irreversible_arr[0]) == "True")
    #print(cin_arr[0]) #C00404;C00001

    ###########################################################################
    #建图，向图中添加节点
    print("build graph")
    import networkx as nx

    G = nx.DiGraph()
    #1.向有向图中添加节点，并设置节点node_id信息，确保生成的图节点顺序相同
    #G.add_node(1, time="5pm")
    # 全部节点 all_node = unique(cin_arr + cout_arr) 
    merge_node = cin_arr + cout_arr
    all_node = list(set(merge_node)) # 用set去重（去重后顺序是乱的）
    all_node.sort(key=merge_node.index) #加上本句顺序就不会乱了
    for i in range(len(all_node)):
        G.add_node(all_node[i], node_id = i)


    #2.向有向图中添加边
    #G.add_edges_from([(2, 3, {"weight": 8})])
    #print(G.edges[2, 3]["weight"])
    for i in range(len(Rid_arr)):

        Rid_symbol = (Rid_arr[i]).replace("R","E") #这个值是这条边的代号 而不是指的某个真正的反应 R00001->E00001
        #情况1: 已存在这条边了, pathway和gene和EC要用append
        if G.has_edge(cin_arr[i], cout_arr[i]):  #注意 G[a][b]['val'] != G[b][a]['val']
            new_edge_attribute = add_edge_attribute(G, cin_arr[i], cout_arr[i], EC_arr[i], Rid_arr[i], Pathway_arr[i], Reaction_formula_arr[i], Gene_symbol_arr[i])
            edge_ec_new_str = new_edge_attribute[0]
            edge_reaction_id_new_str = new_edge_attribute[1]
            edge_pathway_new_str = new_edge_attribute[2]
            edge_reaction_formula_new_str = new_edge_attribute[3]
            edge_gene_symbol_new_str = new_edge_attribute[4]
            G.add_edges_from([(cin_arr[i], cout_arr[i], {"Rid": Rid_symbol, "EC": edge_ec_new_str, "Reaction_id": edge_reaction_id_new_str, "Pathway": edge_pathway_new_str, "Reaction_formula": edge_reaction_formula_new_str, "Gene_symbol": edge_gene_symbol_new_str})])
            
            # True or na
            if (str(Irreversible_arr[i]) == "True") or (str(Irreversible_arr[i]) != "True" and str(Irreversible_arr[i]) != "False"):
                if G.has_edge(cout_arr[i], cin_arr[i]): 
                    irv_new_edge_attribute = add_edge_attribute(G, cout_arr[i], cin_arr[i], EC_arr[i], Rid_arr[i], Pathway_arr[i], Reaction_formula_arr[i], Gene_symbol_arr[i])
                    irv_edge_ec_new_str = irv_new_edge_attribute[0]
                    irv_edge_reaction_id_new_str = irv_new_edge_attribute[1]
                    irv_edge_pathway_new_str = irv_new_edge_attribute[2]
                    irv_edge_reaction_formula_new_str = irv_new_edge_attribute[3]
                    irv_edge_gene_symbol_new_str = irv_new_edge_attribute[4]
                    G.add_edges_from([(cout_arr[i], cin_arr[i], {"Rid": Rid_symbol, "EC": irv_edge_ec_new_str, "reaction_id": irv_edge_reaction_id_new_str, "Pathway": irv_edge_pathway_new_str, "Reaction_formula": irv_edge_reaction_formula_new_str, "Gene_symbol": irv_edge_gene_symbol_new_str})])
                else:
                    G.add_edges_from([(cout_arr[i], cin_arr[i], {"Rid": Rid_symbol, "EC": EC_arr[i], "Reaction_id": Rid_arr[i], "Pathway": Pathway_arr[i], "Reaction_formula": Reaction_formula_arr[i], "Gene_symbol": Gene_symbol_arr[i]})])


        #情况2: 不存在这条边, 第一次为这条边添加值
        else:
            G.add_edges_from([(cin_arr[i], cout_arr[i], {"Rid": Rid_symbol, "EC": EC_arr[i], "Reaction_id": Rid_arr[i], "Pathway": Pathway_arr[i], "Reaction_formula": Reaction_formula_arr[i], "Gene_symbol": Gene_symbol_arr[i]})])
            
            # True or na
            if (str(Irreversible_arr[i]) == "True") or (str(Irreversible_arr[i]) != "True" and str(Irreversible_arr[i]) != "False"):
                if G.has_edge(cout_arr[i], cin_arr[i]): 
                    irv_new_edge_attribute = add_edge_attribute(G, cout_arr[i], cin_arr[i], EC_arr[i], Rid_arr[i], Pathway_arr[i], Reaction_formula_arr[i], Gene_symbol_arr[i])
                    irv_edge_ec_new_str = irv_new_edge_attribute[0]
                    irv_edge_reaction_id_new_str = irv_new_edge_attribute[1]
                    irv_edge_pathway_new_str = irv_new_edge_attribute[2]
                    irv_edge_reaction_formula_new_str = irv_new_edge_attribute[3]
                    irv_edge_gene_symbol_new_str = irv_new_edge_attribute[4]
                    G.add_edges_from([(cout_arr[i], cin_arr[i], {"Rid": Rid_symbol, "EC": irv_edge_ec_new_str, "reaction_id": irv_edge_reaction_id_new_str, "Pathway": irv_edge_pathway_new_str, "Reaction_formula": irv_edge_reaction_formula_new_str, "Gene_symbol": irv_edge_gene_symbol_new_str})])
                else:
                    G.add_edges_from([(cout_arr[i], cin_arr[i], {"Rid": Rid_symbol, "EC": EC_arr[i], "Reaction_id": Rid_arr[i], "Pathway": Pathway_arr[i], "Reaction_formula": Reaction_formula_arr[i], "Gene_symbol": Gene_symbol_arr[i]})])
            

    return G,cin_arr,cout_arr








#########################################################
#test
#问题：运行两次print(G.nodes)发现每次输出结果也不一样
#已解决：原因直接用set去重顺序会乱
# graph_res = read_n_build_graph_func('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/all_pathway_subnet.RData','dp_part_net')
# G = graph_res[0]
# cin_arr = graph_res[1]
# cout_arr = graph_res[2]
# G_nodes = list(G.nodes)
# G_node_dict = list(G.nodes.data("node_id"))

# print(G_nodes)
# print(G_node_dict)


