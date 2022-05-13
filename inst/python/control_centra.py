

############################### 此函数用于生成全部结果 ##########################################




####################################### Directed ###########################################
####################################### Directed ###########################################
####################################### Directed ###########################################


def write_directed(input_net_file, output_path):
    ### init ###
    #input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData'
    #output_path = "RData_output"

    import get_pathway_subnet
    #cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
    cyc_res = get_pathway_subnet.get_pathway_subnet_info_directed(input_net_file)
    G = cyc_res[0]
    cycles_arr = cyc_res[1]

    ### 1 ###
    import get_pathway_subnet
    get_pathway_subnet.write_directed_list_main(input_net_file, output_path, "compound")
    get_pathway_subnet.write_directed_list_main(input_net_file, output_path, "cycle")

    import subnet_edge_expression
    subnet_edge_expression.write_subnet_edge_expression(G, output_path, "directed")


    # ## 2 ###
    ## 2步度
    import cycle_degnode_successors
    cycle_degnode_successors.write_situation_1(G, cycles_arr, output_path, "directed")
    cycle_degnode_successors.write_situation_2(G, cycles_arr, output_path, "directed")
    cycle_degnode_successors.write_situation_3(G, cycles_arr, output_path, "directed")

    ## 1步度
    import cycle_degnode_edgesucs
    cycle_degnode_edgesucs.write_edgesucs_out(G, cycles_arr, output_path, "directed")
    cycle_degnode_edgesucs.write_edgesucs_in(G, cycles_arr, output_path, "directed")

    ### 3 ###
    ## 环整体 表达
    import cycle_expression
    cycle_expression.write_cycle_expression(G, cycles_arr, output_path, "directed")

    ### 4 ###
    ## 环的每条边 表达
    import cycle_edge_expression
    cycle_edge_expression.write_cycle_edge_expression(G, cycles_arr, output_path, "directed")

    # ### 5 ###
    import cycle_deg_edge_expression
    ## 2步度 表达
    cycle_deg_edge_expression.write_situation_comps_1(G, cycles_arr, output_path, "directed")
    cycle_deg_edge_expression.write_situation_comps_2(G, cycles_arr, output_path, "directed")
    cycle_deg_edge_expression.write_situation_comps_3(G, cycles_arr, output_path, "directed")

    ## 1步度 表达
    cycle_deg_edge_expression.write_edgesucs_out_comps(G, cycles_arr, output_path, "directed")
    cycle_deg_edge_expression.write_edgesucs_in_comps(G, cycles_arr, output_path, "directed")

    
# test
# input_net_file = 'E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/subnet_input/hsa_subnet.RData'
# output_path = "E:/scFEA_universal/my_R/Rpackage/cycle_flux/output_files"
# write_directed(input_net_file, output_path)


####################################### unDirected ###########################################
####################################### unDirected ###########################################
####################################### unDirected ###########################################



# def write_undirected():
#     ### init ###
#     import get_pathway_subnet
#     cyc_res = get_pathway_subnet.get_pathway_subnet_info_undirected('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData')
#     G = cyc_res[0]
#     cycles_arr = cyc_res[1]

#     ### 1 ###
#     import get_pathway_subnet
#     get_pathway_subnet.write_undirected_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData', "compound")
#     get_pathway_subnet.write_undirected_list_main('E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/main_input/hsa_net.RData', "cycle")

#     ### 2 ###
#     import cycle_degnode_successors
#     cycle_degnode_successors.write_situation_1(G, cycles_arr, "undirected")
#     cycle_degnode_successors.write_situation_2(G, cycles_arr, "undirected")
#     cycle_degnode_successors.write_situation_3(G, cycles_arr, "undirected")

#     ### 3 ###
#     import cycle_expression
#     cycle_expression.write_cycle_expression(G, cycles_arr, "undirected")

# test
# write_undirected()






