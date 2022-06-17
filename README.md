# install

library(devtools)  
install_github("zxzhangragnar/CycleFlux")  


# parameter

## Refer to the data of folder 'example_input'
net = 'example_input/hsa_net.RData'  
deg_gene = 'example_input/TCGA_deg_genes.RData'  
(no_express: -10 up: 1 gap: -1 no_change: 0)  

## use '?getCycleFlux' to view the meaning of each parameter
library(CycleFlux)  
?getCycleFlux  

'single_graph':
If this parameter is TRUE, an image will be generated in this directory for each cycle that matches the shift definition.

'net_graph':
If this parameter is TRUE, an image containing all cycles that match the shift definition will be generated in this directory.


# function

## step 1
net = 'example_input/hsa_net.RData'  
basic_cycle_result = getBasicCycle(net)  
#save(basic_cycle_result, file="example_input/basic_cycle_result.RData")  

## step2
deg_gene = 'example_input/TCGA_deg_genes.RData'  
#load("example_input/basic_cycle_result.RData")  
shift_result = getCycleFlux(basic_cycle_result, deg_gene, FALSE, FALSE)  
