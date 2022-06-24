# install

library(devtools)  
install_github("zxzhangragnar/CycleFlux")  

## Refer to the data of folder 'example_input'
net: 'example_input/hsa_net.RData'  
(the input metabolic network)
gene_deg: 'example_input/TCGA_deg_genes.RData'  
(no_express: -10 up: 1 gap: -1 no_change: 0)  

## use '?getCycleFlux' to view the meaning of each parameter
library(CycleFlux)  
?getCycleFlux  

# parameter
'draw_single_graph':
If this parameter is TRUE, an image will be generated in this directory for each cycle that matches the shift definition.

'draw_net_graph':
If this parameter is TRUE, an image containing all cycles that match the shift definition will be generated in this directory.


# function
## step 1
load('example_input/hsa_net.RData')  
basic_cycle = getBasicCycle(net)  
#save(basic_cycle, file="example_input/basic_cycle.RData")  

## step2
load('example_input/TCGA_deg_genes.RData')  
#load("example_input/basic_cycle.RData")  
result = getCycleFlux(basic_cycle, gene_deg, FALSE, FALSE)  
