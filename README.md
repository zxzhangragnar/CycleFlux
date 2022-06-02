# install

library(devtools)

install_github("zxzhangragnar/CycleFlux")

library(CycleFlux)

?getCycleFlux

# parameter

## Refer to the data of folder 'example_input'

net = './hsa_net.RData'

stat_gene = './TCGA_stat_genes.RData'

deg_gene = './TCGA_deg_genes.RData'
(no express: 00 up: 01 gap: 10 normal: 11)



## use '?getCycleFlux' to view the meaning of each parameter

prm_1=0.05

prm_2=1

prm_3=1

# function

shift_result = getCycleFlux(net, stat_gene, deg_gene, 0.05, 1, 1)
