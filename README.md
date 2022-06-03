# install

library(devtools)

install_github("zxzhangragnar/CycleFlux")

# parameter

## Refer to the data of folder 'example_input'

net = './hsa_net.RData'

stat_gene = './TCGA_stat_genes.RData'

deg_gene = './TCGA_deg_genes.RData'
(no_express: -10 up: 1 gap: -1 no_change: 0)



## use '?getCycleFlux' to view the meaning of each parameter

library(CycleFlux)

?getCycleFlux

'mode':
IF mode=1, For an edge, if any gene on the edge is up, the edge is up, if the edge contains both the gap gene and the up gene, it is regarded as a gap.
IF mode=2, For an edge, if the gene with the largest mean of tumor value is up(gap), then the edge is up(gap).

'single_graph':
If this parameter is TRUE, an image will be generated in this directory for each cycle that matches the shift definition.

'net_graph':
If this parameter is TRUE, an image containing all cycles that match the shift definition will be generated in this directory.


# function

shift_result = getCycleFlux(net, stat_gene, deg_gene, 1, TRUE, FALSE)
