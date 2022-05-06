

#plot
# cycle_graph <- graph.data.frame(d = graph_cyc_df, directed = FALSE)
# plot(cycle_graph, vertex.label = V(cycle_graph)$name)


##################################### igraph learn ################################################
#############################################################################################

# igraph
# g <- make_ring(10)
# g <- add_edges(g, c(1,2, 2,3, 1,3))
# coreness(g) 		# small core triangle in a ring
# 
# library(igraph)
# 
# df <- data.frame(
#   A = c("Berlin", "Amsterdam", "New York") , 
#   B = c("Munich", "Utrecht", "Chicago")
# )
# df.g <- graph.data.frame(d = df, directed = FALSE)
# 
# plot(df.g, vertex.label = V(df.g)$name)









##################################### igraph learn ################################################
#############################################################################################

#https://stackoverflow.com/questions/54571716/plotting-graph-with-igraph-in-r-edge-length-proportional-to-weight

library(igraph)

# toy data
d = data.frame("n1"=c('A','A','A','A'), "n2"=c('B','C','D','E'), "weight"=c(1,1.5,2,5))

g <- graph_from_data_frame(d)

# we can create a layout object that is just coordinate positions
coords <- layout_(g, as_star())
coords
#>               [,1]          [,2]
#> [1,]  0.000000e+00  0.000000e+00
#> [2,]  1.000000e+00  0.000000e+00
#> [3,]  6.123234e-17  1.000000e+00
#> [4,] -1.000000e+00  1.224647e-16
#> [5,] -1.836970e-16 -1.000000e+00

# Knowing this is a star the N nodes should have N-1 edges. We can scale
# The N-1 positions by the weights
weight.scale <- c(1, d$weight)

coords2 <- weight.scale * coords
coords2
#>               [,1]          [,2]
#> [1,]  0.000000e+00  0.000000e+00
#> [2,]  1.000000e+00  0.000000e+00
#> [3,]  9.184851e-17  1.500000e+00
#> [4,] -2.000000e+00  2.449294e-16
#> [5,] -9.184851e-16 -5.000000e+00

# we can inspect
plot(g, layout = coords2)



##################################### igraph learn ################################################
load("E:/scFEA_universal/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
d = cyc_all_distance_df
#距离矩阵
distance_scale <- d$distance
#ll =matrix(c(0,0,0,1,0,3,0,5),ncol=652,byrow=TRUE)
ll =matrix(distance_scale, ncol=652,byrow=TRUE)
plot(g,layout=ll)
# 按距离聚类
# Clustering with a distance matrix 根据距离矩阵进行聚类
# 详见：https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix
# 距离矩阵



##################################### igraph learn ################################################

#https://stackoverflow.com/questions/60234916/about-how-to-put-weight-as-distance-in-network-analysis-using-tidygraph

# Library
library(igraph)
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(tidygraph)
#> Attaching package: 'tidygraph'
#> The following object is masked from 'package:igraph':
#> 
#>     groups
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(ggraph)
#> Loading required package: ggplot2

set.seed(1)

load("E:/OneDrive - jlu/scFEA_Universal_cloud/my_R/aimA/rdata_cycle_detect/1_cycle_topology/result_topo/cyc_all_distance_df.RData")
#cyc_all_distance_df = cyc_all_distance_df[which(cyc_all_distance_df$from == 0),]
colnames(cyc_all_distance_df) = c("source", "target", "weight")
links = cyc_all_distance_df
# Create data
# links <- data.frame(
#   source = c("A","A", "A", "A", "A","J", "B", "B", "C", "C", "D","I"),
#   target = c("B","B", "C", "D", "J","A","E", "F", "G", "H", "I","I"),
#   weight = sample(1:20, 12, replace = T)
# )
nodes <- data.frame(
  name = c(0:651)
  #name = LETTERS[1:10]
)
network <- graph_from_data_frame(d = links, vertices = nodes, directed = F) 

g <- as_tbl_graph(network)
g %>%
  ggraph(layout = "lgl") +
  geom_edge_link(aes(width = 1),
                 alpha = 0.8,
                 colour = "lightgray") +
  scale_edge_width(range = c(0.1, 1)) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph()




#选出某一条边
g <- make_ring(10)
plot(g)
ei <- get.edge.ids(g, c(1,20))
E(g)[ei]










## 正反方向箭头 颜色不同
g <- make_ring(3)
g <- add_edges(g, c(1,2), color = "cyan")
g <- add_edges(g, c(2,1), color = "green")
plot(g)
