# A hack to allow arrows size in R igraph to match edge width
# https://stackoverflow.com/questions/16942553/a-hack-to-allow-arrows-size-in-r-igraph-to-match-edge-width


library(igraph)

# (almost) your example data
d <- data.frame(start=c("a","a","b","c"),
                end=c("b","b","c","b"),
                size=1:4)

graph <- graph.data.frame(d, directed=TRUE)

## The plotting function
eqarrowPlot <- function(graph, layout, edge.lty=rep(1, ecount(graph)),
                        edge.arrow.size=rep(1, ecount(graph)),
                        vertex.shape="circle",
                        edge.curved=autocurve.edges(graph), ...) {
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape="none")
  for (e in seq_len(ecount(graph))) {
    graph2 <- delete.edges(graph, E(graph)[(1:ecount(graph))[-e]])
    plot(graph2, edge.lty=edge.lty[e], edge.arrow.size=edge.arrow.size[e],
         edge.curved=edge.curved[e], layout=layout, vertex.shape="none",
         vertex.label=NA, add=TRUE, ...)
  }
  plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
       vertex.shape=vertex.shape, add=TRUE, ...)
  invisible(NULL)
}

## Test
eqarrowPlot(graph, layout.auto(graph), edge.arrow.size=E(graph)$size/3,
            edge.width=E(graph)$size*5)






