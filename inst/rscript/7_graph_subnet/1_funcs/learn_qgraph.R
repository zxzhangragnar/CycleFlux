# igraph: Resolving tight overlapping nodes
# https://stackoverflow.com/questions/39290909/igraph-resolving-tight-overlapping-nodes


library(igraph)
library(qgraph)

g <- barabasi.game(355, directed=FALSE)


png("plot1.png", height=6, width=12, units="in", res=250)
par(mfrow=c(1, 3))

plot(g,layout=layout_with_fr,vertex.size=4,vertex.label=NA)
mtext("layout_with_fr", side=1)

e <- get.edgelist(g,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g))
plot(g,layout=l,vertex.size=4,vertex.label=NA)
mtext("qgraph.layout.fruchtermanreingold default", side=1)

l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                       area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
plot(g,layout=l,vertex.size=4,vertex.label=NA)
mtext("qgraph.layout.fruchtermanreingold modified", side=1)

dev.off()
