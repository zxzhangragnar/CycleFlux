library(igraph)

## (almost) your example data
d <- data.frame(start=c("a","a","b","c"),
                end=c("b","b","c","b"))
g <- graph.data.frame(d, directed=TRUE)
E(g)$color = "black"
g <- add_edges(g, c("a","b"), color = "orange")

ei <- get.edge.ids(g, c("a","b"))
E(g)[ei]$color = "cyan"

E(g)$size = 1
E(g)[ei]$size = 5

###
plot(g, edge.width=E(g)$size)




###
#找不存在的边
ei <- get.edge.ids(g, c("e","f"))
are.connected(g, "a","c")
are.connected(g, "a","b")












###############################################################################
library(igraph)
set.seed(123)

cor.matrix <- matrix(runif(100, -1, 1), nrow=10)

t = which(abs(cor.matrix) > 0.6 & lower.tri(cor.matrix),arr.ind=TRUE)
t <- cbind(t, cor.matrix[which(abs(cor.matrix) > 0.6 & lower.tri(cor.matrix),arr.ind=TRUE)]) ##this adds the correlation to the graph as an edge attribute "V3"
t.graph=graph.data.frame(t,directed=F)
E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'magenta','green') #You had this as "V3 > 0.6" which I guess works but it is more readable as 0. that way if you decide to lower the correlation threshold you do not have to change this line too.

#t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
minC <- rep(-Inf, vcount(t.graph))
maxC <- rep(Inf, vcount(t.graph))
minC[1] <- maxC[1] <- 0
l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                    miny=minC, maxy=maxC)      
plot(t.graph, layout=l, 
     rescale=T,
     asp=0,
     edge.arrow.size=0.5, 
     vertex.label.cex=0.8, 
     vertex.label.font=2,
     #vertex.label=t.names,
     vertex.shape="circle", 
     vertex.size=3, 
     vertex.color="deepskyblue2",
     vertex.label.color="black", 
     #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
     edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))
