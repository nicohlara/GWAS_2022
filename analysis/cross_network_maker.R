

#install.packages("igraph")

library(igraph)
setwd('C:/Users/nalara/Documents/GitHub/GWAS_2022/')

crosses <- read.delim("data/cross_info/cross_matrix.csv", sep=",")
row.names(crosses) <- crosses$X
crosses <- as.matrix(crosses[,-1])
network <- graph_from_adjacency_matrix(crosses, mode='undirected', diag=F)

#par(mfrow=c(2,2), mar=c(1,1,1,1))
#plot(network, layout=layout.circle, main="circle")

rotation <- c(seq(0, -pi, length=7),seq(pi, 0, length=6))
distance <- c(1.25,0.5,0.1,0,0,0.75,1.75,0,0,0.75,0.5,0,1.5)


png(width=3500, height=3000, filename = 'output/cross_network.png')

par(bg=NA)
plot(network, 
	layout=layout.circle, 
	
	#main="Parent Crosses", 
	
	edge.color="#495696", 
	edge.curved=.5, 
	edge.width=3, 
	
	#vertex.shape="none", 
	vertex.color="#FFE18E",
	vertex.frame.color="#ECBE74",
	vertex.label.cex=7,
	vertex.label.degree = rotation ,
	vertex.label.dist=distance ,
	vertex.label.color="#B58439"
)
dev.copy(png, "output/cross_network.png")
dev.off()
