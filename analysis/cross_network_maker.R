

#install.packages("igraph")

library(igraph)

crosses <- read.delim("/Users/nico/Documents/GitHub/GWAS_2022/data/cross_info/cross_matrix.csv", sep=",")
row.names(crosses) <- crosses$X
crosses <- as.matrix(crosses[,-1])
network <- graph_from_adjacency_matrix(crosses, mode='undirected', diag=F)

#par(mfrow=c(2,2), mar=c(1,1,1,1))
#plot(network, layout=layout.circle, main="circle")

rotation <- c(seq(0, -pi, length=7),seq(pi, 0, length=6))
distance <- c(1.25,0.5,0.1,0,0,0.75,1.75,0,0,0.75,0.5,0,1.5)


png(width=1500, height=1500, filename = '/Users/nico/Documents/GitHub/GWAS_2022/output/cross_network.png')

par(bg=NA)
plot(network, 
	layout=layout.circle, 
	
	#main="Parent Crosses", 
	
	edge.color="#67A247", 
	edge.curved=.5, 
	edge.width=3, 
	
	#vertex.shape="none", 
	vertex.color="#E9CA4744",
	vertex.frame.color="#E9CA4744",
	vertex.label.cex=.9,
	vertex.label.degree = rotation ,
	vertex.label.dist=distance ,
	vertex.label.color="black"
)
dev.copy(png, "/Users/nico/Documents/GitHub/GWAS_2022/output/cross_network.png")
dev.off()