#Create plots from GWAS data
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-7-11

library(tidyverse)
library(gaston)
#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#Read in GWAS file
GWAS_results <- read.delim("output/data/SunRILs_gwas_2022.csv", sep=",")

manhattan_plot <- function(graph_name, SNP, p_value) {
	dataframe <- data.frame(id = SNP, p = p_value)
	dataframe <- dataframe[complete.cases(dataframe),]
	#Add new columns for staticstical metrics, chromosomes, etc.
	dataframe$LOG <- -log10(dataframe$p)
	dataframe$chr <- str_replace(str_replace(dataframe$id, "^S", ""), "_\\d*$", "")
	dataframe$pos <- as.numeric(str_replace(dataframe$id, "^S\\d[ABD]_", ""))
	dataframe$FDR <- p.adjust(dataframe$p, method = "fdr")
	dataframe$bon <- p.adjust(dataframe$p, method = "bonferroni")
	#filter out significant markers for future analysis
	####Should be bonferoni 0.05, testing
	significant_markers <- filter(dataframe, bon < .05)
	#significant_markers <- filter(dataframe, bon < .6)
	#plot out model
	par(cex = 1.5, cex.main = 1.5, pch = 1, lwd=3, mai = c(4,4,3,1))
	manhattan(dataframe, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = graph_name)
	abline(h = 5.6, col = "#659157")
	#Give a main dataframe all significant markers
	print(head(significant_markers))
	return(significant_markers)

}

sigMarkers <- data.frame()
for (i in c(2:length(GWAS_results[1,]))) {
	print(i)
	name <- substring(names(GWAS_results[i]), 3)
	print(name)
	png(width=2500, height=1500, pointsize = 36, filename = paste('output/plots/', name, '.png', sep=""))
	sig_m <- manhattan_plot(name, GWAS_results$id, GWAS_results[,i])
	dev.off()
	if (!(dim(sig_m)[1] == 0)) {
		sig_m$trait <- name
		sigMarkers <- rbind(sigMarkers, sig_m)
	}
}

###Export all significant markers along with their associated trait
write_csv(sigMarkers, "output/SunRILs_gwas_sig_markers_2022.csv")
