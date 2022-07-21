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
#GWAS_kinston_heading <- read.delim("output/data/SunRILs_Kinston_heading_gwas_2022.csv", sep=",")
#GWAS_raleigh_heading <- read.delim("output/data/SunRILs_Raleigh_heading_gwas_2022.csv", sep=",")

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
	par(cex = 1.5, cex.main = 1.5, pch = 1, lwd=3, mai = c(3,3,2,1))
	manhattan(dataframe, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = graph_name)
	abline(h = 5.6, col = "#659157")
	#Give a main dataframe all significant markers
	print(head(significant_markers))
	return(significant_markers)
}
#graph_name_list <- c("Mutation location", "Heading Date", "Height", "Awns", "Winter Dormancy Release", "Powdery Mildew", "Spikelet per Spike", "Infertile spikelet per Spike")
#graph_name_list <- c("Mutation location","Seed Weight", "Mean Seed Area", "Median Seed Area", "Seeds per Spikelet")


sigMarkers <- data.frame()
for (i in c(2:length(GWAS_results[1,]))) {
	print(i)
	#name <- substring(names(GWAS_results[i]), 3)
	name <- graph_name_list[i]
	print(name)
	png(width=2500, height=1500, pointsize = 15, filename = paste('output/plots/', name, '.png', sep=""))
	sig_m <- manhattan_plot(name, GWAS_results$id, GWAS_results[,i])
	dev.off()
	if (!(dim(sig_m)[1] == 0)) {
		sig_m$trait <- name
		sigMarkers <- rbind(sigMarkers, sig_m)
	}
}
###Export all significant markers along with their associated trait
write_csv(sigMarkers, "output/SunRILs_gwas_sig_markers_2022.csv")






#Calculating significance threshold
#According to Kaler and Purcell, 2019
###Y = a + bX 
#Y = sig. threshold for -log10 P-value
#a = intercept of regression
#b = slope of regression coefficient
#X = marker based heritability
Y = a + bX), where Y was the significant threshold
(−log 10 P-value), a was the intercept, and b was the slope
of the regression coefficient for the marker-based herit-
ability (X)


###Hard coded stuff for a few specific cases
#png(width=2500, height=1500, pointsize = 15, filename = paste('output/plots/', "Kinston Heading Date", '.png', sep=""))
#sig_m <- manhattan_plot("Kinston Heading Date", GWAS_kinston_heading$id, GWAS_kinston_heading$p_days_to_head)
#dev.off()

#png(width=2500, height=1500, pointsize = 15, filename = paste('output/plots/', "Raleigh Heading Date", '.png', sep=""))
#sig_m <- manhattan_plot("Raleigh Heading Date", GWAS_raleigh_heading$id, GWAS_raleigh_heading$p_days_to_head)
#dev.off()
#sig_m$trait <- "Kinston_heading"
#sigMarker_ral_kin_separate <- rbind(sigMarker_ral_kin_separate, sig_m)
#write_csv(sigMarker_ral_kin_separate, "output/sigMarker_ral_kin_separate_heading")
