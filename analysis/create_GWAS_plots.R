#Create plots from GWAS data
#Created by Nicolas A. H. Lara
#Last edit: 2022-7-26


#library(tidyverse)
library(gaston)
library(RAINBOWR)
library(RAINBOWR)

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")

#Read in GWAS file
#GWAS_results <- read.delim("output/data/SunRILs_gwas_2022.csv", sep=",")
#GWAS_kinston_heading <- read.delim("output/data/SunRILs_Kinston_heading_gwas_2022.csv", sep=",")
#GWAS_raleigh_heading <- read.delim("output/data/SunRILs_Raleigh_heading_gwas_2022.csv", sep=",")

manhattan_plot <- function(graph_name, SNP, p_value) {
	dataframe <- data.frame(id = SNP, p = p_value)
	dataframe <- dataframe[complete.cases(dataframe),]
	###Add new columns for statistical metrics, chromosomes, etc.
	dataframe$LOG <- -log10(dataframe$p)
	thresh_df <-data.frame(marker=substr(dataframe$id,2,15), chrom=substr(dataframe$id,2,3), pos=substr(dataframe$id,5,15),	LOG=dataframe$LOG)
	sig_thresh <- CalcThreshold(thresh_df,method=c("BH","Bonf"))
	dataframe$chr <- str_replace(str_replace(dataframe$id, "^S", ""), "_\\d*$", "")
	dataframe$pos <- as.numeric(str_replace(dataframe$id, "^S\\d[ABD]_", ""))
	dataframe$FDR <- p.adjust(dataframe$p, method = "fdr")
	dataframe$bon <- p.adjust(dataframe$p, method = "bonferroni")
	###filter out significant markers for future analysis
	significant_markers <- filter(dataframe, bon < .05)
	###plot out model
	par(cex = 1.5, cex.main = 1.5, pch = 1, lwd=3, mai = c(3,3,2,1))
	gaston::manhattan(dataframe, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = graph_name)
	abline(h = sig_thresh[[1]], col = "#d48b50")
	abline(h = sig_thresh[[2]], col = "#aed9a0")
	legend("topright", legend=c("Bonferroni Threshold","Benjamini-Hochberg Threshold"), fill=c("#aed9a0","#d48b50"),cex=0.8)
	###Give a main dataframe all significant markers
	#print(head(significant_markers))
	return(significant_markers)
}

