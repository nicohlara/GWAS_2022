#Create plots from GWAS data
#Created by Nicolas A. H. Lara
#Last edit: 2022-7-22


#library(tidyverse)
library(gaston)
library(RAINBOWR)

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#Read in GWAS file
#GWAS_results <- read.delim("output/data/SunRILs_gwas_2022.csv", sep=",")
#GWAS_kinston_heading <- read.delim("output/data/SunRILs_Kinston_heading_gwas_2022.csv", sep=",")
#GWAS_raleigh_heading <- read.delim("output/data/SunRILs_Raleigh_heading_gwas_2022.csv", sep=",")

clean_genotype <- function(genotype) {
	#filter out parents and thin on LD
	genotype <- select.inds(genotype, grepl("^UX", id))
	genotype <- LD.thin(genotype, threshold = .8, max.dist = 350e6) #made for humans, Noah suggested changing to ~350MB 
	genotype_matrix <- as.matrix(genotype)
	#reformat ids to get rid of extra text
	genotype@ped$id <- gsub("-NWG", "", genotype@ped$id)
	genotype@ped$id <- gsub("-A+", "", genotype@ped$id)
	genotype@ped$id <- gsub("-NEG", "", genotype@ped$id)
	#filter for resequenced and duplicated lines
	genotype <- select.inds(genotype, !grepl("-A-", id))
	genotype <- select.inds(genotype, !duplicated(id))
	return(genotype)
}

GWAS_loop <- function(genotype, phenotype) {
	plot_df <- data.frame()
	len <- length(genotype@snps$id)
	count <- 0
	print(Sys.time())
	for (i in c(1:len)) {
		if (round((i/len)*100) > count) {
			count <- round((i/len)*100)
			print(paste(count, "%", sep=""))
	}	
	a <- as.matrix(genotype[,i])
	b <- data.frame(Entry = rownames(a), Marker = as.vector(a[,1]))
	c <- colnames(a)
	if (length(unique(b$Marker)) > 1) {
		d <- merge(phenotype,b, by="Entry")
		p_vals <- data.frame(id = c)
		for (j in c(4:length(phenotype[1,]))) {
			if (length(unique(d$Location)) > 1) {
				###TREATING LOCATION AND FAMILY RANDOM EFFECTS
				mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Location) + (1|Cross_ID)))
			} else {
				###Running without location
				mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Cross_ID)))
			}
			p <- anova(mm)["Marker",6]
			p_vals[paste("p_", names(phenotype[j]), sep="")] <- p
		}
		plot_df <- rbind(plot_df, p_vals)
		#plot_df <- rbind(plot_df, data.frame(id = c, p = p_val))
	}
	return(plot_df)
}

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
