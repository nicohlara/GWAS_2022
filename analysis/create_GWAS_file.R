###Perform GWAS on data from Spring 2022
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created and modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-9-1

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

#install.packages("pegas")
library(tidyverse)
library(gaston)
library(lmerTest)
library(RAINBOWR)
source("analysis/GWAS_functions.R")

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

###SYNCING UP PHENOTYPE DATA

phenotype_sync <- filter(phenotype, Entry %in% genotype@ped$id)

###SNB was only taken in Raleigh and can't be run though location based model
###Rust was only taken in Kinston and can't be run though location based model
#phenotype_run <- subset(phenotype_sync, select=-c(SNB, Rust))
#phenotype_run <- subset(phenotype_sync, Location == "Raleigh", select=-c(Awns, Rust, Powdery_mildew))

###RUN THE GWAS AND SAVE RESULTS TO A NEW DATATABLE
########can add specific marker locations into the model and remove them from the tested model list
########can also do iterative model and take out the most significant marker each time then repeat GWAS with that as fixed effect

plot_df <- GWAS_loop(genotype, phenotype_run)

###Create a new directory for runs and plots if needed
dir <- paste("output/run_", Sys.Date(), sep="")
dir.create(dir)
#write_csv(plot_df, paste(dir, "/SunRILs_gwas_2022.csv", sep=""))
write_csv(plot_df, paste(dir, "/SunRILs_gwas_2022_Raleigh.csv", sep=""))

###PLOT DATA###
GWAS_results <- plot_df
#GWAS_results <- read.delim("output/run_2022-07-26/SunRILs_gwas_2022.csv", sep=",")

###Set names for graphs
graph_name_list <- substring(names(GWAS_results), 3)
#graph_name_list <- c("Mutation location", "Awns", "Winter Dormancy Release", "Powdery Mildew", "Height", "Spikelet per Spike", "Infertile spikelet per Spike", "Seed Weight", "Seed Area (mean)", "Seed Area (median)",  "Heading Date", "Seeds per Spikelet")
#graph_name_list <- c("Mutation location","Seed Weight", "Mean Seed Area", "Median Seed Area", "Seeds per Spikelet")
#dir.create(paste(dir, "/plots", sep=""))
dir.create(paste(dir, "/Raleigh_plots", sep=""))
sigMarkers <- data.frame()
for (i in c(2:length(GWAS_results[1,]))) {
	name <- graph_name_list[i]
	png(width=2500, height=1500, pointsize = 15, filename = paste(dir, "/Raleigh_plots/", name, '.png', sep=""))
	sig_m <- manhattan_plot(name, GWAS_results$id, GWAS_results[,i])
	dev.off()
	if (!(dim(sig_m)[1] == 0)) {
		sig_m$trait <- name
		sigMarkers <- rbind(sigMarkers, sig_m)
	}
}
###Export all significant markers along with their associated trait
#write_csv(sigMarkers, paste(dir,"/SunRILs_gwas_sig_markers_2022.csv", sep=""))
write_csv(sigMarkers, paste(dir,"/SunRILs_gwas_sig_markers_2022_Raleigh.csv", sep=""))