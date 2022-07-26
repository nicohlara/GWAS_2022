###Perform GWAS on data from Spring 2022
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-7-26

#install.packages("RAINBOWR")
library(tidyverse)
library(gaston)
library(lmerTest)
library(RAINBOWR)
source("analysis/GWAS_functions.R")

#set working directory
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
###READ IN IMPUTED VCF FILE, CLEAN UP DATA
##IMPUTED USING BEAGLE5###only has some populations
#genotype <- read.vcf("data/used_vcf/SunRILs_2021_postimp_filt.vcf.gz", convert.chr = F)
##IMPUTED USING NOAH'S ALGORITHM###
genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader_filtered_imp.vcf.gz", convert.chr = F)

#unique(genotype@snps$chr)

###perform many cleaning steps on genotype dataset
genotype <- clean_genotype(genotype)

###SYNCING UP GENOTYPE AND PHENOTYPE DATA
genotype_sync <- select.inds(genotype, id %in% phenotype$Entry)
phenotype_sync <- filter(phenotype, Entry %in% genotype_sync@ped$id)

###SNB was only taken in Raleigh and can't be run though location based model
phenotype_run <- subset(phenotype_sync, select=-c(SNB))
#phenotype_run <- subset(phenotype_sync, Location == "Kinston", select=-c(Awns, SNB, Height, Powdery_mildew, ave_SpS, ave_infert, WeightOf1000Particles, SampleAreaAverage, WKAreaMedian, seeds_per_spikelet))

###RUN THE GWAS AND SAVE RESULTS TO A NEW DATATABLE
########can add specific marker locations into the model and remove them from the tested model list
########can also do iterative model and take out the most significant marker each time then repeat GWAS with that as fixed effect
plot_df <- GWAS_loop(genotype_sync, phenotype_run)

###Create a new directory for runs and plots if needed
dir <- paste("output/run_", Sys.Date(), sep="")
dir.create(dir)
write_csv(plot_df, paste(dir, "/SunRILs_gwas_2022.csv", sep=""))
#write_csv(plot_df, paste(dir, "/SunRILs_gwas_2022_Kinston.csv", sep=""))

###PLOT DATA###
GWAS_results <- plot_df

###Set names for graphs
graph_name_list <- substring(names(GWAS_results), 3)
#graph_name_list <- c("Mutation location", "Awns", "Winter Dormancy Release", "Powdery Mildew", "Height", "Spikelet per Spike", "Infertile spikelet per Spike", "Seed Weight", "Seed Area (mean)", "Seed Area (median)",  "Heading Date", "Seeds per Spikelet")
#graph_name_list <- c("Mutation location","Seed Weight", "Mean Seed Area", "Median Seed Area", "Seeds per Spikelet")

sigMarkers <- data.frame()
for (i in c(2:length(GWAS_results[1,]))) {
	name <- graph_name_list[i]
	png(width=2500, height=1500, pointsize = 15, filename = paste('output/plots/', name, '.png', sep=""))
	sig_m <- manhattan_plot(name, GWAS_results$id, GWAS_results[,i])
	dev.off()
	if (!(dim(sig_m)[1] == 0)) {
		sig_m$trait <- name
		sigMarkers <- rbind(sigMarkers, sig_m)
	}
}
###Export all significant markers along with their associated trait
write_csv(sigMarkers, paste(dir,"/SunRILs_gwas_sig_markers_2022.csv", sep=""))
#write_csv(sigMarkers, paste(dir,"/SunRILs_gwas_sig_markers_2022_Kinston.csv", sep=""))