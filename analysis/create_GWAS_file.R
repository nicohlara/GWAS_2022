###Perform GWAS on data from Spring 2022
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-7-11

library(tidyverse)
library(gaston)
library(lmerTest)

#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
###READ IN IMPUTED VCF FILE, CLEAN UP DATA
##USED:
#---only has some populations
#genotype <- read.vcf("data/used_vcf/SunRILs_2021_postimp_filt.vcf.gz", convert.chr = F)
#Only has first chromosome?
#genotype <- read.vcf("data/used_vcf/SunRILs_combGenos_fakeHeader_imp.vcf.gz", convert.chr = F)

#To try:
genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader_filtered_imp.vcf.gz", convert.chr = F)
#genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader_filtered.vcf.gz", convert.chr = F)
#genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader.vcf.gz", convert.chr = F)

##Combine by hand:
#genotype_file_list <- list.files(path = "data/population_VCF_files/", pattern = "*.vcf.gz")
#genotype_files <- lapply(genotype_file_list, read.vcf)

unique(genotype@snps$chr)

#filter out parents and thin on LD
genotype <- select.inds(genotype, grepl("^UX", id))
genotype <- LD.thin(genotype, threshold = .8)
genotype_matrix <- as.matrix(genotype)
#reformat ids to get rid of extra text
genotype@ped$id <- gsub("-NWG", "", genotype@ped$id)
genotype@ped$id <- gsub("-A+", "", genotype@ped$id)
genotype@ped$id <- gsub("-NEG", "", genotype@ped$id)
#filter for resequenced and duplicated lines
genotype <- select.inds(genotype, !grepl("-A-", id))
genotype <- select.inds(genotype, !duplicated(id))

###SYNCING UP GENOTYPE AND PHENOTYPE DATA
genotype_sync <- select.inds(genotype, id %in% phenotype$Entry)
phenotype_sync <- filter(phenotype, Entry %in% genotype_sync@ped$id)

###CREATE BASIC MODEL FOR MANHATTAN PLOT
plot_df <- data.frame()
len <- length(genotype_sync@snps$id)
for (i in c(1:len)) {
	print(paste(round((i/len)*100, digits=2), "%", sep=""))
	a <- as.matrix(genotype_sync[,i])
	b <- data.frame(Entry = rownames(a), Marker = as.vector(a[,1]))
	c <- colnames(a)
	if (length(unique(b$Marker)) > 1) {
		d <- merge(phenotype_sync,b, by="Entry")
		p_vals <- data.frame(id = c)
		for (j in c(4:length(phenotype_sync[1,]))) {
			###STATISTICAL TESTS
			#VERY BASIC LINEAR TEST
			#mm <- lm(data = d, d[,j] ~ Marker)
			#p <- summary(mm)$coefficients["Marker",4]
			#TREATING LOCATION AND FAMILY RANDOM EFFECTS, takes a while to run
			mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Location) + (1|Cross_ID)))
			p <- anova(mm)["Marker",6]
			p_vals[paste("p_", names(phenotype_sync[j]), sep="")] <- p
		}
		plot_df <- rbind(plot_df, p_vals)
		#plot_df <- rbind(plot_df, data.frame(id = c, p = p_val))
	}
}
write_csv(plot_df, "output/data/SunRILs_gwas_2022.csv")
