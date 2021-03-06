###Find populations with recombination at significant loci
##Created by Nicolas A. H. Lara
##Last edit: 2022-7-28

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(tidyverse)
library(gaston)
source("analysis/GWAS_functions.R")

###READ IN AND PROCESS DATAFILES
sig_markers <- read.delim("output/run_2022-07-26/SunRILs_gwas_sig_markers_2022.csv", sep=",")
genotype <- read.vcf("data/SunRILs_combGenos_fakeHeader_filtered_imp.vcf.gz", convert.chr = F)

##Clean up genotype
genotype <- clean_genotype(genotype)
##Add family marker for selection
#genotype@ped$family 
a <- data.frame(ID = genotype@ped$famid) %>% separate(ID, c("family"), extra = "drop")
genotype@ped$family <- a[,1]

###LOOPING THROUGH SIGNIFICANT MARKERS AND GETTING HET STATE FOR EACH FAMILY
fam_sig_id <- data.frame()
for (traits in unique(sig_markers$trait)) {
	print(traits)
	##subset trait
	markers <- subset(sig_markers, trait == traits)
	print(length(markers$trait))
	##select snps
	geno_sub <- select.snps(genotype, id %in% markers$id)
	for (fam in unique(geno_sub@ped$family)) {
		print(fam)
		fam_group <- as.data.frame(as.matrix(select.inds(geno_sub, family == fam)))
		for (mark in colnames(fam_group)) {
			variant <- unique(fam_group[[mark]])
			temp_df <- data.frame(family = fam, trait = traits, marker = mark, p = subset(sig_markers, id == mark)$p, variants = paste(variant, collapse=", "))
			if (length(variant) == 1) {temp_df["hom_het_state"] <- "homozygous"} else {temp_df["hom_het_state"] <- "heterozygous"}
			fam_sig_id <- rbind(fam_sig_id, temp_df)
		}
	}
}

write_csv(fam_sig_id, "output/SunRILs_family_significant_id_2022.csv")
het_sig_id <-  subset(fam_sig_id, hom_het_state == "heterozygous")
write_csv(het_sig_id, "output/SunRILs_segregating_sign_markers.csv")