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
sig_markers <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022.csv", sep=",")
Kin_sig <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022_Kinston.csv", sep=",")
Ral_sig <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022_Raleigh.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

sig_markers <- Ral_sig

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
			sub_geno <- select.inds(geno_sub, family==fam)
			sub_geno <- select.snps(sub_geno, id==mark)
			
			#variant <- unique(fam_group[[mark]])
			#temp_df <- data.frame(family = fam, trait = traits, marker = mark, p = subset(sig_markers, id == mark)$p, variants = paste(variant, collapse=", "))
			#if (length(variant) == 1) {temp_df["hom_het_state"] <- "homozygous"} else {temp_df["hom_het_state"] <- "heterozygous"}
			
			temp_df <- data.frame(family=fam, trait=traits, marker=mark, MAF=sub_geno@snps$maf, p=subset(sig_markers, (id == mark & trait == traits))$p)
			fam_sig_id <- rbind(fam_sig_id, temp_df)
		}
	}
}
fam_sig_id <- subset(fam_sig_id, MAF > .25)

write_csv(fam_sig_id, "output/SunRILs_significant_id_2022.csv")
# het_sig_id <-  subset(fam_sig_id, hom_het_state == "heterozygous")
# write_csv(het_sig_id, "output/SunRILs_segregating_sign_markers.csv")
# het_sig_id <- read.delim("output/SunRILs_segregating_sign_markers.csv", sep=",")
# maf_merge <- data.frame()
# for (i in c(1:length(het_sig_id$marker))) {
	# sub_geno <- select.inds(genotype, family==het_sig_id[i,1])
	# sub_geno <- select.snps(sub_geno, id==het_sig_id[i,3])
	# temp_df <- het_sig_id[i,]
	# temp_df$MAF <- sub_geno@snps$maf
	# maf_merge <- rbind(maf_merge, temp_df)
# }
# write_csv(maf_merge, "output/SunRILs_MAF_sig_id_2022.csv")

subset_markers <- subset(fam_sig_id, MAF > .25)


significant_markers_all <- read.delim("output/SunRILs_family_significant_id_2022.csv", sep=",")
significant_markers_Kinston <- read.delim("output/SunRILs_family_significant_id_2022_Kinston.csv", sep=",")
significant_markers_Raleigh <- read.delim("output/SunRILs_family_significant_id_2022_Raleigh.csv", sep=",")


significant_markers_all$Location <- "All"
names(significant_markers_all)[names(significant_markers_all) == "MAF....sub_geno.snps.maf"] <- "MAF"
significant_markers_Kinston$Location <- "Kinston"
significant_markers_Raleigh$Location <- "Raleigh"
significant_markers <- rbind(significant_markers_all, significant_markers_Kinston, significant_markers_Raleigh)
significant_markers <- subset(significant_markers, MAF > .15)
significant_markers <- arrange(significant_markers, trait, marker)
write_csv(significant_markers, "output/SunRILs_sig_markers_2022.csv")