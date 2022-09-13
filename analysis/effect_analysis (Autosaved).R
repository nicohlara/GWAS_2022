#####Using the qtl2 package, create a cross2 object and perform analysis on it
#Created by Nicolas A. H. Lara
#Last edit: 2022-9-13

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
library(tidyverse)
library(gaston)
library(qtl2)

phenotype <- read.delim("output/data/2022_phenotype.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

###creating cross2 files
#filter to significant markers
sig_markers <- read.csv("output/SunRILs_sig_mark_2022.csv", sep=",")
#filter by height
sig_markers <- subset(sig_markers, trait=="Height")# & p < 1e-10)
geno_subset <- select.snps(genotype, id %in% sig_markers$marker)
#geno_subset <- select.inds(geno_subset, family %in% c("UX2029", "UX2000", "UX1995"))
geno <- as.matrix(geno_subset)
geno[geno==0] <- "AA"
geno[geno==1] <- "AB"
geno[geno==2] <- "BB"
write.csv(geno, "output/data/SunCross/SunRILs_geno.csv")
gmap <- subset(geno_subset@snps, select=c("id", "chr", "pos"))
names(gmap)[names(gmap)=="id"] <- "marker"
write.csv(gmap, "output/data/SunCross/SunRILs_gmap.csv", row.names=FALSE)
pheno <- subset(phenotype, select=c("Entry", "Location", "Height"))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995"))
pheno_K <- subset(pheno, Location == "Kinston", select=c("Entry", "Height"))
pheno_K <- rename(pheno_K, Kinson = Height)
pheno_R <- subset(pheno, Location == "Raleigh", select=c("Entry", "Height"))
pheno_R <- rename(pheno_R, Raleigh = Height)
pheno <- merge(pheno_R, pheno_K, by="Entry")
write.csv(pheno, "output/data/SunCross/SunRILs_pheno.csv", row.names=FALSE)
#phenocovar <- data.frame(pheno=c("Kinston", "Raleigh"), location = c("Kinston", "Raleigh"))
#write.csv(phenocovar, "/Users/nico/Desktop/cross2_file_test/SunRILs_phenocovar.csv", row.names=FALSE)
covar <- distinct(subset(phenotype, select=c("Entry", "Cross_ID")))#, Cross_ID %in% c("UX2029", "UX2000", "UX1995")))
write.csv(covar, "output/data/SunCross/SunRILs_covar.csv", row.names=FALSE)
SunCross <- read_cross2("output/data/SunCross/SunRILs.yaml")

# ###Performing analysis
# pr <- calc_genoprob(cross=SunCross, map=NULL)
# c2eff_me <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"])#iron$pheno[,"liver"])
# plot_coef(c2eff_me, map=SunCross$gmap, columns=1:2)
# #c2effB <- scan1coef(pr[,"6A"], SunCross$pheno[,"Raleigh"], contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
# # include BLUPs
# c2blup <- scan1blup(pr[,"6A"], SunCross$pheno[,"Raleigh"])
# plot_coef(c2eff_me, SunCross$gmap["6A"], columns=1:3)
# plot(c2blup, SunCross$gmap["6A"], columns=1:3, add=TRUE, lty=2, legend = "topright")


# dir.create("output/effects")
SunProb <- calc_genoprob(cross=SunCross, map=NULL)
SunEffect <- data.frame()
for (i in chr_names(SunCross)) {
	#do analysis
	SunEff <- scan1coef(SunProb[,i], SunCross$pheno[,"Raleigh"])
	SunBlup <- scan1blup(SunProb[,i], SunCross$pheno[,"Raleigh"])
	#plot out and save
	png(width=2500, height=1500, pointsize = 15, filename = paste("output/effects/", i, '_effect.png', sep=""))
	plot_coef(SunEff, SunCross$gmap[i], columns=1:2, main=paste(i, "Effect Score", sep=""))
	plot(SunBlup, SunCross$gmap[i], columns=1:2, add=TRUE, lty=2, legend = "topright")
	dev.off()
	#add to main file
	SunEffect <- rbind(SunEffect, SunEff)
}
write.csv(SunEffect, "output/effects/SunRILs_height_effects.csv")


