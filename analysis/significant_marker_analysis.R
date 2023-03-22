###Find populations with recombination at significant loci
##Created by Nicolas A. H. Lara
##Last edit: 2022-7-28

###SET WORKING DIRECTORY
#setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")

###LOAD PACKAGES
install.packages('snp_plotter')
library(tidyverse)
library(gaston)
source("analysis/GWAS_functions.R")

###READ IN AND PROCESS DATAFILES
sig_markers <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022.csv", sep=",")
Kin_sig <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022_Kinston.csv", sep=",")
Ral_sig <- read.delim("output/run_2022-08-18/SunRILs_gwas_sig_markers_2022_Raleigh.csv", sep=",")
genotype <- read.bed.matrix("data/SunRILs_filtered_2022")

#sig_markers <- Ral_sig
sig_markers <- subset(sig_markers, trait == 'Height' & LOG > 10)

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

fam_sig_id <- subset(fam_sig_id, MAF > .12)


SNP <- fam_sig_id[order(fam_sig_id$marker),]
SNP <- separate(data=SNP, col = marker, into = c('chr', 'id'), sep='\\_')
SNP$chr <- gsub('S','',SNP$chr)
SNP_convert <- data.frame()
chroms <- unique(SNP$chr)
bins <- c(1,1,1,1,2,3,3)
# for (i in 1:length(chroms)) {
	# chrom <- chroms[i]
	# max <- max(subset(genotype@snps, chr == chrom)$pos)
	# a <- subset(SNP, chr == chrom)
	# print(chrom)
	# convert <- round(as.numeric(unique(a$id))/max, 2)*100
	# print(convert)
	# if (bins[i] > 1) {
		# br <- cut(convert, breaks=bins[i])
		# print(br)
	# }
	# b <- data.frame(family = $family, SNP = )
	# print('')	
# }
one_bins <- c('short', 'short', 'short', 'long')
for (i in 1:length(chroms)) {
	chrom <- chroms[i]
	print(chrom)
	a <- subset(SNP, chr == chrom)
	a$id <- as.numeric(a$id)
	if (bins[i] == 1) {a$QTL <- one_bins[i]} #{print(unique(cut(as.numeric(a$id), breaks=bins[i])))}
	if (chrom == '4A') {a$QTL <- ifelse(a$id < 6.38e+08, "short","long")}
	if (chrom == '5A') {a$QTL <- ifelse(a$id >= 3.04e+08, ifelse(a$id >= 4.45e+08,'long','centromeric'),"short")}
	if (chrom == '6A') {a$QTL <- ifelse(a$id >= 2.14e+08, ifelse(a$id >= 3.67e+08,'long','centromeric'),"short")}
	#if (chrom == '5A') {if (a$id < 3.04e+08) {a$QTL <- "short"} else if (a$id >= 4.45e+08) {a$QTL <- "long"} else {a$QTL <- 'centromeric'}}
	#if (chrom == '6A') {if (a$id < 2.14e+08) {a$QTL <- "short"} else if (a$id >= 3.67e+08) {a$QTL <- "long"} else {a$QTL <- 'centromeric'}}
	SNP_convert <- rbind(SNP_convert, a)
}

plot_snp <- distinct(data.frame(family = SNP_convert$family, QTL = paste(SNP_convert$chr, SNP_convert$QTL, sep="_")))
subset_list <- c('Parents', 'UX1443', 'UX1444', 'UX1990', 'UX2002', 'UX2007', 'UX2008', 'UX2020', 'UX2028')

plsn <- matrix(nrow=length(unique(plot_snp$family)), ncol = length(unique(plot_snp$QTL)), dimnames=list(sort(unique(plot_snp$family), decreasing=TRUE), unique(plot_snp$QTL)))
for (i in rownames(plsn)) {
	for (j in colnames(plsn)) {
		if (dim(plot_snp[plot_snp$family == i & plot_snp$QTL == j,])[1] > 0) {plsn[i,j] = 1} else {plsn[i,j] = 0}
	}
}
plsn <- plsn[!(row.names(plsn) %in% subset_list),]

png(width=1500, height=1500, pointsize=30, filename = "output/sig_QTL_heatmap.png")
#par(mar=c(30,1,1,2))
heatmap(plsn, Rowv=NA, Colv=NA, scale='none', margins=c(8,4))
graphics.off()


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


###ID pops segregating for 6A PH QTL
SixA_qtl <- subset(sig_markers, chr == '6A' & trait == 'Height')
geno_sub <- select.snps(genotype, id %in% SixA_qtl$id)
subset_list <- c('Parents', 'UX1443', 'UX1444', 'UX1990', 'UX2002', 'UX2007', 'UX2008', 'UX2020', 'UX2028')
geno_sub <- select.inds(geno_sub, !(family %in% subset_list))

SixA_fam <- data.frame()
for (fam in unique(geno_sub@ped$family)) {
  print(fam)
  fam_group <- as.data.frame(as.matrix(select.inds(geno_sub, family == fam)))
  for (mark in colnames(fam_group)) {
    sub_geno <- select.inds(geno_sub, family==fam)
    sub_geno <- select.snps(sub_geno, id==mark)
    temp_df <- data.frame(family=fam, trait=traits, marker=mark, MAF=sub_geno@snps$maf, p=subset(sig_markers, (id == mark & trait == traits))$p)
    SixA_fam <- rbind(SixA_fam, temp_df)
  }
}

SixA_mat <- matrix(nrow = length(unique(SixA_fam$family)), ncol = length(unique(SixA_fam$marker)), dimnames=list(sort(unique(SixA_fam$family), decreasing=T), unique(SixA_fam$marker)))
for (i in rownames(SixA_mat)) {
  for (j in colnames(SixA_mat)) {
    SixA_mat[i,j] = subset(SixA_fam, family == i & marker == j)$MAF
  }
}
par(pin=c(5,5))
heatmap(SixA_mat, Rowv=NA, Colv=NA)#, scale='none', margins=c(8,4))
