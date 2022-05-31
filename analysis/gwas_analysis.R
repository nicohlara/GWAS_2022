#Perform GWAS on phenotypic data from Spring 2022 field data
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Last edit: 2022-5-26

#Load in packages
library(tidyverse)
library(gaston)
library(rrBLUP)

#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")

###READ IN PHENOTYPE FILE, CLEAN COLUMN NAMES, COMBINE LOCATIONS
K22_pheno <- read.delim("data/Kin22-SunRils-T1-T30 Final.xlsx - 2022-05-12-11-55-18_Kin22-SunRi.csv", sep=",")
R22_pheno <- read.delim("data/R22-SunRil T1-30 Final.xlsx - Sheet1.csv", sep=",")
#clean up dataframes, rename columns
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
K22_pheno <- header.true(K22_pheno)
names(K22_pheno)[14] <- "N_rust"
#names(K22_pheno)[names(K22_pheno) == "Tray"] <- "Kin_tray"
names(K22_pheno)[15] <- "G_rust"
names(K22_pheno)[names(K22_pheno) == "WDR 3/7"] <- "WDR"
#names(K22_pheno)[names(K22_pheno) == "flowering"] <- "Kin_flowering"
names(K22_pheno)[names(K22_pheno) == "Powdery Mildew"] <- "Powdery_mildew"
#names(R22_pheno)[names(R22_pheno) == "WDR"] <- "Ral_WDR"
#names(R22_pheno)[names(R22_pheno) == "flowering"] <- "Ral_flowering"
names(R22_pheno)[names(R22_pheno) == "Powdery.Mildew"] <- "Powdery_mildew"
#names(R22_pheno)[names(R22_pheno) == "Tray"] <- "Ral_tray"
names(R22_pheno)[names(R22_pheno) == "Plot_Id"] <- "Plot_ID"
names(R22_pheno)[names(R22_pheno) == "Cross_Num"] <- "Cross_ID"
names(R22_pheno)[names(R22_pheno) == "Cross_id"] <- "Cross"
#names(R22_pheno)[names(R22_pheno) == "SNB"] <- "Ral_SNB"
R22_pheno["Location"] <- "Raleigh"
K22_pheno["Location"] <- "Kinston"
R22_pheno <- R22_pheno[!(R22_pheno$Entry=="Barley"),]
K22_pheno <- K22_pheno[!(K22_pheno$Entry=="Barley"),]
K22_pheno$SNB <- ""
#Combine phenotype files
K22_pheno_subset <- K22_pheno[c("Cross_ID", "Entry", "WDR", "flowering", "Powdery_mildew", "SNB", "Location")]
R22_pheno_subset <- R22_pheno[c("Cross_ID", "Entry", "WDR", "flowering", "Powdery_mildew", "SNB", "Location")]
t_ph <- rbind(K22_pheno_subset, R22_pheno_subset)
write.table(t_ph, file='output/all_phenotypes.tsv', quote=FALSE, sep='\t', row.names = FALSE)

###READ IN IMPUTED VCF FILE, CLEAN UP DATA
sunVCF <- read.vcf("data/SunRILs_2021_postimp_filt.vcf.gz", convert.chr = F)
#filter out parents and thin on LD
sunVCF <- select.inds(sunVCF, grepl("^UX", id))
sunVCF <- LD.thin(sunVCF, threshold = .8)
sunM <- as.matrix(sunVCF)
#reformat ids to get rid of extra text
sunVCF@ped$id <- gsub("-NWG", "", sunVCF@ped$id)
sunVCF@ped$id <- gsub("-A+", "", sunVCF@ped$id)
sunVCF@ped$id <- gsub("-NEG", "", sunVCF@ped$id)
#filter for resequenced and duplicated lines
sunVCF <- select.inds(sunVCF, !grepl("-A-", id))
sunVCF <- select.inds(sunVCF, !duplicated(id))



###CHECK DATA TO ENSURE RELATIVE NORMALITY
#Phenotypic data
hist(as.numeric(t_ph$WDR), main="Winter dormancy release 2022")
#have to convert flowering time to days
days <- c()
#for (i in length(t_ph$flowering)) {
for (i in c(1:length(t_ph$flowering))) {
	if (t_ph$Location[i] == "Kinston") {
		a <- difftime(as.Date(t_ph$flowering[i], format="%m/%d/%Y"),ISOdate(2021,10,28))
	} else {
		a <- difftime(as.Date(t_ph$flowering[i], format="%m/%d/%Y"),ISOdate(2021,11,4))
	}
	days <- c(days, as.numeric(a) -.5)
}
t_ph$days_to_head <- days
hist(as.numeric(days), breaks = 40, main = "Flowering date 2022")
hist(as.numeric(t_ph$Powdery_mildew), main="Powdery mildew 2022")
hist(as.numeric(t_ph$SNB), main="Septoria nodorum blotch 2022")
#location testing
sum(t_ph$Location=="Kinston")
sum(t_ph$Location=="Raleigh")
#Perform ANOVA to look at variation in heading date
anova(lm(data = t_ph, days_to_head ~ Location + Cross_ID + Location:Cross_ID))

#Genotypic data
#check family structure using heatmap
sunMTest <- sunM[sample(c(1:nrow(sunM)), 1000), sample(c(1:ncol(sunM)), 5000)]
sunGRMTest <- A.mat(sunMTest)
heatmap(sunGRMTest, symm = T)


###CLEAN UP PHENOTYPIC DATA FOR COMBINING WITH GENOTYPE
#check for duplication and lack of replication
genoCounts <- group_by(t_ph, Entry) %>% count()
parLines <- filter(genoCounts, n > 10)$Entry
single_crossLines <- filter(genoCounts, 1 == n)$Entry
crossLines <- filter(genoCounts, n == 2 )$Entry
multi_crossLines <- filter(genoCounts, 2 < n & 10 > n )$Entry 

#since we have 23 lines with more than 1 replicate at each location, we average them
t_ph_group <- group_by(t_ph, Location, Cross_ID, Entry) %>% summarize(days_to_head = mean(days_to_head))
#filter out any unreplicated lines
t_ph_group <- filter(t_ph_group, !(Entry %in% single_crossLines))
t_ph_group <- as.data.frame(t_ph_group)
#rownames(t_ph_group) <- t_ph_group$Entry


###SYNCING UP GENOTYPE AND PHENOTYPE DATA
#process to sync up phenotype and genotype files
sunVCF_sync <- select.inds(sunVCF, id %in% t_ph_group$Entry)
t_ph_group <- filter(t_ph_group, Entry %in% sunVCF_sync@ped$id)



###Incredibly janky filtering for matrix work, makes model run but results aren't great looking
matrix_df <- data.frame()
#matrix by individual IN family, works better than by family
for (i in unique(t_ph_group$Cross_ID)) {
	a <- filter(t_ph_group, Cross_ID == i)
	b <- sample(unique(a$Entry), length(unique(a$Entry))/2)
	df <- filter(a, (Entry %in% b & Location == "Kinston") | (!(Entry %in% b) & !(Location == "Kinston")))
	matrix_df <- rbind(matrix_df, df)
}
#matrix by family, 
#matrix_df <- data.frame()
#rand_cross_id <- sample(unique(t_ph_group$Cross_ID), length(unique(t_ph_group$Cross_ID))/2)
#matrix_df <- filter(t_ph_group, (Cross_ID %in% rand_cross_id & Location == "Kinston") | (!(Cross_ID %in% rand_cross_id) & (Location == "Raleigh")))
#write_csv(matrix_df, "output/matrix_df.csv")
#rownames(matrix_df) <- matrix_df$Entry


#hist(t_ph_group$days_to_head)
###MAKE AND SHOW MODELS
#make covariate model
#familyLocCov <- model.matrix(~ Location + Entry, t_ph_group)
#familyLocCov <- model.matrix(~gsub("-\\d*$", "", t_ph_group$Entry) + t_ph_group$Location)
familyLocCov <- model.matrix(~ Location + Entry, matrix_df)
familyLocCov <- model.matrix(~gsub("-\\d*$", "", matrix_df$Entry) + matrix_df$Location)

######NEED TO SYNC UP LENGTH OF SUNVCF AND T_PHGN
testGWAS <- association.test(sunVCF_sync, matrix_df$days_to_head, X = familyLocCov)
#testGWAS <- association.test(sunVCF_sync, t_ph_group$days_to_head, X = familyLocCov)
testGWAS$LOG <- -log10(testGWAS$p)
testGWAS$chr <- str_replace(str_replace(testGWAS$id, "^S", ""), "_\\d*$", "")
testGWAS$pos <- as.numeric(str_replace(testGWAS$id, "^S\\d[ABD]_", ""))
testGWAS$FDR <- p.adjust(testGWAS$p, method = "fdr")
testGWAS$bon <- p.adjust(testGWAS$p, method = "bonferroni")

#plot out model
manhattan(testGWAS, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = "SunRILs Heading Date GWAS")
abline(h = 5.6, col = "#659157")

#Export significant markers for heading date
sigMarkers <- filter(testGWAS, bon < .05)
write_csv(sigMarkers, "output/SunRILs_heading_days_gwasOutput_2022.csv")















