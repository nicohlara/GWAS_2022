#Perform GWAS on phenotypic data from Spring 2022 field data
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-7-11

#install.packages("tidyverse")
#install.packages("gaston")
##install.packages("rrBLUP")
#install.packages("lmerTest")

#Load in packages
library(tidyverse)
library(gaston)
#library(rrBLUP)
#install.packages("lmerTest")
#library(lme4)
library(lmerTest)

#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
###READ IN IMPUTED VCF FILE, CLEAN UP DATA
#sunVCF <- read.vcf("data/SunRILs_2021_postimp_filt.vcf.gz", convert.chr = F)
sunVCF <- read.vcf("data/SunRILs_combGenos_fakeHeader_imp.vcf.gz")

###READ IN PHENOTYPE FILE, CLEAN COLUMN NAMES, COMBINE LOCATIONS
#K22_pheno <- read.delim("data/Kin22-SunRils-T1-T30 Final.xlsx - 2022-05-12-11-55-18_Kin22-SunRi.csv", sep=",")
K22_pheno <- read.delim("data/Kin22-SunRils-T1-T30.csv", sep=",")
R22_pheno <- read.delim("data/R22-SunRil-T1-30.csv", sep=",")
K22_SpS <- read.delim("data/SpS_data.xlsx - Kinston.csv", sep=",")
R22_SpS <- read.delim("data/SpS_data.xlsx - Raleigh.csv", sep=",")
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
#R22_pheno["Height"] <- 0
R22_pheno <- R22_pheno[!(R22_pheno$Entry=="Barley"),]
K22_pheno <- K22_pheno[!(K22_pheno$Entry=="Barley"),]
K22_pheno$SNB <- ""
#CLEAN DATA/MISLABELLED LINES
K22_pheno <- mutate(K22_pheno, Entry = if_else(Entry == "UX1994-6.4", "UX1994-64", Entry))
K22_pheno <- mutate(K22_pheno, Entry = if_else(Entry == "UX2029-4.6", "UX2029-46", Entry))
#which(R22_pheno$Entry == "UX2029-4.6")
#Process Spikelet per spike (SpS) data files
K22_SpS$ave_SpS <- apply(K22_SpS[,grepl("SpS", names(K22_SpS))], 1, mean, na.rm = TRUE)
K22_SpS$ave_infert <- apply(K22_SpS[,grepl("Inf", names(K22_SpS))], 1, mean, na.rm = TRUE)
R22_SpS$ave_SpS <- apply(R22_SpS[,grepl("SpS", names(R22_SpS))], 1, mean, na.rm = TRUE)
R22_SpS$ave_infert <- apply(R22_SpS[,grepl("Inf", names(R22_SpS))], 1, mean, na.rm = TRUE)
#Add SpS averages to phenotype files
K22_pheno <- merge(K22_pheno, K22_SpS[,c("Tray", "ave_SpS", "ave_infert")], by= "Tray", all.x = TRUE)
R22_pheno <- merge(R22_pheno, R22_SpS[,c("Tray", "ave_SpS", "ave_infert")], by= "Tray", all.x = TRUE)
#Combine phenotype files
K22_pheno_subset <- K22_pheno[c("Cross_ID", "Entry", "Awns", "WDR", "flowering", "Powdery_mildew", "SNB", "Location", "Height", "ave_SpS", "ave_infert")]
R22_pheno_subset <- R22_pheno[c("Cross_ID", "Entry", "Awns", "WDR", "flowering", "Powdery_mildew", "SNB", "Location", "Height",  "ave_SpS", "ave_infert")]
t_ph <- rbind(K22_pheno_subset, R22_pheno_subset)
t_ph$Height <- as.numeric(t_ph$Height)
#ADD COLUMN CONVERTING FLOWERING TIME FROM DATE TO DAYS
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
write.table(t_ph, file='output/all_phenotypes.tsv', quote=FALSE, sep='\t', row.names = FALSE)
#Change awn from + to 1,0 for coding purposes
awn_vector <- ifelse(t_ph$Awns=="+", 1,0)
t_ph$Awns <- awn_vector



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

#check coverage of genotype to phenotype
#vcf_id <- sunVCF@ped$famid
#t_ph$Entry
#intersect(sunVCF@ped$famid,t_ph$Entry)

#check preimpute
#sunVCF2 <- read.vcf("data/initialRun/filt_and_imp_VCF/SunRILs_2021_preimp_filt_imp.vcf.gz", convert.chr = F)
#sunVCF2 <- select.inds(sunVCF2, grepl("^UX", id))

###CHECK DATA TO ENSURE RELATIVE NORMALITY
#Phenotypic data
hist(as.numeric(t_ph$WDR), main="Winter dormancy release 2022")
hist(as.numeric(days), breaks = 40, main = "Flowering date 2022")
hist(as.numeric(t_ph$Powdery_mildew), main="Powdery mildew 2022")
hist(as.numeric(t_ph$SNB), main="Septoria nodorum blotch 2022")
hist(t_ph$Height, main="Height")
#location testing
sum(t_ph$Location=="Kinston")
sum(t_ph$Location=="Raleigh")
#Perform ANOVA to look at variation in heading date
anova(lm(data = t_ph, days_to_head ~ Location + Cross_ID + Location:Cross_ID))
anova(lm(data = t_ph, Height ~ Location + Cross_ID + Location:Cross_ID))
#Genotypic data
#check family structure using heatmap
sunMTest <- sunM[sample(c(1:nrow(sunM)), 1000), sample(c(1:ncol(sunM)), 5000)]
sunGRMTest <- A.mat(sunMTest)
heatmap(sunGRMTest, symm = T)

###CLEAN UP PHENOTYPIC DATA FOR COMBINING WITH GENOTYPE
#check for duplication and lack of replication
genoCounts <- group_by(t_ph, Entry) %>% count() %>% ungroup()
parLines <- filter(genoCounts, n > 10)$Entry
single_crossLines <- filter(genoCounts, 1 == n)$Entry
crossLines <- filter(genoCounts, n == 2 )$Entry
multi_crossLines <- filter(genoCounts, 2 < n & 10 > n )$Entry 

#since we have 23 lines with more than 1 replicate at each location, we average them
t_ph_group <- t_ph %>% group_by(Location, Entry, Cross_ID) %>% summarize(days_to_head = mean(days_to_head), Height = mean(Height), Awns = mean(Awns), WDR = mean(as.numeric(WDR)), Powdery_mildew = mean(as.numeric(Powdery_mildew)), ave_SpS = mean(ave_SpS), ave_infert = mean(ave_infert)) %>% ungroup()
# %>% mutate(loc_fam = paste(Location, Cross_ID, sep="_"))
#filter out any unreplicated lines
t_ph_group <- filter(t_ph_group, !(Entry %in% single_crossLines))
t_ph_group <- as.data.frame(t_ph_group)
#rownames(t_ph_group) <- t_ph_group$Entry

###SYNCING UP GENOTYPE AND PHENOTYPE DATA
#process to sync up phenotype and genotype files
sunVCF_sync <- select.inds(sunVCF, id %in% t_ph_group$Entry)
t_ph_group <- filter(t_ph_group, Entry %in% sunVCF_sync@ped$id)

###CREATE BASIC MODEL FOR MANHATTAN PLOT
#plot_df <- data.frame(id = character(), p = numeric())
plot_df <- data.frame()
len <- length(sunVCF_sync@snps$id)
for (i in c(1:len)) {
	print(paste(round((i/len)*100, digits=2), "%", sep=""))
	a <- as.matrix(sunVCF_sync[,i])
	b <- data.frame(Entry = rownames(a), Marker = as.vector(a[,1]))
	c <- colnames(a)
	if (length(unique(b$Marker)) > 1) {
		d <- merge(t_ph_group,b, by="Entry")
		p_vals <- data.frame(id = c)
		for (j in c(4:length(t_ph_group[1,]))) {
			###STATISTICAL TESTS
			#VERY BASIC LINEAR TEST
			#mm <- lm(data = d, d[,j] ~ Marker)
			#p <- summary(mm)$coefficients["Marker",4]
			#TREATING LOCATION AND FAMILY RANDOM EFFECTS, takes a while to run
			mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Location) + (1|Cross_ID)))
			p <- anova(mm)["Marker",6]
			p_vals[paste("p_", names(t_ph_group[j]), sep="")] <- p
		}
		plot_df <- rbind(plot_df, p_vals)
		#plot_df <- rbind(plot_df, data.frame(id = c, p = p_val))
	}
}
write_csv(plot_df, "output/SunRILs_gwas_2022.csv")

manhattan_plot <- function(graph_name, SNP, p_value) {
	dataframe <- data.frame(id = SNP, p = p_value)
	dataframe <- dataframe[complete.cases(dataframe),]
	#Add new columns for staticstical metrics, chromosomes, etc.
	dataframe$LOG <- -log10(dataframe$p)
	dataframe$chr <- str_replace(str_replace(dataframe$id, "^S", ""), "_\\d*$", "")
	dataframe$pos <- as.numeric(str_replace(dataframe$id, "^S\\d[ABD]_", ""))
	dataframe$FDR <- p.adjust(dataframe$p, method = "fdr")
	dataframe$bon <- p.adjust(dataframe$p, method = "bonferroni")
	#filter out significant markers for future analysis
	####Should be bonferoni 0.05, testing
	#significant_markers <- filter(dataframe, bon < .05)
	significant_markers <- filter(dataframe, bon < .6)
	print(head(significant_markers))
	#plot out model
	par(cex = 1.5, cex.main = 1.5, pch = 1, lwd=3, mai = c(4,4,3,1))
	manhattan(dataframe, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = graph_name)
	abline(h = 5.6, col = "#659157")
	#give a main dataframe all significant markers
	return(significant_markers)
}

sigMarkers <- data.frame()
for (i in c(2:length(plot_df[1,]))) {
	print(i)
	name <- substring(names(plot_df[i]), 3)
	print(name)
	png(width=2500, height=1500, pointsize = 36, filename = paste('output/plots/', name, '.png', sep=""))
	sig_m <- manhattan_plot(name, plot_df$id, plot_df[,i])
	dev.off()
	sig_m$trait <- name
	sigMarkers <- rbind(sigMarkers, sig_m)
}

###Export all significant markers along with their associated trait
write_csv(sigMarkers, "output/SunRILs_gwas_sig_markers_2022.csv")

#manhattan_plot("Awn GWAS", plot_df$id, plot_df$p_Awns)
#manhattan_plot("Height GWAS", plot_df$id, plot_df$p_Height)
#manhattan_plot("Days to head GWAS", plot_df$id, plot_df$p_days_to_head)

#write_csv(plot_df, "output/SunRILs_gwasOutput_2022.csv")













#------------Hopefully don't need anymore from here down------------#
###Incredibly janky filtering for matrix work, makes model run but results aren't great looking
matrix_df_ind <- data.frame()
#matrix by individual IN family, works better than by family
for (i in unique(t_ph_group$Cross_ID)) {
	a <- filter(t_ph_group, Cross_ID == i)
	b <- sample(unique(a$Entry), length(unique(a$Entry))/2)
	df <- filter(a, (Entry %in% b & Location == "Kinston") | (!(Entry %in% b) & !(Location == "Kinston")))
	matrix_df_ind <- rbind(matrix_df_ind, df)
}
#matrix by family
matrix_df_fam <- data.frame()
rand_cross_id <- sample(unique(t_ph_group$Cross_ID), length(unique(t_ph_group$Cross_ID))/2)
matrix_df_fam <- filter(t_ph_group, (Cross_ID %in% rand_cross_id & Location == "Kinston") | (!(Cross_ID %in% rand_cross_id) & (Location == "Raleigh")))
rownames(matrix_df_fam) <- matrix_df_fam$Entry




#hist(t_ph_group$days_to_head)
###MAKE AND SHOW MODELS
#make covariate model
#familyLocCov <- model.matrix(~ Location + Entry, t_ph_group)
#familyLocCov <- model.matrix(~gsub("-\\d*$", "", t_ph_group$Entry) + t_ph_group$Location)
#testGWAS <- association.test(sunVCF_sync, t_ph_group$days_to_head, X = familyLocCov)
###USING RANDOMLY FILTERED INDIVIDUAL MODEL
familyLocCov <- model.matrix(~ Location + Entry, matrix_df_ind)
familyLocCov <- model.matrix(~gsub("-\\d*$", "", matrix_df_ind$Entry) + matrix_df_ind$Location)
testGWAS_ind <- association.test(sunVCF_sync, matrix_df_ind$days_to_head, X = familyLocCov)
###USING RANDOMLY FILTERED FAMILY MODEL
familyLocCov <- model.matrix(~ Location + Entry, matrix_df_fam)
familyLocCov <- model.matrix(~gsub("-\\d*$", "", matrix_df_fam$Entry) + matrix_df_fam$Location)
testGWAS_fam <- association.test(sunVCF_sync, matrix_df_fam$days_to_head, X = familyLocCov)
###USING ONLY KINSTON
t_ph_group_K <- filter(t_ph_group, Location == "Kinston")
testGWAS_K <- association.test(sunVCF_sync, t_ph_group_K$days_to_head)
###USING ONLY RALEIGH
t_ph_group_R <- filter(t_ph_group, Location == "Raleigh")
testGWAS_R <- association.test(sunVCF_sync, t_ph_group_R$days_to_head)

tests <- list(testGWAS_ind, testGWAS_fam, testGWAS_K, testGWAS_R)

manhattan_plot <- function(graph_name, dataframe) {
	dataframe$LOG <- -log10(dataframe$p)
	dataframe$chr <- str_replace(str_replace(dataframe$id, "^S", ""), "_\\d*$", "")
	dataframe$pos <- as.numeric(str_replace(dataframe$id, "^S\\d[ABD]_", ""))
	dataframe$FDR <- p.adjust(dataframe$p, method = "fdr")
	dataframe$bon <- p.adjust(dataframe$p, method = "bonferroni")
	#plot out model
	manhattan(dataframe, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = graph_name)
	abline(h = 5.6, col = "#659157")
}
#manhatten_plot("SunRILs Heading Date GWAS", tests[[1]])
manhattan_plot("SunRILs Heading Date GWAS random individuals per location", tests[[1]])

#Run statistical tests
testGWAS$LOG <- -log10(testGWAS$p)
testGWAS$chr <- str_replace(str_replace(testGWAS$id, "^S", ""), "_\\d*$", "")
testGWAS$pos <- as.numeric(str_replace(testGWAS$id, "^S\\d[ABD]_", ""))
testGWAS$FDR <- p.adjust(testGWAS$p, method = "fdr")
testGWAS$bon <- p.adjust(testGWAS$p, method = "bonferroni")

#plot out model
manhattan(testGWAS, chrom.col = c("#659157", "#69A2B0", "#FFCAB1"), main = "SunRILs Heading Date GWAS")
abline(h = 5.6, col = "#659157")

#Export significant markers for heading date
