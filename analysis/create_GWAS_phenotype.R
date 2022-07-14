###Create phenotype files from field data 2022
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modelled after code written by Noah Dewitt in 2021
#Last edit: 2022-7-11

library(tidyverse)

#set working directory
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")

###READ IN PHENOTYPE FILES, CLEAN COLUMN NAMES, COMBINE LOCATIONS
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
#Change awn from + to 1,0 for coding purposes
awn_vector <- ifelse(t_ph$Awns=="+", 1,0)
t_ph$Awns <- awn_vector

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

write_csv(t_ph_group, "output/data/2022_phenotype.csv")







###CHECK PHENOTYPE DATA TO ENSURE RELATIVE NORMALITY
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

