###Create phenotype files from field data 2022
#Locations: Kinston and Midpines, NC
#Created by Nicolas A. H. Lara
#Created modeled after code written by Noah Dewitt in 2021
#Last edit: 2022-9-14

### LOAD IN PACKAGES
library(tidyverse)

### SET WORKING DIRECTORY ###
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")
#setwd("C:/Users/nalara/Documents/GitHub/GWAS_2022/")


### READ IN PHENOTYPE FILES ###
#growing season traits
K22_pheno <- read.delim("data/phenotype/Kin22-SunRils-T1-T30.csv", sep=",")
R22_pheno <- read.delim("data/phenotype/R22-SunRil-T1-30.csv", sep=",")
#spikelet per spike/infertile spikelet traits
K22_SpS <- read.delim("data/phenotype/SpS_data.xlsx - Kinston.csv", sep=",")
R22_SpS <- read.delim("data/phenotype/SpS_data.xlsx - Raleigh.csv", sep=",")
#seed traits from VIBE
VIBE_file_list <- list.files(path = "data/phenotype/VIBE_reports/", pattern = "*.csv")
VIBE_files <- lapply(paste("data/phenotype/VIBE_reports/", VIBE_file_list, sep=""), read.csv)
VIBE <- bind_rows(VIBE_files)
#outliers to filter
outliers <- read.delim("output/data/outliers.csv", sep=",")

### CLEAN UP FILES, RENAME COLUMNS, ETC. ###
#clean up dataframes, rename columns for field phenotyps
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
K22_pheno <- header.true(K22_pheno)
R22_pheno["Location"] <- "Raleigh"
K22_pheno["Location"] <- "Kinston"
names(K22_pheno)[14] <- "N_rust"
names(K22_pheno)[15] <- "G_rust"
names(K22_pheno)[names(K22_pheno) == "WDR 3/7"] <- "WDR"
names(K22_pheno)[names(K22_pheno) == "Powdery Mildew"] <- "Powdery_mildew"
names(R22_pheno)[names(R22_pheno) == "Powdery.Mildew"] <- "Powdery_mildew"
names(R22_pheno)[names(R22_pheno) == "Plot_Id"] <- "Plot_ID"
names(R22_pheno)[names(R22_pheno) == "Cross_Num"] <- "Cross_ID"
names(R22_pheno)[names(R22_pheno) == "Cross_id"] <- "Cross"
R22_pheno <- R22_pheno[!(R22_pheno$Entry=="Barley"),]
K22_pheno <- K22_pheno[!(K22_pheno$Entry=="Barley"),]
K22_pheno$SNB <- NA
#clean data/rename mislabelled lines
K22_pheno <- mutate(K22_pheno, Entry = if_else(Entry == "UX1994-6.4", "UX1994-64", Entry))
K22_pheno <- mutate(K22_pheno, Entry = if_else(Entry == "UX2029-4.6", "UX2029-46", Entry))
#which(R22_pheno$Entry == "UX2029-4.6")
##Clean up rust scoring from Kinston
K22_pheno$N_rust <- replace(K22_pheno$N_rust, K22_pheno$N_rust=="R",0)
K22_pheno$G_rust <- replace(K22_pheno$G_rust, K22_pheno$G_rust=="R",0)
K22_pheno$G_rust <- replace(K22_pheno$G_rust, K22_pheno$G_rust=="RR",0)
K22_pheno$N_rust<-gsub("R","",as.character(K22_pheno$N_rust))
K22_pheno$N_rust<-gsub("S","",as.character(K22_pheno$N_rust))
K22_pheno$N_rust<-gsub(":","",as.character(K22_pheno$N_rust))
K22_pheno$N_rust<-gsub("M","",as.character(K22_pheno$N_rust))
K22_pheno$N_rust<-gsub("/","",as.character(K22_pheno$N_rust))
K22_pheno$G_rust<-gsub("R/","",as.character(K22_pheno$G_rust))
K22_pheno$N_rust <- as.numeric(K22_pheno$N_rust)
K22_pheno$G_rust <- as.numeric(K22_pheno$G_rust)
K22_pheno$Rust <- rowMeans(K22_pheno[, c("N_rust", "G_rust")], na.rm=TRUE)
R22_pheno$Rust <- NA



### PROCESS SPIKELET PER SPIKE (SpS) DATA FILES ###
K22_SpS$Location <- "Kinston"
R22_SpS$Location <- "Raleigh"
SpS_data <- rbind(K22_SpS, R22_SpS)
#Average awn and awnless lines
SpS_data <- subset(SpS_data, select=-c(Awn_dif, Notes)) %>% group_by(Location, Tray, Cross_ID, Entry) %>% summarise(across(everything(), mean)) %>% ungroup()
SpS_data$spikes_per_plot <- rowSums(!is.na(SpS_data[,grepl("SpS", names(SpS_data))]))
SpS_data$ave_SpS <- apply(SpS_data[,grepl("SpS", names(SpS_data))], 1, mean, na.rm = TRUE)
SpS_data$ave_infert <- apply(SpS_data[,grepl("Inf", names(SpS_data))], 1, mean, na.rm = TRUE)

### PROCESS VIBE DATA FILES###
#subset data to metrics of interest
#VIBE <- VIBE[,c("SampleName","NumberOfParticles", "SampleWeight", "WeightOf1000Particles", "SampleAreaAverage", "WKAreaMedian")]
VIBE <- subset(VIBE, !(SampleName == "test" | SampleName == ""))
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "Kin22-UX2023-2_3-9-4_Ap", "Kin22-UX2023-2_3-9-1_Ap", SampleName))
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "Kin22_UX1993-70_29-17-2_Ap", "Kin22_UX1993-70_26-17-2_Ap", SampleName))
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "_Kin22-UX2023-8_3-20-2_Ap", "Kin22-UX2023-8_3-20-2_Ap", SampleName))
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "R22_UX1992-82_52-17-1_Ap", "R22_UX1992-82_25-17-1_Ap", SampleName))
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "R22_UX2029-204_14-18-47_Am", "R22_UX2029-204_14-18-4_Am", SampleName))
VIBE <- data.frame(lapply(VIBE, function(x) {gsub("22-", "22_", x)}))
#get parents for SampleName splitting
parents <- unique(subset(R22_pheno, Cross_ID == "Parent")$Entry)
VIBE_parents <- VIBE[grep(paste(parents, collapse="|"), VIBE$SampleName),]
VIBE_parents <- VIBE_parents %>% separate(SampleName, c("Location","Entry", "Tray", "Awns"), sep="_")
VIBE_parents$Cross_ID <- "Parent"
#get population entries for SampleName splitting
VIBE_population <- VIBE[!grepl(paste(parents, collapse="|"), VIBE$SampleName),] 
VIBE_population <- VIBE_population %>% separate(SampleName, c("Location", "Entry", "Tray", "Awns"), sep="_")
VIBE_population$Cross_ID <- str_extract(VIBE_population$Entry, "UX[0-9]+")
VIBE <- rbind(VIBE_parents, VIBE_population)
VIBE <- data.frame(lapply(VIBE, function(x) {gsub("Kin22", "Kinston", x)}))
VIBE <- data.frame(lapply(VIBE, function(x) {gsub("R22", "Raleigh", x)}))
VIBE$Tray <- lapply(VIBE$Tray, function(x) {gsub("-", ".", x)})
#Convert columns to numeric and average between awn and awnless in mixed lines
VIBE_num_col <- colnames(subset(VIBE, select=-c(Location, Entry, Tray, Awns, Cross_ID)))
VIBE[, VIBE_num_col] <- apply(VIBE[, VIBE_num_col], 2, as.numeric)
VIBE <- VIBE %>% group_by(Location, Entry, Tray, Cross_ID) %>% summarise(across(everything(), mean)) %>% data.frame %>% ungroup()
#VIBE$Tray <- as.character(VIBE$Tray)

### COMBINE DATAFILES INTO ONE TOTAL PHENOTYPE FILE ###
#Combine phenotype files by location
K22_pheno_subset <- K22_pheno[c("Location", "Tray", "Cross_ID", "Entry", "Awns", "WDR", "flowering", "Powdery_mildew", "SNB", "Height", "Rust")]
R22_pheno_subset <- R22_pheno[c("Location", "Tray", "Cross_ID", "Entry", "Awns", "WDR", "flowering", "Powdery_mildew", "SNB", "Height", "Rust")]
t_ph <- rbind(K22_pheno_subset, R22_pheno_subset)
t_ph <- merge(t_ph, SpS_data[,c("Location","Tray", "spikes_per_plot", "ave_SpS", "ave_infert")], by= c("Location", "Tray"), all.x = TRUE)

#Remove outliers from VIBE file
VIBE <- subset(VIBE, !(Location %in% outliers$Location & Tray %in% outliers$Tray & Cross_ID %in% outliers$Cross_ID & Entry %in% outliers$Entry))

#add VIBE files to phenotype files
VIBE_subset <- VIBE[c("Location", "Tray", "NumberOfParticles", "WeightOf1000Particles","SampleAreaAverage", "WKAreaMedian")]
t_ph <- merge(t_ph, VIBE_subset, by=c("Location", "Tray"), all.x=TRUE)





#add column converting flowering time date to days
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
#change columns to numeric when necessary
num_col <- colnames(subset(t_ph, select=-c(Location, Cross_ID, Entry)))
t_ph[,num_col] <- apply(t_ph[,num_col], 2, as.numeric)
#sapply(t_ph, class)
t_ph$seeds_per_spikelet <- (t_ph$NumberOfParticles/t_ph$spikes_per_plot)/t_ph$ave_SpS

##thin out outliers, as determined by SpS/VIBE analysis of outliers. See data processing file
##This was already done in earlier step
#t_ph <- subset(t_ph, !(Entry %in% outliers$Entry & Location == "Kinston" ))


### CLEAN UP PHENOTYPIC DATA FOR COMBINING WITH GENOTYPE ###
t_ph = subset(t_ph, select=-c(flowering, Tray, spikes_per_plot, NumberOfParticles))

##export raw phenotypic data for viewing
write_csv(t_ph, 'output/data/2022_raw_phenotype.csv')

#check for duplication and lack of replication
genoCounts <- group_by(t_ph, Entry) %>% count() %>% ungroup()
#parLines <- filter(genoCounts, n > 10)$Entry
single_crossLines <- filter(genoCounts, 1 == n)$Entry
#crossLines <- filter(genoCounts, n == 2 )$Entry
#multi_crossLines <- filter(genoCounts, 2 < n & 10 > n )$Entry 

#since we have 23 lines with more than 1 replicate at each location, we average them
#t_ph_group <- t_ph %>% group_by(Location, Entry, Cross_ID) %>% summarize(days_to_head = mean(days_to_head), Height = mean(Height), Awns = mean(Awns), WDR = mean(as.numeric(WDR)), Powdery_mildew = mean(as.numeric(Powdery_mildew)), ave_SpS = mean(ave_SpS), ave_infert = mean(ave_infert)) %>% ungroup()
t_ph_group <- t_ph %>% group_by(Location, Cross_ID, Entry) %>% summarise(across(everything(), mean)) %>% ungroup()

#filter out any unreplicated lines
t_ph_group <- filter(t_ph_group, !(Entry %in% single_crossLines))
t_ph_group <- as.data.frame(t_ph_group)

#Output csv file
write_csv(t_ph_group, "output/data/2022_phenotype.csv")




###CHECK PHENOTYPE DATA TO ENSURE RELATIVE NORMALITY
#Phenotypic data
hist(as.numeric(t_ph$WDR), main="Winter dormancy release 2022")
hist(as.numeric(days), breaks = 40, main = "Flowering date 2022")
hist(as.numeric(t_ph$Powdery_mildew), main="Powdery mildew 2022")
hist(as.numeric(t_ph$SNB), main="Septoria nodorum blotch 2022")
hist(as.numeric(t_ph$ave_SpS), main="Spikelet per Spike")
hist(as.numeric(t_ph$ave_infert), main="Infertile Spikelet per Spike")
hist(t_ph$Height, main="Height")
#VIBE data
hist(as.numeric(VIBE$SampleWeight))
hist(as.numeric(VIBE$NumberOfParticles))


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


#looking for gremlins: helpful functions
#subset(t_ph ,duplicated(t_ph))
#TrayCounts <- group_by(VIBE, Tray) %>% count() %>% ungroup()
#filter(TrayCounts, n > 2)
#subset(t_ph, Tray == "29.17.2")
#VIBE_tray <- paste(VIBE$Location, substr(VIBE$Tray, 1,2), sep = "_")
#SpS_tray <- paste(SpS_data$Location, substr(SpS_data$Tray, 1,2), sep = "_")
# VIBE_tray <- paste(VIBE$Location, VIBE$Tray, sep = "_")
# SpS_tray <- paste(SpS_data$Location, SpS_data$Tray, sep = "_")

# VIBE_tray <- unique(VIBE_tray)
# SpS_tray <- unique(SpS_tray)
# for (i in SpS_tray) {
	# if (i %in% VIBE_tray) {
		# #print(paste(i, "VIBE'd", sep=" "))
	# } else {
		# print(paste(i, "missing", sep=" "))
	# }	
# }
# for (i in VIBE_tray) {
	# if (i %in% SpS_tray) {
		# #print(paste(i, "VIBE'd", sep=" "))
	# } else {
		# print(paste(i, "mislabeled", sep=" "))
	# }	
# }
# orderVIBE <- VIBE[order(unlist(VIBE$Tray)),]
#PopCounts <- group_by(t_ph, Location, Cross_ID) %>% count() %>% ungroup()
