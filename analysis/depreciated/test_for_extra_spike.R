### Test uncertain extra spike data to normalize for analysis ###
#Created by Nicolas A. H. Lara
#Last edit: 2022-7-20

### LOAD IN PACKAGES
library(tidyverse)

### SET WORKING DIRECTORY ###
setwd("/Users/nico/Documents/GitHub/GWAS_2022/")

### READ IN RAW DATA FILES ###
#spikelet per spike/infertile spikelet traits
K22_SpS <- read.delim("data/SpS_data.xlsx - Kinston.csv", sep=",")
R22_SpS <- read.delim("data/SpS_data.xlsx - Raleigh.csv", sep=",")
#seed traits from VIBE
VIBE_file_list <- list.files(path = "data/VIBE_reports/", pattern = "*.csv")
VIBE_files <- lapply(paste("data/VIBE_reports/", VIBE_file_list, sep=""), read.csv)
VIBE <- bind_rows(VIBE_files)

### PROCESS SPIKELET PER SPIKE (SpS) DATA FILES ###
#Average awn and awnless lines
K22_SpS$Location <- "Kinston"
R22_SpS$Location <- "Raleigh"
SpS_data <- rbind(K22_SpS, R22_SpS)
SpS_data <- subset(SpS_data, select=-c(Awn_dif, Notes)) %>% group_by(Location, Tray, Cross_ID, Entry) %>% summarise(across(everything(), mean)) %>% ungroup()
SpS_data$spikes_per_plot <- rowSums(!is.na(SpS_data[,grepl("SpS", names(SpS_data))]))
SpS_data$ave_SpS <- apply(SpS_data[,grepl("SpS", names(SpS_data))], 1, mean, na.rm = TRUE)
SpS_data$ave_infert <- apply(SpS_data[,grepl("Inf", names(SpS_data))], 1, mean, na.rm = TRUE)

### PROCESS VIBE DATA FILES###
#subset data to metrics of interest
#VIBE <- VIBE[,c("SampleName","NumberOfParticles", "SampleWeight", "WeightOf1000Particles", "SampleAreaAverage", "WKAreaMedian")]
VIBE <- mutate(VIBE, SampleName = if_else(SampleName == "_Kin22-UX2023-8_3-20-2_Ap", "Kin22-UX2023-8_3-20-2_Ap", SampleName))
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


#add VIBE file to SpS file
VIBE_subset <- VIBE[c("Location", "Tray", "SampleWeight","NumberOfParticles", "WeightOf1000Particles","SampleAreaAverage", "WKAreaMedian")]
seed_data <- merge(SpS_data, VIBE_subset, by=c("Location", "Tray"))


#A measurement error occurred during early sampling
#This applies only to bags that contained 6 heads, which was a very small subset of bags in affected trays
#This error affects about 6 trays (Kinston 2,3,8,9,14, and half of 26 and 29) covering parts (but not all of) of UX1991, UX1992, UX1993, and UX2023 populations from Kinston
#This error does not affect Raleigh data
#These were threshed without distinguishing between bags with 5 or 6 spikes
#Of bags with 6 spikes, only 5 were counted
#In later trays, either all 6 spikes were counted or the 6th spike was discarded if it was broken or atypical
#affected_pop <- subset(SpS_data, Thresh_dif_chance == 1)
#affected_pop <- affected_pop %>% separate(Tray, c("Tray","Range", "Row"))
#fam_counts <- group_by(affected_pop, Tray) %>% count() %>% ungroup()

error_find_df <- subset(seed_data, Thresh_dif_chance == 1)
num_col <- colnames(subset(error_find_df, select=-c(Location, Tray, Cross_ID, Entry)))
error_find_df[,num_col] <- apply(error_find_df[,num_col], 2, as.numeric)

outliers <- data.frame()
### RUNNING HISTOGRAMS TO LOOK FOR OBVIOUS DEVIATIONS ###
hist(error_find_df$spikes_per_plot, main="R calculated spikes per plot")
#no big outliers here, since all bags with 6 SpS datapoints have been accounted for already
hist(error_find_df$SampleWeight, main="Sample Weight")
#This looks pretty normal, there's not an obvious bump
hist(error_find_df$NumberOfParticles, main="Number of seeds")
#here we do see a small bump outside of normal, these are the likely candidates for data smoothing. Let's store those
outlier_df <- subset(error_find_df, NumberOfParticles >500)
outlier_df$Trait <- "Seed number"
outliers <- rbind(outliers, outlier_df[,c("Location", "Tray", "Cross_ID", "Entry", "Trait")])

error_find_df$seeds_per_spikelet <- error_find_df$NumberOfParticles/(error_find_df$ave_SpS * error_find_df$spikes_per_plot)
hist(error_find_df$seeds_per_spikelet, main="Seeds per spikelet")
#another couple of bumps outside of normal, one at 6 and one near 8
outlier_df <- subset(error_find_df, seeds_per_spikelet >5.2)
outlier_df$Trait <- "Seeds per spikelet"
outliers <- rbind(outliers, outlier_df[,c("Location", "Tray", "Cross_ID", "Entry", "Trait")])

#Let's check some weight metrics. First, weight per spikelet
error_find_df$weight_per_spikelet <- error_find_df$SampleWeight/(error_find_df$ave_SpS * error_find_df$spikes_per_plot)
hist(error_find_df$weight_per_spikelet, main="Weight per spikelet")
#No obvious outliers there

print(sort(outliers$Tray))
write_csv(unique(outliers[,c("Location", "Tray", "Cross_ID", "Entry")]), "output/data/outliers.csv")