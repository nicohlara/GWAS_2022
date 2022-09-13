###Evaluate rust scores from 2022 data
#Locations: Kinston 
#Created by Nicolas A. H. Lara
#Last edit: 2022-8-18

### LOAD IN PACKAGES
library(tidyverse)

K22_pheno <- read.delim("data/Kin22-SunRils-T1-T30.csv", sep=",", na.strings="")

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
K22_pheno <- header.true(K22_pheno)
K22_pheno <- K22_pheno[!(K22_pheno$Entry=="Barley"),]
names(K22_pheno)[14] <- "N_rust"
names(K22_pheno)[15] <- "G_rust"
rust_df <- subset(K22_pheno, select=c(Cross_ID, Tray, Entry, N_rust, G_rust))

###COMPARING COMPLETE CASES OF DATA BETWEEN RATING SCHEMES
comparison_rust <- rust_df[complete.cases(rust_df),]
###at this stage, we have a dataframe with only lines that both G and N evaluated. 
hist(as.numeric(subset(comparison_rust, N_rust=="S")$G_rust))
###N's 'S' rating doesn't correspond well to anything in G's scores
hist(as.numeric(subset(comparison_rust, N_rust=="R")$G_rust))
###N's 'R' score overwhelmingly corresponds to G's '0' score (or 'R'), so we will replace 'R' with 0
comparison_rust$N_rust <- replace(comparison_rust$N_rust, comparison_rust$N_rust=="R",0)
subset(comparison_rust, G_rust=="R")
###Likewise, G's 'R' score corresponds overwhelmingly to N's 0 or R scores, so replace again:
comparison_rust$G_rust <- replace(comparison_rust$G_rust, comparison_rust$G_rust=="R",0)
comparison_rust$G_rust <- replace(comparison_rust$G_rust, comparison_rust$G_rust=="RR",0)
###Some scores are preceded by 'R' or 'S', R < 30 and S > 30. Remove both:
comparison_rust$N_rust<-gsub("R","",as.character(comparison_rust$N_rust))
comparison_rust$N_rust<-gsub("S","",as.character(comparison_rust$N_rust))
comparison_rust$N_rust<-gsub(":","",as.character(comparison_rust$N_rust))
comparison_rust$G_rust<-gsub("R/","",as.character(comparison_rust$G_rust))
plot(comparison_rust$N_rust, comparison_rust$G_rust)
###The plotted out cleaned data doesn't have much correlation, which isn't great. We'll average what we can in the final run

rust_df$N_rust <- replace(rust_df$N_rust, rust_df$N_rust=="R",0)
rust_df$G_rust <- replace(rust_df$G_rust, rust_df$G_rust=="R",0)
rust_df$G_rust <- replace(rust_df$G_rust, rust_df$G_rust=="RR",0)
rust_df$N_rust<-gsub("R","",as.character(rust_df$N_rust))
rust_df$N_rust<-gsub("S","",as.character(rust_df$N_rust))
rust_df$N_rust<-gsub(":","",as.character(rust_df$N_rust))
rust_df$N_rust<-gsub("M","",as.character(rust_df$N_rust))
rust_df$N_rust<-gsub("/","",as.character(rust_df$N_rust))
#rust_df$N_rust<-gsub("",NA,as.character(rust_df$N_rust))
rust_df$G_rust<-gsub("R/","",as.character(rust_df$G_rust))

rust_df$N_rust <- as.numeric(rust_df$N_rust)
rust_df$G_rust <- as.numeric(rust_df$G_rust)

rust_df$rust_score <- rowMeans(rust_df[, c(4:5)], na.rm=TRUE)
