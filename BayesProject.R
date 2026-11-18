# Script for the data analysis part // JV

library(data.table)
library(tidyr)

# This has all cBioPortal TCGA FireHose legacy cancer studies

# We will only focus on LUAD and LUSC: Both non-small cell lung cancer
data <- fread("combined_study_clinical_data.tsv")

grep("disease", colnames(data))

lung_data <- data[grepl("luad|lusc", data$`Study ID`), ]
setnames(lung_data, make.unique(names(lung_data)))  #fix duplicate columns
#only lung data:

cols_to_keep<-names(lung_data)[colMeans(is.na(lung_data)) < 0.2]
lung_data <- subset(lung_data,select=cols_to_keep)

# " The censoring variable 
#should contain the values 'left', 'none', and 'right' (or equivalenty -1, 0, and 1)"
lung_data$cens<-ifelse(lung_data$`Overall Survival Status`=="1:DECEASED",0,1)

#data cleaning
model_cols <- c(
  "Patient ID",
  "Study ID",
  "Overall Survival (Months)",
  "Overall Survival Status",
  "Diagnosis Age",
  "Neoplasm Disease Stage American Joint Committee on Cancer Code.1",
  "Sex",
  "Patient Smoking History Category",
  "Fraction Genome Altered",
  "Prior Cancer Diagnosis Occurence",
  "Neoadjuvant Therapy Type Administered Prior To Resection Text",
  "Year Cancer Initial Diagnosis",
  "cens" #cencoring value
)
lung_model<-subset(lung_data,select = model_cols)
lung_model<- lung_model %>% drop_na()   #can figure out imputations but this still leaves 790 samples

#remove non sensible values
lung_model<-subset(lung_model, lung_model$`Overall Survival (Months)`!=0)
         