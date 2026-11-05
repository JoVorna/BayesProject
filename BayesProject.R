# Script for the data analysis part // JV

library(data.table)

#data <- fread("clinical_PANCAN_patient_with_followup.tsv")

data <- fread("combined_study_clinical_data.tsv")

grep("disease", colnames(data))

View(data)

# Get wanted columns only
