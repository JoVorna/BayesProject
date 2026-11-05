# Script for the data analysis part // JV

library(data.table)


# This has all cBioPortal TCGA FireHose legacy cancer studies

# We will only focus on LUAD and LUSC: Both non-small cell lung cancer
data <- fread("combined_study_clinical_data.tsv")

grep("disease", colnames(data))

View(data)

# colnames(data)

# Get wanted and interesting colnames

data <- subset(data, select = c(
  "Study ID",
  "Patient ID",
  "Sample ID",
  "Cancer Type",
  "Oncotree Code",
  "Diagnosis Age",
  "Sex",
  "Alcohol Exposure Intensity",
  "Alcohol History Documented",
  "Neoplasm Disease Stage American Joint Committee on Cancer Code",
  "Neoplasm American Joint Committee on Cancer Clinical Group Stage",
  "Neoplasm American Joint Committee on Cancer Clinical Distant Metastasis M Stage",
  "Neoplasm American Joint Committee on Cancer Clinical Regional Lymph Node N Stage",
  "Neoplasm Histologic Grade",
  "Neoplasm Histologic Type Name",
  "Performance Status",
  "Karnofsky Performance Score",
  "Tumor burden pathologic",
  "Metastatic disease confirmed",
  "Lymphovascular invasion indicator",
  "Mutation Count",
  "TMB (nonsynonymous)",
  "Fraction Genome Altered",
  "KRAS Mutation",
  "BRAF Gene Analysis Result",
  "ALK Translocation Status",
  "Tumor Infiltrating Lymphocytes",
  "Overall Survival (Months)",
  "Overall Survival Status"
))



grep("Months", colnames(data))

# Get wanted columns only

luad_cases <- which(grepl("luad", c(data$`Study ID`)))
lusc_cases <- which(grepl("lusc", c(data$`Study ID`)))


luad <- data[luad_cases, ]
lusc <- data[lusc_cases, ]

               

      