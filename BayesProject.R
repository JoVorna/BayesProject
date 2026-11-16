# Script for the data analysis part // JV

library(data.table)


# This has all cBioPortal TCGA FireHose legacy cancer studies

# We will only focus on LUAD and LUSC: Both non-small cell lung cancer
data <- fread("combined_study_clinical_data.tsv")

grep("disease", colnames(data))

lung_data <- data[grepl("luad|lusc", data$`Study ID`), ]


names<-names(lung_data)[colMeans(is.na(lung_data)) < 0.2]
lung_data<-subset(lung_data,select=names)




               

      