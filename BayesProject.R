# Script for the data analysis part // JV

library(data.table)
library(tidyr)
library(brms)
library(posterior)
library(bayesplot)
library(ggplot2)

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


# About the clinical file:

# https://groups.google.com/g/cbioportal/c/irEXZRj9Who
# Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
# Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
# Current reformed smoker for > 15 years (greater than 15 years) = 3
# Current reformed smoker for ≤15 years (less than or equal to 15 years) = 4
# Current reformed smoker, duration not specified = 5


# Add a column called "Smoking_history" which has 3 categories only
# non-smoker, current-smoker and recurrent-smoker based on the Info above


lung_model$Smoking_history <- c("")

for(i in 1:nrow(lung_model)) {
  if (lung_model$`Patient Smoking History Category`[i] == 1) {
    lung_model$Smoking_history[i] <- "non-smoker"
  } else if (lung_model$`Patient Smoking History Category`[i] == 2) {
    lung_model$Smoking_history[i] <- "current-smoker"
  } else {
    lung_model$Smoking_history[i] <- "recurrent-smoker"
  }
}

# Group these together so that letters are removed
unique(lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`)

# regex // Claude code, works now

lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1` <-
  sub("(Stage [IV]+)[AB]$", "\\1",
      lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`)

unique(lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`)

# Rename awkward named cols

lung_model$os_months <- lung_model$`Overall Survival (Months)`
lung_model$Age <- lung_model$`Diagnosis Age`
lung_model$Stage <- lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`


# 2 different lung cancer types: Lung adenocarcinoma (LUAD) and Lung squamous cell carcinoma (LUSC)

luad_model <- subset(lung_model, lung_model$`Study ID` == "luad_tcga")
lusc_model <- subset(lung_model, lung_model$`Study ID` == "lusc_tcga")


# Now fitting


#?brm
colnames(lung_model)


# Use os_months and cens

# Coefficients are then the variables (cols) such as Age, Sex, Stage, anything we want



lognormal <- brm(
  os_months | cens(cens) ~ 1 + Smoking_history + Age + Sex + Stage,
  data   = lung_model,
  family = lognormal(),
  prior  = prior(normal(3, 2), class = b),
  chains = 4, cores = 4,
  backend = "cmdstanr"
)

summary(lognormal)
mcmc_trace(lognormal)


weibull <- brm(
  os_months | cens(cens) ~ 1 + Smoking_history,
  data   = lung_model,
  family = "weibull",
  prior  = prior(normal(3, 2), class = b),
  chains = 4, cores = 4,
  backend = "cmdstanr"
)

mcmc_trace(weibull)


# Priorit betoille (coefficients) eli esimerkiksi tässä SMoking_historyssa, kuinka paljon luulemme että tupakointi status 0-2 vaikuttaa selviytytyvyyteen
# Siitä sitten joku jakauma





# Assignment 8 code:




lung_model$Age_c <- (lung_model$Age - 60) / 10


# Relevel smoking history so that the "reference" is non-smoker

lung_model$Smoking_history <- factor(
  lung_model$Smoking_history,
  levels = c("non-smoker", "current-smoker", "recurrent-smoker")
)

levels(lung_model$Smoking_history)   # first one = reference

colnames(lung_model)

lung_model$Genome_altered <- lung_model$`Fraction Genome Altered`
lung_model$prior_diagnosis <- lung_model$`Prior Cancer Diagnosis Occurence`



# center = F / T
formula_lognormal <- bf(
  os_months | cens(cens) ~ 1 + Smoking_history + Age_c + Stage + Sex,
  family = "lognormal",
  center = TRUE
)




get_prior(formula_lognormal, data = lung_model)



lognormal_fit_priors <- c(
  prior(
    normal(4, 1),
    class = "Intercept"#,
    #coef = "Intercept"
  ),
  prior(
    normal(-0.7, 0.5),
    class = "b",
    coef = "Smoking_historycurrentMsmoker"
  ),
  prior(
    normal(-0.3, 0.5),
    class = "b",
    coef = "Smoking_historyrecurrentMsmoker"
  ),
  prior(
    normal(0, 1),
    class = "sigma"
  ),
  prior(
    normal(-0.2, 0.5),
    class = "b",
    coef = "StageStageII"
  ),
  prior(
    normal(-0.3, 0.5),
    class = "b",
    coef = "StageStageIII"
  ),
  prior(
    normal(-0.4, 0.5),
    class = "b",
    coef = "StageStageIV"
  ),
  prior(
    normal(-0.1, 0.3),
    class = "b",
    coef = "Age_c"
    ),
  prior(
    normal(0, 0.5),
    class = "b",
    coef = "SexMale")
  )




# 

# lognormal_fit_priors <- c(
#   prior(
#     normal(4, 1),
#     class = "b",
#     coef = "Intercept"
#   ),
#   prior <- c(
#     set_prior(
#       "normal(-0.3, 1)",
#       class = "b",
#       coef = "SexMale"
#     )
#   )
# )

lognormal_fit <- brm(
  formula = formula_lognormal,
  prior = lognormal_fit_priors,
  data = lung_model,
  family = lognormal(),
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  cores = 8
)


# Here you can try different data (like luad or lusc) instead of whole lung data
# lognormal_fit <- brm(
#   formula = formula_lognormal,
#   prior = lognormal_fit_priors,
#   data = luad_model,
#   family = lognormal(),
#   backend = "cmdstanr",
#   chains = 4,
#   cores = 6
# )

#library(bayesplot); library(posterior); library(brms)


# Claude code 100% from now on
# Similar visualization funcs available in assignments (USE THOSE LATER ON!)

mcmc_trace(lognormal_fit, regex_pars = "^b_|sigma")
mcmc_rank_overlay(lognormal_fit, regex_pars = "^b_|sigma")
summarise_draws(as_draws_df(lognormal_fit)) |> dplyr::select(variable, rhat, ess_bulk)





library("ggplot2")
draws <- as_draws_df(lognormal_fit)
b_cols <- grep("^b_", names(draws), value = TRUE)
b_cols <- setdiff(b_cols, "b_Intercept")

tr <- exp(as.matrix(draws[, b_cols]))
colnames(tr) <- sub("^b_", "", colnames(tr))


mcmc_intervals(tr, prob = 0.5, prob_outer = 0.95) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Time ratio (>1 = longer survival)")

mcmc_areas(tr, prob = 0.95) + geom_vline(xintercept = 1, linetype = "dashed")

