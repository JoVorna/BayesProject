######################################################################
#####################  Bayesian Data Analysis Project ################
######################################################################


# Set your own path for the data.
setwd("/Users/joonasvornanen/MRD2026/MRD2026/BayesProject")


# Load needed libraries here
library(data.table)
library(tidyr)
library(brms)
library(posterior)
library(bayesplot)
library(ggplot2)
library(cmdstanr)
library("priorsense")

# This has all cBioPortal TCGA FireHose legacy cancer studies
# We will only focus on Lung cancer



######################################################################
#####################  Data loading and cleansing ####################
######################################################################



data <- fread("combined_study_clinical_data.tsv")

grep("disease", colnames(data))

# Luad and lusc are both lung cancer

lung_data <- data[grepl("luad|lusc", data$`Study ID`), ]
setnames(lung_data, make.unique(names(lung_data)))

#only lung data:
cols_to_keep<-names(lung_data)[colMeans(is.na(lung_data)) < 0.2]
lung_data <- subset(lung_data,select=cols_to_keep)

# " The censoring variable 
#should contain the values 'left', 'none', and 'right' (or equivalenty -1, 0, and 1)"
lung_data$cens<-ifelse(lung_data$`Overall Survival Status`=="1:DECEASED",0,1)

#data cleaning / take only the columns we consider using later on (Not necessarily everything is used as coefficients though)
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

# regex // Claude code for cleaning regex

lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1` <-
  sub("(Stage [IV]+)[AB]$", "\\1",
      lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`)

unique(lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`)

# Rename awkward named cols

lung_model$os_months <- lung_model$`Overall Survival (Months)`
lung_model$Age <- lung_model$`Diagnosis Age`
lung_model$Stage <- lung_model$`Neoplasm Disease Stage American Joint Committee on Cancer Code.1`


# Relevel smoking history so that the "reference" is non-smoker for 

lung_model$Smoking_history <- factor(
  lung_model$Smoking_history,
  levels = c("non-smoker", "current-smoker", "recurrent-smoker")
)

levels(lung_model$Smoking_history)   # first one = reference

colnames(lung_model)

lung_model$Genome_altered <- lung_model$`Fraction Genome Altered`

unique(lung_model$`Prior Cancer Diagnosis Occurence`)
lung_model$prior_diagnosis <- factor(
  ifelse(startsWith(lung_model$`Prior Cancer Diagnosis Occurence`, "Yes"), "Yes", "No"),
  levels = c("No", "Yes")
)


# 2 different lung cancer types: Lung adenocarcinoma (LUAD) and Lung squamous cell carcinoma (LUSC)
# Could be used in the future, but for this purpose we only stick with general lung_model


###luad_model <- subset(lung_model, lung_model$`Study ID` == "luad_tcga")
###lusc_model <- subset(lung_model, lung_model$`Study ID` == "lusc_tcga")


# Now model fitting:


# Use os_months and cens

# Coefficients are then the variables (cols) such as Age, Sex, Stage, anything we want



######################################################################
#####################  Lognormal MODEL ###############################
######################################################################


formula_lognormal <- bf(
  os_months | cens(cens) ~ 1 + Smoking_history + Age + Stage + Sex + prior_diagnosis + Genome_altered,
  family = "lognormal",
  center = FALSE
)

# Look for prior names (these are weirdly named)
get_prior(formula_lognormal, data = lung_model)


# Set priors for coefficients

# These are based on our own thinking and some prior knowledge about lung cancer

lognormal_fit_priors <- c(
  prior(
    normal(4, 1),
    class = "b",
    coef = "Intercept"
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
    student_t(3, 1.5, 1),
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
    coef = "Age"
    ),
  prior(
    normal(-0.2, 0.5),
    class = "b",
    coef = "SexMale")
  ,
prior(
  normal(-0.5, 0.5),
  class = "b",
  coef = "prior_diagnosisYes")
,
prior(
  normal(-0.3, 0.5),
  class = "b",
  coef = "Genome_altered"
  )
)


# Fit lognormal model with brm()
lognormal_fit <- brm(
  formula = formula_lognormal,
  prior = lognormal_fit_priors,
  data = lung_model,
  family = lognormal(),
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  cores = 4
)
# Summary of the fit model (Rhat, ESS etc.)

summary(lognormal_fit)

# Trace plot to see converge

mcmc_trace(lognormal_fit)

# For a few parameters only

# Do the posterior predictive check, KDE and Kaplan-Meyer overlay (This KM-plot is widely used in survival analysis, the one this whole data is based on)

pp_check(lognormal_fit, ndraws = 100) + scale_x_log10()
pp_check(lognormal_fit, type = "km_overlay", ndraws = 100)





######################################################################
#####################  Weibull MODEL #################################
######################################################################

# Same things, model formula, priors, fit:


formula_weibull <- bf(
  os_months | cens(cens) ~ 1 + Smoking_history + Age + Stage + Sex + prior_diagnosis + Genome_altered,
  family = "weibull",
  center = FALSE
)

#get_prior(formula_weibull, data = lung_model)


# Priors for Weibull

weibull_fit_priors <- c(
  prior(
    normal(4, 1),
    class = "b",
    coef = "Intercept"
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
    coef = "Age"
  ),
  prior(
    normal(-0.2, 0.5),
    class = "b",
    coef = "SexMale")
  ,
  prior(
    normal(-0.5, 0.5),
    class = "b",
    coef = "prior_diagnosisYes")
  ,
  prior(
    normal(-0.3, 0.5),
    class = "b",
    coef = "Genome_altered"
  ),
  prior(
    student_t(3, 1.5, 1),
    class = "shape"
  )
)

#?student_t

weibull_fit <- brm(
  formula = formula_weibull,
  prior = weibull_fit_priors, # For the test purpose use the same priors
  data = lung_model,
  family = lognormal(),
  backend = "cmdstanr",
  chains = 4,
  iter = 4000,
  cores = 8
)

# Summary weibull

summary(weibull_fit)


# Trace for convergence checking

mcmc_trace(weibull_fit)


# Posterior predictive check for Weibull + KM overlay, as for lognormal

pp_check(weibull_fit, ndraws = 100) + scale_x_log10()
pp_check(weibull_fit, type = "km_overlay", ndraws = 100)

#?ppc_km_overlay()


######################################################################
#####################  Prior sensitivity analysis ####################
######################################################################


# Lognormal
powerscale_plot_dens(lognormal_fit)


# Weibull
powerscale_plot_dens(weibull_fit)



######################################################################
#####################  Loo comparison #################################
######################################################################

# Fit loo for lognormal model

loo_fit_lognormal <- add_criterion(
  lognormal_fit,
  criterion = "loo"
)

# Fit loo for weibull

loo_fit_weibull<- add_criterion(weibull_fit,criterion = "loo")


# Compare and print out the results
loo_compare(
  loo_fit_lognormal,
  loo_fit_weibull
)


#sessionInfo()
# End of script.
