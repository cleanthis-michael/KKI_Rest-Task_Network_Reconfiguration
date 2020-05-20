#### KKI GNG Task Analyses ####
## May 2020, CM ##
## This script compares behavioral metrics of performance between ADHD and TD children on GNGs and GNGr, including meanRT, sdRT, cvRT, mu, sigma, tau, commission rates, and omission rates ##


rm(list=ls())

# Set current directory as the working directory
library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
print( getwd() )


#### A. Load in the csv file with the behavioral data, merge it with the csv file containing demographic information (age, gender, diagnosis), and select only behavioral metrics combining across runs (named 'overall' in the file containing the behavioral data) ####
require(readxl)
require(dplyr)

demog_data <- read.csv('subject_demographics.csv', sep = ',', stringsAsFactors = FALSE)

all_behav_data <- read_excel('./behavioral_data/GNG_BehavSummary.xlsx')
all_behav_data$sub[all_behav_data$sub < 1000] <- paste('0', all_behav_data$sub, sep = '')
all_behav_data$sub <- paste('sub-', all_behav_data$sub, sep = '')

behav_data <- merge(all_behav_data, demog_data, by = ('sub'))
behav_data <- behav_data %>% dplyr::select(sub, diag:gender, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate')
behav_data$`Overall_MeanRT_diff` <- behav_data$'Overall_Repeat_Mean RT' - behav_data$'Overall_Simple_Mean RT'
behav_data$`Overall_SDRT_diff` <- behav_data$'Overall_Repeat_Std RT' - behav_data$'Overall_Simple_Std RT'
behav_data$`Overall_CVRT_diff` <- behav_data$'Overall_Repeat_CV RT' - behav_data$'Overall_Simple_CV RT'
behav_data$`Overall_SDRT_diff` <- behav_data$'Overall_Repeat_Std RT' - behav_data$'Overall_Simple_Std RT'
behav_data$`Overall_Mu_diff` <- behav_data$'Overall_Repeat_Mu' - behav_data$'Overall_Simple_Mu'
behav_data$`Overall_Sigma_diff` <- behav_data$'Overall_Repeat_Sigma' - behav_data$'Overall_Simple_Sigma'
behav_data$`Overall_Tau_diff` <- behav_data$'Overall_Repeat_Tau' - behav_data$'Overall_Simple_Tau'
behav_data$`Overall_Commission_diff` <- behav_data$'Overall_Repeat_Commision Rate' - behav_data$'Overall_Simple_Commision Rate'
behav_data$`Overall_Omission_diff` <- behav_data$'Overall_Repeat_Ommision Rate' - behav_data$'Overall_Simple_Ommision Rate'


# Create an RData file with all the relevant data for only the relevant subjects #
save(behav_data, file = "behav_data.RData") 


#### B. Running t-tests comparing behavioral performance between ADHD versus TD children without controlling for covariates (e.g., age) ####

#### B1. Loading the RData file containing the behavioral data ####
behav_data <- get(load('behav_data.RData'))


#### B2. Creating an empty data frame to save the results ####
t_results_df <- data.frame(comparison = NA,
                         mean_x = NA,
                         sd_x = NA,
                         mean_y = NA,
                         sd_y = NA,
                         t.value = NA,
                         df = NA,
                         p.value = NA)

# Create empty lists that will be updated after each test and, at the end, will be used to update the t_results_df before saving it as a csv file #
comparison_list = c()
mean_x_list = c()
sd_x_list = c()
mean_y_list = c()
sd_y_list = c()
t.value_list = c()
df_list = c()
p.value_list = c()

options(scipen = 999)  # set scientific notation to FALSE so that all decimal points are depicted (rather than presenting them in terms of e.g. e-05)


### B3. Running the t-tests ###
## GNGs - ADHD vs. TD ##
GNGs_meanRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_meanRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_meanRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_meanRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_meanRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_meanRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_meanRT_ADHD_vs_TD$p.value)

GNGs_sdRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_sdRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_sdRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_sdRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_sdRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_sdRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_sdRT_ADHD_vs_TD$p.value)

GNGs_cvRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_cvRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_cvRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_cvRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_cvRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_cvRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_cvRT_ADHD_vs_TD$p.value)

GNGs_mu_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Mu`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_mu_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_mu_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_mu_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_mu_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_mu_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_mu_ADHD_vs_TD$p.value)

GNGs_sigma_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_sigma_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_sigma_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_sigma_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_sigma_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_sigma_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_sigma_ADHD_vs_TD$p.value)

GNGs_tau_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Tau`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_tau_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_tau_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_tau_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_tau_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_tau_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_tau_ADHD_vs_TD$p.value)

GNGs_commission_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Commision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Commision`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_commission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_commission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_commission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_commission_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_commission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_commission_ADHD_vs_TD$p.value)

GNGs_omission_ADHD_vs_TD <- t.test(behav_data$`Overall_Simple_Ommision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Simple_Ommision`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_omission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGs_omission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGs_omission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_omission_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGs_omission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGs_omission_ADHD_vs_TD$p.value)


## GNGr - ADHD vs. TD ##
GNGr_meanRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_meanRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_meanRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_meanRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_meanRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_meanRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_meanRT_ADHD_vs_TD$p.value)

GNGr_sdRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_sdRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_sdRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_sdRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_sdRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_sdRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_sdRT_ADHD_vs_TD$p.value)

GNGr_cvRT_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_cvRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_cvRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_cvRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_cvRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_cvRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_cvRT_ADHD_vs_TD$p.value)

GNGr_mu_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_mu_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_mu_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_mu_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_mu_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_mu_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_mu_ADHD_vs_TD$p.value)

GNGr_sigma_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_sigma_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_sigma_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_sigma_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_sigma_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_sigma_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_sigma_ADHD_vs_TD$p.value)

GNGr_tau_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_tau_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_tau_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_tau_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_tau_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_tau_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_tau_ADHD_vs_TD$p.value)

GNGr_commission_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Commision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Commision`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_commission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_commission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_commission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_commission_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_commission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_commission_ADHD_vs_TD$p.value)

GNGr_omission_ADHD_vs_TD <- t.test(behav_data$`Overall_Repeat_Ommision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Ommision`[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGr_omission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, GNGr_omission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, GNGr_omission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGr_omission_ADHD_vs_TD$statistic)
df_list <- c(df_list, GNGr_omission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, GNGr_omission_ADHD_vs_TD$p.value)


## GNGs vs. GNGr - TD ##
GNGs_vs_GNGr_meanRT_TD <- t.test(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_meanRT_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_meanRT_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_meanRT_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_meanRT_TD$p.value)

GNGs_vs_GNGr_sdRT_TD <- t.test(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_sdRT_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_sdRT_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_sdRT_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_sdRT_TD$p.value)

GNGs_vs_GNGr_cvRT_TD <- t.test(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_cvRT_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_cvRT_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_cvRT_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_cvRT_TD$p.value)

GNGs_vs_GNGr_mu_TD <- t.test(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_mu_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_mu_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_mu_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_mu_TD$p.value)

GNGs_vs_GNGr_sigma_TD <- t.test(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_sigma_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_sigma_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_sigma_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_sigma_TD$p.value)

GNGs_vs_GNGr_tau_TD <- t.test(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_tau_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_tau_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_tau_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_tau_TD$p.value)

GNGs_vs_GNGr_commission_TD <- t.test(behav_data$`Overall_Simple_Commision`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Commision`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_commission_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_commission_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_commission_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_commission_TD$p.value)

GNGs_vs_GNGr_omission_TD <- t.test(behav_data$`Overall_Simple_Ommision`[behav_data$diag == 'TD'], behav_data$`Overall_Repeat_Ommision`[behav_data$diag == 'TD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_omission_TD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'TD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'TD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'TD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_omission_TD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_omission_TD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_omission_TD$p.value)


## GNGs vs. GNGr - ADHD ##
GNGs_vs_GNGr_meanRT_ADHD <- t.test(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_meanRT_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mean RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mean RT`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_meanRT_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_meanRT_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_meanRT_ADHD$p.value)

GNGs_vs_GNGr_sdRT_ADHD <- t.test(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_sdRT_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Std RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Std RT`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_sdRT_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_sdRT_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_sdRT_ADHD$p.value)

GNGs_vs_GNGr_cvRT_ADHD <- t.test(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_cvRT_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_CV RT`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_CV RT`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_cvRT_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_cvRT_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_cvRT_ADHD$p.value)

GNGs_vs_GNGr_mu_ADHD <- t.test(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_mu_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Mu`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Mu`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_mu_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_mu_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_mu_ADHD$p.value)

GNGs_vs_GNGr_sigma_ADHD <- t.test(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_sigma_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Sigma`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Sigma`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_sigma_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_sigma_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_sigma_ADHD$p.value)

GNGs_vs_GNGr_tau_ADHD <- t.test(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_tau_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Tau`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Tau`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_tau_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_tau_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_tau_ADHD$p.value)

GNGs_vs_GNGr_commission_ADHD <- t.test(behav_data$`Overall_Simple_Commision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Commision`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_commission_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Commision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Commision Rate`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_commission_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_commission_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_commission_ADHD$p.value)

GNGs_vs_GNGr_omission_ADHD <- t.test(behav_data$`Overall_Simple_Ommision`[behav_data$diag == 'ADHD'], behav_data$`Overall_Repeat_Ommision`[behav_data$diag == 'ADHD'], paired = TRUE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'GNGs_vs_GNGr_omission_ADHD')
mean_x_list <- c(mean_x_list, mean(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'ADHD']))
sd_x_list <- c(sd_x_list, sd(behav_data$`Overall_Simple_Ommision Rate`[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, mean(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'ADHD']))
sd_y_list <- c(sd_y_list, sd(behav_data$`Overall_Repeat_Ommision Rate`[behav_data$diag == 'ADHD']))
t.value_list <- c(t.value_list, GNGs_vs_GNGr_omission_ADHD$statistic)
df_list <- c(df_list, GNGs_vs_GNGr_omission_ADHD$parameter)
p.value_list <- c(p.value_list, GNGs_vs_GNGr_omission_ADHD$p.value)


## Difference Scores - ADHD vs. TD ##
diff_meanRT_ADHD_vs_TD <- t.test(behav_data$Overall_MeanRT_diff[behav_data$diag == 'ADHD'], behav_data$Overall_MeanRT_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_meanRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_meanRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_MeanRT_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_meanRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_MeanRT_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_meanRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_meanRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_meanRT_ADHD_vs_TD$p.value)

diff_sdRT_ADHD_vs_TD <- t.test(behav_data$Overall_SDRT_diff[behav_data$diag == 'ADHD'], behav_data$Overall_SDRT_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_sdRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_sdRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_SDRT_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_sdRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_SDRT_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_sdRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_sdRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_sdRT_ADHD_vs_TD$p.value)

diff_cvRT_ADHD_vs_TD <- t.test(behav_data$Overall_CVRT_diff[behav_data$diag == 'ADHD'], behav_data$Overall_CVRT_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_cvRT_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_cvRT_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_CVRT_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_cvRT_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_CVRT_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_cvRT_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_cvRT_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_cvRT_ADHD_vs_TD$p.value)

diff_mu_ADHD_vs_TD <- t.test(behav_data$Overall_Mu_diff[behav_data$diag == 'ADHD'], behav_data$Overall_Mu_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_mu_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_mu_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_Mu_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_mu_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_Mu_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_mu_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_mu_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_mu_ADHD_vs_TD$p.value)

diff_sigma_ADHD_vs_TD <- t.test(behav_data$Overall_Sigma_diff[behav_data$diag == 'ADHD'], behav_data$Overall_Sigma_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_sigma_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_sigma_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_Sigma_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_sigma_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_Sigma_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_sigma_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_sigma_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_sigma_ADHD_vs_TD$p.value)

diff_tau_ADHD_vs_TD <- t.test(behav_data$Overall_Tau_diff[behav_data$diag == 'ADHD'], behav_data$Overall_Tau_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_tau_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_tau_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_Tau_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_tau_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_Tau_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_tau_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_tau_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_tau_ADHD_vs_TD$p.value)

diff_commission_ADHD_vs_TD <- t.test(behav_data$Overall_Commission_diff[behav_data$diag == 'ADHD'], behav_data$Overall_Commission_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_commission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_commission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_Commission_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_commission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_Commission_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_commission_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_commission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_commission_ADHD_vs_TD$p.value)

diff_omission_ADHD_vs_TD <- t.test(behav_data$Overall_Omission_diff[behav_data$diag == 'ADHD'], behav_data$Overall_Omission_diff[behav_data$diag == 'TD'], paired = FALSE, var.equal = FALSE)
comparison_list <- c(comparison_list, 'diff_omission_ADHD_vs_TD')
mean_x_list <- c(mean_x_list, diff_omission_ADHD_vs_TD$estimate['mean of x'])
sd_x_list <- c(sd_x_list, sd(behav_data$Overall_Omission_diff[behav_data$diag == 'ADHD']))
mean_y_list <- c(mean_y_list, diff_omission_ADHD_vs_TD$estimate['mean of y'])
sd_y_list <- c(sd_y_list, sd(behav_data$Overall_Omission_diff[behav_data$diag == 'TD']))
t.value_list <- c(t.value_list, diff_omission_ADHD_vs_TD$statistic)
df_list <- c(df_list, diff_omission_ADHD_vs_TD$parameter)
p.value_list <- c(p.value_list, diff_omission_ADHD_vs_TD$p.value)


### B4. Updating the t_results_df with all the results from the t-tests and saving it as a csv file ###
t_results_df = data.frame(comparison = comparison_list, mean_x = mean_x_list, sd_x = sd_x_list, mean_y = mean_y_list, sd_y = sd_y_list, t.value = t.value_list, df = df_list, p.value = p.value_list)

write.csv(t_results_df, file = 'GNGtask_t-test_results.csv', row.names = FALSE)



#### C. Running regression analyses with diagnosis (categorical) and age (continuous) as predictors and behavioral performance as outcomes ####

#### C1. Loading the RData file containing the behavioral data ####
behav_data <- get(load('behav_data.RData'))
names(behav_data)[names(behav_data) == "Overall_Simple_Mean RT"] <- 'Overall_Simple_MeanRT'
names(behav_data)[names(behav_data) == "Overall_Simple_Std RT"] <- 'Overall_Simple_SDRT'
names(behav_data)[names(behav_data) == "Overall_Simple_CV RT"] <- 'Overall_Simple_CVRT'
names(behav_data)[names(behav_data) == "Overall_Simple_Commision Rate"] <- 'Overall_Simple_Commission'
names(behav_data)[names(behav_data) == "Overall_Simple_Ommision Rate"] <- 'Overall_Simple_Omission'

names(behav_data)[names(behav_data) == "Overall_Repeat_Mean RT"] <- 'Overall_Repeat_MeanRT'
names(behav_data)[names(behav_data) == "Overall_Repeat_Std RT"] <- 'Overall_Repeat_SDRT'
names(behav_data)[names(behav_data) == "Overall_Repeat_CV RT"] <- 'Overall_Repeat_CVRT'
names(behav_data)[names(behav_data) == "Overall_Repeat_Commision Rate"] <- 'Overall_Repeat_Commission'
names(behav_data)[names(behav_data) == "Overall_Repeat_Ommision Rate"] <- 'Overall_Repeat_Omission'


#### C2. Creating an empty data frame to save the results ####
reg_results_df <- data.frame(outcome = NA,
                             age_beta = NA,
                             age_t.value = NA,
                             age_p.value = NA,
                             diag_beta = NA,
                             diag_t.value = NA,
                             diag_p.value = NA,
                             r_squared = NA,
                             r_squared_adj = NA,
                             r_squared_f.value = NA,
                             r_squared_p.value = NA)

# Create empty lists that will be updated after each test and, at the end, will be used to update the reg_results_df before saving it as a csv file #
outcome_list <- c()
age_beta_list <- c()
age_t.value_list <- c()
age_p.value_list <- c()
diag_beta_list <- c()
diag_t.value_list <- c()
diag_p.value_list <- c()
r_squared_list <- c()
r_squared_adj_list <- c()
r_squared_f.value_list <- c()
r_squared_p.value_list <- c()

options(scipen = 999)  # set scientific notation to FALSE so that all decimal points are depicted (rather than presenting them in terms of e.g. e-05)


# Define a function that pulls out the p-value of the overall regression analysis #
lmp <- function (modelobject) {
  p <- pf(modelobject[1],modelobject[2],modelobject[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


### C3. Running the regression analyses ###
## GNGs ##
GNGs_meanRT <- summary(lm(Overall_Simple_MeanRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_meanRT")
age_beta_list <- c(age_beta_list, GNGs_meanRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_meanRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_meanRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_meanRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_meanRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_meanRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_meanRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_meanRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_meanRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_meanRT$fstatistic))

GNGs_sdRT <- summary(lm(Overall_Simple_SDRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_sdRT")
age_beta_list <- c(age_beta_list, GNGs_sdRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_sdRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_sdRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_sdRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_sdRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_sdRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_sdRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_sdRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_sdRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_sdRT$fstatistic))

GNGs_cvRT <- summary(lm(Overall_Simple_CVRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_cvRT")
age_beta_list <- c(age_beta_list, GNGs_cvRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_cvRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_cvRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_cvRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_cvRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_cvRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_cvRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_cvRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_cvRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_cvRT$fstatistic))

GNGs_mu <- summary(lm(Overall_Simple_Mu ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_mu")
age_beta_list <- c(age_beta_list, GNGs_mu$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_mu$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_mu$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_mu$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_mu$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_mu$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_mu$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_mu$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_mu$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_mu$fstatistic))

GNGs_sigma <- summary(lm(Overall_Simple_Sigma ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_sigma")
age_beta_list <- c(age_beta_list, GNGs_sigma$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_sigma$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_sigma$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_sigma$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_sigma$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_sigma$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_sigma$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_sigma$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_sigma$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_sigma$fstatistic))

GNGs_tau <- summary(lm(Overall_Simple_Tau ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_tau")
age_beta_list <- c(age_beta_list, GNGs_tau$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_tau$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_tau$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_tau$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_tau$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_tau$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_tau$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_tau$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_tau$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_tau$fstatistic))

GNGs_commission <- summary(lm(Overall_Simple_Commission ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_commission")
age_beta_list <- c(age_beta_list, GNGs_commission$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_commission$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_commission$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_commission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_commission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_commission$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_commission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_commission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_commission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_commission$fstatistic))

GNGs_omission <- summary(lm(Overall_Simple_Omission ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGs_omission")
age_beta_list <- c(age_beta_list, GNGs_omission$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGs_omission$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGs_omission$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGs_omission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGs_omission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGs_omission$coefficients[12])
r_squared_list <- c(r_squared_list, GNGs_omission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGs_omission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGs_omission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGs_omission$fstatistic))


## GNGr ##
GNGr_meanRT <- summary(lm(Overall_Repeat_MeanRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_meanRT")
age_beta_list <- c(age_beta_list, GNGr_meanRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_meanRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_meanRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_meanRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_meanRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_meanRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_meanRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_meanRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_meanRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_meanRT$fstatistic))

GNGr_sdRT <- summary(lm(Overall_Repeat_SDRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_sdRT")
age_beta_list <- c(age_beta_list, GNGr_sdRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_sdRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_sdRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_sdRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_sdRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_sdRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_sdRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_sdRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_sdRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_sdRT$fstatistic))

GNGr_cvRT <- summary(lm(Overall_Repeat_CVRT ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_cvRT")
age_beta_list <- c(age_beta_list, GNGr_cvRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_cvRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_cvRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_cvRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_cvRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_cvRT$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_cvRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_cvRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_cvRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_cvRT$fstatistic))

GNGr_mu <- summary(lm(Overall_Repeat_Mu ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_mu")
age_beta_list <- c(age_beta_list, GNGr_mu$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_mu$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_mu$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_mu$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_mu$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_mu$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_mu$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_mu$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_mu$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_mu$fstatistic))

GNGr_sigma <- summary(lm(Overall_Repeat_Sigma ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_sigma")
age_beta_list <- c(age_beta_list, GNGr_sigma$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_sigma$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_sigma$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_sigma$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_sigma$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_sigma$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_sigma$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_sigma$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_sigma$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_sigma$fstatistic))

GNGr_tau <- summary(lm(Overall_Repeat_Tau ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_tau")
age_beta_list <- c(age_beta_list, GNGr_tau$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_tau$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_tau$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_tau$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_tau$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_tau$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_tau$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_tau$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_tau$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_tau$fstatistic))

GNGr_commission <- summary(lm(Overall_Repeat_Commission ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_commission")
age_beta_list <- c(age_beta_list, GNGr_commission$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_commission$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_commission$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_commission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_commission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_commission$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_commission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_commission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_commission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_commission$fstatistic))

GNGr_omission <- summary(lm(Overall_Repeat_Omission ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "GNGr_omission")
age_beta_list <- c(age_beta_list, GNGr_omission$coefficients[2])
age_t.value_list <- c(age_t.value_list, GNGr_omission$coefficients[8])
age_p.value_list <- c(age_p.value_list, GNGr_omission$coefficients[11])
diag_beta_list <- c(diag_beta_list, GNGr_omission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, GNGr_omission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, GNGr_omission$coefficients[12])
r_squared_list <- c(r_squared_list, GNGr_omission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, GNGr_omission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, GNGr_omission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(GNGr_omission$fstatistic))


## Difference Scores ##
diff_meanRT <- summary(lm(Overall_MeanRT_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_meanRT")
age_beta_list <- c(age_beta_list, diff_meanRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_meanRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_meanRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_meanRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_meanRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_meanRT$coefficients[12])
r_squared_list <- c(r_squared_list, diff_meanRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_meanRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_meanRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_meanRT$fstatistic))

diff_sdRT <- summary(lm(Overall_SDRT_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_sdRT")
age_beta_list <- c(age_beta_list, diff_sdRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_sdRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_sdRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_sdRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_sdRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_sdRT$coefficients[12])
r_squared_list <- c(r_squared_list, diff_sdRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_sdRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_sdRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_sdRT$fstatistic))

diff_cvRT <- summary(lm(Overall_CVRT_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_cvRT")
age_beta_list <- c(age_beta_list, diff_cvRT$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_cvRT$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_cvRT$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_cvRT$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_cvRT$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_cvRT$coefficients[12])
r_squared_list <- c(r_squared_list, diff_cvRT$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_cvRT$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_cvRT$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_cvRT$fstatistic))

diff_mu <- summary(lm(Overall_Mu_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_mu")
age_beta_list <- c(age_beta_list, diff_mu$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_mu$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_mu$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_mu$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_mu$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_mu$coefficients[12])
r_squared_list <- c(r_squared_list, diff_mu$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_mu$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_mu$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_mu$fstatistic))

diff_sigma <- summary(lm(Overall_Sigma_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_sigma")
age_beta_list <- c(age_beta_list, diff_sigma$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_sigma$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_sigma$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_sigma$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_sigma$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_sigma$coefficients[12])
r_squared_list <- c(r_squared_list, diff_sigma$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_sigma$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_sigma$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_sigma$fstatistic))

diff_tau <- summary(lm(Overall_Tau_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_tau")
age_beta_list <- c(age_beta_list, diff_tau$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_tau$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_tau$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_tau$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_tau$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_tau$coefficients[12])
r_squared_list <- c(r_squared_list, diff_tau$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_tau$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_tau$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_tau$fstatistic))

diff_commission <- summary(lm(Overall_Commission_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_commission")
age_beta_list <- c(age_beta_list, diff_commission$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_commission$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_commission$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_commission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_commission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_commission$coefficients[12])
r_squared_list <- c(r_squared_list, diff_commission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_commission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_commission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_commission$fstatistic))

diff_omission <- summary(lm(Overall_Omission_diff ~ age + diag, data = behav_data))
outcome_list <- c(outcome_list, "diff_omission")
age_beta_list <- c(age_beta_list, diff_omission$coefficients[2])
age_t.value_list <- c(age_t.value_list, diff_omission$coefficients[8])
age_p.value_list <- c(age_p.value_list, diff_omission$coefficients[11])
diag_beta_list <- c(diag_beta_list, diff_omission$coefficients[3])
diag_t.value_list <- c(diag_t.value_list, diff_omission$coefficients[9])
diag_p.value_list <- c(diag_p.value_list, diff_omission$coefficients[12])
r_squared_list <- c(r_squared_list, diff_omission$r.squared)
r_squared_adj_list <- c(r_squared_adj_list, diff_omission$adj.r.squared)
r_squared_f.value_list <- c(r_squared_f.value_list, diff_omission$fstatistic[1])
r_squared_p.value_list <- c(r_squared_p.value_list, lmp(diff_omission$fstatistic))


### C4. Updating the reg_results_df with all the results from the regression analyses and saving it as a csv file ###
reg_results_df = data.frame(outcome = outcome_list, age_beta = age_beta_list, age_t.value = age_t.value_list, age_p.value = age_p.value_list, diag_beta = diag_beta_list, diag_t.value = diag_t.value_list, diag_p.value = diag_p.value_list, r_squared = r_squared_list, r_squared_adj = r_squared_adj_list, r_squared.f_value = r_squared_f.value_list, r_squared_p.value = r_squared_p.value_list)

write.csv(reg_results_df, file = 'GNGtask_regression_results.csv', row.names = FALSE)

