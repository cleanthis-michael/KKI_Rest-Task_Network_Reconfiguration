#### CM Flux Analysis ####
## February 2020, MEMitchell
## March 2020, Adapted by CM and MEM for CM Flux Analyses
## this script relies on the calculation of graph metrics from the 'Graph_Construction_KKI_CM.R' script


###### CM Flux Questions ######

### H1: higher global efficiency and lower modularity from Rest -> GNGs -> GNGr (Cohen & D'Esposito, 2016)
### H1a: lower global efficiency and higher modularity in ADHD relative to TD across task states (e.g., Lin et al., 2014) (and lower relative changes across tasks)

### H2-a1: higher PC and NDI in task-relevant networks (FPN, DAN, CON) from Rest -> GNGs -> GNGr
### H2-a2: lower PC and NDI in task-relevant networks in ADHD relative to TD

### H2-b1: lower PC and NDI in DMN from Rest -> GNGs -> GNGr (i.e., DMN will become increasingly anticorrelated from task-relevant networks with increasing cognitive demands)
### H2-b2: higher PC and and NDI in DMN in ADHD relative to TD

### H3: higher PC and NDI in sensorimotor networks (somatomotor/visual) from Rest -> GNGs -> GNGr
### H3a: higher PC and NDI in somatomotor and visual networks in ADHD relative to TD


rm(list=ls())

library(plyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ppcor)
library(psych)
library(miceadds)
#install.packages("lmerTest")
#install.packages('ppcor')
#install.packages('psych')
#install.packages('miceadds')  # to load Rdata files

library(rstudioapi)
#current_path <- getActiveDocumentContext()$path 
#setwd(dirname(current_path))
print( getwd() )



### Create a data frame that will summarize all the correlation and p-values for all comparisons (i.e., whole-brain & nodal metrics with every behavioral metric for GNGs, GNGr, and difference scores) ###
all_corr_data <- data.frame(comparison = NA,
                            threshold = NA,
                            diagnosis = NA,
                            r.value = NA,
                            p.value = NA)


#### Read in the whole-brain metrics across all thresholds ####
whole_brain_data_000 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh000/CMflux_wholebrain_metrics_04.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_005 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh005/CMflux_wholebrain_metrics_04.17.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_010 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh010/CMflux_wholebrain_metrics_04.17.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_015 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh015/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_020 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh020/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_025 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh025/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_030 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh030/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_035 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh035/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_040 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh040/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)
whole_brain_data_045 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh045/CMflux_wholebrain_metrics_03.16.2020.csv", sep = ",", stringsAsFactors = F)

diag_df <- read.csv('subject_demographics.csv', sep = ',', stringsAsFactors = FALSE)
whole_brain_data_000 <- merge(whole_brain_data_000, diag_df, by = ('sub'))
whole_brain_data_005 <- merge(whole_brain_data_005, diag_df, by = ('sub'))
whole_brain_data_010 <- merge(whole_brain_data_010, diag_df, by = ('sub'))
whole_brain_data_015 <- merge(whole_brain_data_015, diag_df, by = ('sub'))
whole_brain_data_020 <- merge(whole_brain_data_020, diag_df, by = ('sub'))
whole_brain_data_025 <- merge(whole_brain_data_025, diag_df, by = ('sub'))
whole_brain_data_030 <- merge(whole_brain_data_030, diag_df, by = ('sub'))
whole_brain_data_035 <- merge(whole_brain_data_035, diag_df, by = ('sub'))
whole_brain_data_040 <- merge(whole_brain_data_040, diag_df, by = ('sub'))
whole_brain_data_045 <- merge(whole_brain_data_045, diag_df, by = ('sub'))


#### Read in the nodal metrics within networks across all thresholds ####

# Pull density column out of the whole-brain data frame #
require(dplyr)
dat_dens_000 <- whole_brain_data_000 %>% select(sub, task, density)
dat_dens_005 <- whole_brain_data_005 %>% select(sub, task, density)
dat_dens_010 <- whole_brain_data_010 %>% select(sub, task, density)
dat_dens_015 <- whole_brain_data_015 %>% select(sub, task, density)
dat_dens_020 <- whole_brain_data_020 %>% select(sub, task, density)
dat_dens_025 <- whole_brain_data_025 %>% select(sub, task, density)
dat_dens_030 <- whole_brain_data_030 %>% select(sub, task, density)
dat_dens_035 <- whole_brain_data_035 %>% select(sub, task, density)
dat_dens_040 <- whole_brain_data_040 %>% select(sub, task, density)
dat_dens_045 <- whole_brain_data_045 %>% select(sub, task, density)


### Participation Coefficient ###
pc_netavg_000 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh000/CMflux_particcoeffsNETAVG_04.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_005 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh005/CMflux_particcoeffsNETAVG_04.17.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_010 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh010/CMflux_particcoeffsNETAVG_04.17.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_015 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh015/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_020 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh020/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_025 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh025/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_030 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh030/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_035 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh035/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_040 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh040/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg_045 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh045/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)

pc_netavg_000_Wdens = merge(pc_netavg_000, dat_dens_000, by = c("sub", "task"))
pc_netavg_000_Wdiag = merge(pc_netavg_000_Wdens, diag_df, by = ('sub'))
pc_netavg_005_Wdens = merge(pc_netavg_005, dat_dens_005, by = c("sub", "task"))
pc_netavg_005_Wdiag = merge(pc_netavg_005_Wdens, diag_df, by = ('sub'))
pc_netavg_010_Wdens = merge(pc_netavg_010, dat_dens_010, by = c("sub", "task"))
pc_netavg_010_Wdiag = merge(pc_netavg_010_Wdens, diag_df, by = ('sub'))
pc_netavg_015_Wdens = merge(pc_netavg_015, dat_dens_015, by = c("sub", "task"))
pc_netavg_015_Wdiag = merge(pc_netavg_015_Wdens, diag_df, by = ('sub'))
pc_netavg_020_Wdens = merge(pc_netavg_020, dat_dens_020, by = c("sub", "task"))
pc_netavg_020_Wdiag = merge(pc_netavg_020_Wdens, diag_df, by = ('sub'))
pc_netavg_025_Wdens = merge(pc_netavg_025, dat_dens_025, by = c("sub", "task"))
pc_netavg_025_Wdiag = merge(pc_netavg_025_Wdens, diag_df, by = ('sub'))
pc_netavg_030_Wdens = merge(pc_netavg_030, dat_dens_030, by = c("sub", "task"))
pc_netavg_030_Wdiag = merge(pc_netavg_030_Wdens, diag_df, by = ('sub'))
pc_netavg_035_Wdens = merge(pc_netavg_035, dat_dens_035, by = c("sub", "task"))
pc_netavg_035_Wdiag = merge(pc_netavg_035_Wdens, diag_df, by = ('sub'))
pc_netavg_040_Wdens = merge(pc_netavg_040, dat_dens_040, by = c("sub", "task"))
pc_netavg_040_Wdiag = merge(pc_netavg_040_Wdens, diag_df, by = ('sub'))
pc_netavg_045_Wdens = merge(pc_netavg_045, dat_dens_045, by = c("sub", "task"))
pc_netavg_045_Wdiag = merge(pc_netavg_045_Wdens, diag_df, by = ('sub'))


### Node Dissociation Index ###
ndi_netavg_000 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh000/CMflux_ndisNETAVG_04.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_005 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh005/CMflux_ndisNETAVG_04.17.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_010 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh010/CMflux_ndisNETAVG_04.17.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_015 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh015/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_020 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh020/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_025 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh025/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_030 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh030/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_035 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh035/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_040 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh040/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg_045 <- read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh045/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)

ndi_netavg_000_Wdens = merge(ndi_netavg_000, dat_dens_000, by = c("sub", "task"))
ndi_netavg_000_Wdiag = merge(ndi_netavg_000_Wdens, diag_df, by = ('sub'))
ndi_netavg_005_Wdens = merge(ndi_netavg_005, dat_dens_005, by = c("sub", "task"))
ndi_netavg_005_Wdiag = merge(ndi_netavg_005_Wdens, diag_df, by = ('sub'))
ndi_netavg_010_Wdens = merge(ndi_netavg_010, dat_dens_010, by = c("sub", "task"))
ndi_netavg_010_Wdiag = merge(ndi_netavg_010_Wdens, diag_df, by = ('sub'))
ndi_netavg_015_Wdens = merge(ndi_netavg_015, dat_dens_015, by = c("sub", "task"))
ndi_netavg_015_Wdiag = merge(ndi_netavg_015_Wdens, diag_df, by = ('sub'))
ndi_netavg_020_Wdens = merge(ndi_netavg_020, dat_dens_020, by = c("sub", "task"))
ndi_netavg_020_Wdiag = merge(ndi_netavg_020_Wdens, diag_df, by = ('sub'))
ndi_netavg_025_Wdens = merge(ndi_netavg_025, dat_dens_025, by = c("sub", "task"))
ndi_netavg_025_Wdiag = merge(ndi_netavg_025_Wdens, diag_df, by = ('sub'))
ndi_netavg_030_Wdens = merge(ndi_netavg_030, dat_dens_030, by = c("sub", "task"))
ndi_netavg_030_Wdiag = merge(ndi_netavg_030_Wdens, diag_df, by = ('sub'))
ndi_netavg_035_Wdens = merge(ndi_netavg_035, dat_dens_035, by = c("sub", "task"))
ndi_netavg_035_Wdiag = merge(ndi_netavg_035_Wdens, diag_df, by = ('sub'))
ndi_netavg_040_Wdens = merge(ndi_netavg_040, dat_dens_040, by = c("sub", "task"))
ndi_netavg_040_Wdiag = merge(ndi_netavg_040_Wdens, diag_df, by = ('sub'))
ndi_netavg_045_Wdens = merge(ndi_netavg_045, dat_dens_045, by = c("sub", "task"))
ndi_netavg_045_Wdiag = merge(ndi_netavg_045_Wdens, diag_df, by = ('sub'))


### Clustering Coefficient (AKA Transitivity) ###
clustcoeff_netavg_000 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh000/CMflux_netclustcoeff_04.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_005 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh005/CMflux_netclustcoeff_04.17.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_010 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh010/CMflux_netclustcoeff_04.17.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_015 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh015/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_020 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh020/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_025 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh025/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_030 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh030/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_035 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh035/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_040 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh040/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg_045 <-  read.csv("./thresh_testing_5CoreNetworks/graph_metrics_thresh045/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)

clustcoeff_netavg_000_Wdens = merge(clustcoeff_netavg_000, dat_dens_000, by = c("sub", "task"))
clustcoeff_netavg_000_Wdiag = merge(clustcoeff_netavg_000_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_005_Wdens = merge(clustcoeff_netavg_005, dat_dens_005, by = c("sub", "task"))
clustcoeff_netavg_005_Wdiag = merge(clustcoeff_netavg_005_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_010_Wdens = merge(clustcoeff_netavg_010, dat_dens_010, by = c("sub", "task"))
clustcoeff_netavg_010_Wdiag = merge(clustcoeff_netavg_010_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_015_Wdens = merge(clustcoeff_netavg_015, dat_dens_015, by = c("sub", "task"))
clustcoeff_netavg_015_Wdiag = merge(clustcoeff_netavg_015_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_020_Wdens = merge(clustcoeff_netavg_020, dat_dens_020, by = c("sub", "task"))
clustcoeff_netavg_020_Wdiag = merge(clustcoeff_netavg_020_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_025_Wdens = merge(clustcoeff_netavg_025, dat_dens_025, by = c("sub", "task"))
clustcoeff_netavg_025_Wdiag = merge(clustcoeff_netavg_025_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_030_Wdens = merge(clustcoeff_netavg_030, dat_dens_030, by = c("sub", "task"))
clustcoeff_netavg_030_Wdiag = merge(clustcoeff_netavg_030_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_035_Wdens = merge(clustcoeff_netavg_035, dat_dens_035, by = c("sub", "task"))
clustcoeff_netavg_035_Wdiag = merge(clustcoeff_netavg_035_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_040_Wdens = merge(clustcoeff_netavg_040, dat_dens_040, by = c("sub", "task"))
clustcoeff_netavg_040_Wdiag = merge(clustcoeff_netavg_040_Wdens, diag_df, by = ('sub'))
clustcoeff_netavg_045_Wdens = merge(clustcoeff_netavg_045, dat_dens_045, by = c("sub", "task"))
clustcoeff_netavg_045_Wdiag = merge(clustcoeff_netavg_045_Wdens, diag_df, by = ('sub'))




###### Brain (Graph Metric) - Behavior (GNG Performance) Correlations ######
library(readxl)
behav_data <- read_excel('./behavioral_data/GNG_BehavSummary.xlsx')
behav_data$sub[behav_data$sub < 1000] <- paste('0', behav_data$sub, sep = '')
behav_data$sub <- paste('sub-', behav_data$sub, sep = '')

mean(behav_data$`Overall_Simple_Mean RT`, na.rm = T)
sd(behav_data$`Overall_Simple_Mean RT`, na.rm = T)
mean(behav_data$`Overall_Simple_Std RT`, na.rm = T)
sd(behav_data$`Overall_Simple_Std RT`, na.rm = T)
mean(behav_data$`Overall_Simple_CV RT`, na.rm = T)
sd(behav_data$`Overall_Simple_CV RT`, na.rm = T)
mean(behav_data$`Overall_Simple_Mu`, na.rm = T)
sd(behav_data$`Overall_Simple_Mu`, na.rm = T)
mean(behav_data$`Overall_Simple_Sigma`, na.rm = T)
sd(behav_data$`Overall_Simple_Sigma`, na.rm = T)
mean(behav_data$`Overall_Simple_Tau`, na.rm = T)
sd(behav_data$`Overall_Simple_Tau`, na.rm = T)
mean(behav_data$`Overall_Simple_Commision Rate`, na.rm = T)
sd(behav_data$`Overall_Simple_Commision Rate`, na.rm = T)
mean(behav_data$`Overall_Simple_Ommision Rate`, na.rm = T)
sd(behav_data$`Overall_Simple_Ommision Rate`, na.rm = T)

mean(behav_data$`Overall_Repeat_CV RT`, na.rm = T)
sd(behav_data$`Overall_Repeat_CV RT`, na.rm = T)
mean(behav_data$`Overall_Repeat_Commision Rate`, na.rm = T)
sd(behav_data$`Overall_Repeat_Commision Rate`, na.rm = T)
mean(behav_data$`Overall_Repeat_Ommision Rate`, na.rm = T)
sd(behav_data$`Overall_Repeat_Ommision Rate`, na.rm = T)




### Creating the data frames that will be used for the correlation analyses ###
## GNGs ##
# Threshold 0.00 #
dat_wholebrain_GNGs_000 <- filter(whole_brain_data_000, task == 'task-GNGs')
behav_data_GNGs_thresh000 <- merge(dat_wholebrain_GNGs_000, behav_data, by = 'sub')
behav_data_GNGs_thresh000 <- merge(behav_data_GNGs_thresh000, pc_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh000 <- merge(behav_data_GNGs_thresh000, ndi_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh000 <- merge(behav_data_GNGs_thresh000, clustcoeff_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh000 <- behav_data_GNGs_thresh000 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh000, file = "brainbehav_data_GNGs_thresh000.RData")


# Threshold 0.05 #
dat_wholebrain_GNGs_005 <- filter(whole_brain_data_005, task == 'task-GNGs')
behav_data_GNGs_thresh005 <- merge(dat_wholebrain_GNGs_005, behav_data, by = 'sub')
behav_data_GNGs_thresh005 <- merge(behav_data_GNGs_thresh005, pc_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh005 <- merge(behav_data_GNGs_thresh005, ndi_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh005 <- merge(behav_data_GNGs_thresh005, clustcoeff_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh005 <- behav_data_GNGs_thresh005 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh005, file = "brainbehav_data_GNGs_thresh005.RData")


# Threshold 0.10 #
dat_wholebrain_GNGs_010 <- filter(whole_brain_data_010, task == 'task-GNGs')
behav_data_GNGs_thresh010 <- merge(dat_wholebrain_GNGs_010, behav_data, by = 'sub')
behav_data_GNGs_thresh010 <- merge(behav_data_GNGs_thresh010, pc_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh010 <- merge(behav_data_GNGs_thresh010, ndi_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh010 <- merge(behav_data_GNGs_thresh010, clustcoeff_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh010 <- behav_data_GNGs_thresh010 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh010, file = "brainbehav_data_GNGs_thresh010.RData")


# Threshold 0.15 #
dat_wholebrain_GNGs_015 <- filter(whole_brain_data_015, task == 'task-GNGs')
behav_data_GNGs_thresh015 <- merge(dat_wholebrain_GNGs_015, behav_data, by = 'sub')
behav_data_GNGs_thresh015 <- merge(behav_data_GNGs_thresh015, pc_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh015 <- merge(behav_data_GNGs_thresh015, ndi_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh015 <- merge(behav_data_GNGs_thresh015, clustcoeff_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh015 <- behav_data_GNGs_thresh015 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh015, file = "brainbehav_data_GNGs_thresh015.RData")


# Threshold 0.20 #
dat_wholebrain_GNGs_020 <- filter(whole_brain_data_020, task == 'task-GNGs')
behav_data_GNGs_thresh020 <- merge(dat_wholebrain_GNGs_020, behav_data, by = 'sub')
behav_data_GNGs_thresh020 <- merge(behav_data_GNGs_thresh020, pc_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh020 <- merge(behav_data_GNGs_thresh020, ndi_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh020 <- merge(behav_data_GNGs_thresh020, clustcoeff_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh020 <- behav_data_GNGs_thresh020 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh020, file = "brainbehav_data_GNGs_thresh020.RData")


# Threshold 0.25 #
dat_wholebrain_GNGs_025 <- filter(whole_brain_data_025, task == 'task-GNGs')
behav_data_GNGs_thresh025 <- merge(dat_wholebrain_GNGs_025, behav_data, by = 'sub')
behav_data_GNGs_thresh025 <- merge(behav_data_GNGs_thresh025, pc_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh025 <- merge(behav_data_GNGs_thresh025, ndi_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh025 <- merge(behav_data_GNGs_thresh025, clustcoeff_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh025 <- behav_data_GNGs_thresh025 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh025, file = "brainbehav_data_GNGs_thresh025.RData")


# Threshold 0.30 #
dat_wholebrain_GNGs_030 <- filter(whole_brain_data_030, task == 'task-GNGs')
behav_data_GNGs_thresh030 <- merge(dat_wholebrain_GNGs_030, behav_data, by = 'sub')
behav_data_GNGs_thresh030 <- merge(behav_data_GNGs_thresh030, pc_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh030 <- merge(behav_data_GNGs_thresh030, ndi_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh030 <- merge(behav_data_GNGs_thresh030, clustcoeff_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh030 <- behav_data_GNGs_thresh030 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh030, file = "brainbehav_data_GNGs_thresh030.RData")


# Threshold 0.35 #
dat_wholebrain_GNGs_035 <- filter(whole_brain_data_035, task == 'task-GNGs')
behav_data_GNGs_thresh035 <- merge(dat_wholebrain_GNGs_035, behav_data, by = 'sub')
behav_data_GNGs_thresh035 <- merge(behav_data_GNGs_thresh035, pc_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh035 <- merge(behav_data_GNGs_thresh035, ndi_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh035 <- merge(behav_data_GNGs_thresh035, clustcoeff_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh035 <- behav_data_GNGs_thresh035 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh035, file = "brainbehav_data_GNGs_thresh035.RData")


# Threshold 0.40 #
dat_wholebrain_GNGs_040 <- filter(whole_brain_data_040, task == 'task-GNGs')
behav_data_GNGs_thresh040 <- merge(dat_wholebrain_GNGs_040, behav_data, by = 'sub')
behav_data_GNGs_thresh040 <- merge(behav_data_GNGs_thresh040, pc_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh040 <- merge(behav_data_GNGs_thresh040, ndi_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh040 <- merge(behav_data_GNGs_thresh040, clustcoeff_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh040 <- behav_data_GNGs_thresh040 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh040, file = "brainbehav_data_GNGs_thresh040.RData")


# Threshold 0.45 #
dat_wholebrain_GNGs_045 <- filter(whole_brain_data_045, task == 'task-GNGs')
behav_data_GNGs_thresh045 <- merge(dat_wholebrain_GNGs_045, behav_data, by = 'sub')
behav_data_GNGs_thresh045 <- merge(behav_data_GNGs_thresh045, pc_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh045 <- merge(behav_data_GNGs_thresh045, ndi_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGs_thresh045 <- merge(behav_data_GNGs_thresh045, clustcoeff_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGs_thresh045 <- behav_data_GNGs_thresh045 %>% dplyr::select(sub:ge, 'Overall_Simple_Mean RT':'Overall_Simple_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGs_thresh045, file = "brainbehav_data_GNGs_thresh045.RData")


## GNGr ##
# Threshold 0.00 #
dat_wholebrain_GNGr_000 <- filter(whole_brain_data_000, task == 'task-GNGr')
behav_data_GNGr_thresh000 <- merge(dat_wholebrain_GNGr_000, behav_data, by = 'sub')
behav_data_GNGr_thresh000 <- merge(behav_data_GNGr_thresh000, pc_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh000 <- merge(behav_data_GNGr_thresh000, ndi_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh000 <- merge(behav_data_GNGr_thresh000, clustcoeff_netavg_000_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh000 <- behav_data_GNGr_thresh000 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh000, file = "brainbehav_data_GNGr_thresh000.RData")


# Threshold 0.05 #
dat_wholebrain_GNGr_005 <- filter(whole_brain_data_005, task == 'task-GNGr')
behav_data_GNGr_thresh005 <- merge(dat_wholebrain_GNGr_005, behav_data, by = 'sub')
behav_data_GNGr_thresh005 <- merge(behav_data_GNGr_thresh005, pc_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh005 <- merge(behav_data_GNGr_thresh005, ndi_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh005 <- merge(behav_data_GNGr_thresh005, clustcoeff_netavg_005_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh005 <- behav_data_GNGr_thresh005 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh005, file = "brainbehav_data_GNGr_thresh005.RData")


# Threshold 0.10 #
dat_wholebrain_GNGr_010 <- filter(whole_brain_data_010, task == 'task-GNGr')
behav_data_GNGr_thresh010 <- merge(dat_wholebrain_GNGr_010, behav_data, by = 'sub')
behav_data_GNGr_thresh010 <- merge(behav_data_GNGr_thresh010, pc_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh010 <- merge(behav_data_GNGr_thresh010, ndi_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh010 <- merge(behav_data_GNGr_thresh010, clustcoeff_netavg_010_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh010 <- behav_data_GNGr_thresh010 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh010, file = "brainbehav_data_GNGr_thresh010.RData")


# Threshold 0.15 #
dat_wholebrain_GNGr_015 <- filter(whole_brain_data_015, task == 'task-GNGr')
behav_data_GNGr_thresh015 <- merge(dat_wholebrain_GNGr_015, behav_data, by = 'sub')
behav_data_GNGr_thresh015 <- merge(behav_data_GNGr_thresh015, pc_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh015 <- merge(behav_data_GNGr_thresh015, ndi_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh015 <- merge(behav_data_GNGr_thresh015, clustcoeff_netavg_015_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh015 <- behav_data_GNGr_thresh015 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh015, file = "brainbehav_data_GNGr_thresh015.RData")


# Threshold 0.20 #
dat_wholebrain_GNGr_020 <- filter(whole_brain_data_020, task == 'task-GNGr')
behav_data_GNGr_thresh020 <- merge(dat_wholebrain_GNGr_020, behav_data, by = 'sub')
behav_data_GNGr_thresh020 <- merge(behav_data_GNGr_thresh020, pc_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh020 <- merge(behav_data_GNGr_thresh020, ndi_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh020 <- merge(behav_data_GNGr_thresh020, clustcoeff_netavg_020_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh020 <- behav_data_GNGr_thresh020 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh020, file = "brainbehav_data_GNGr_thresh020.RData")


# Threshold 0.25 #
dat_wholebrain_GNGr_025 <- filter(whole_brain_data_025, task == 'task-GNGr')
behav_data_GNGr_thresh025 <- merge(dat_wholebrain_GNGr_025, behav_data, by = 'sub')
behav_data_GNGr_thresh025 <- merge(behav_data_GNGr_thresh025, pc_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh025 <- merge(behav_data_GNGr_thresh025, ndi_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh025 <- merge(behav_data_GNGr_thresh025, clustcoeff_netavg_025_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh025 <- behav_data_GNGr_thresh025 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh025, file = "brainbehav_data_GNGr_thresh025.RData")


# Threshold 0.30 #
dat_wholebrain_GNGr_030 <- filter(whole_brain_data_030, task == 'task-GNGr')
behav_data_GNGr_thresh030 <- merge(dat_wholebrain_GNGr_030, behav_data, by = 'sub')
behav_data_GNGr_thresh030 <- merge(behav_data_GNGr_thresh030, pc_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh030 <- merge(behav_data_GNGr_thresh030, ndi_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh030 <- merge(behav_data_GNGr_thresh030, clustcoeff_netavg_030_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh030 <- behav_data_GNGr_thresh030 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh030, file = "brainbehav_data_GNGr_thresh030.RData")


# Threshold 0.35 #
dat_wholebrain_GNGr_035 <- filter(whole_brain_data_035, task == 'task-GNGr')
behav_data_GNGr_thresh035 <- merge(dat_wholebrain_GNGr_035, behav_data, by = 'sub')
behav_data_GNGr_thresh035 <- merge(behav_data_GNGr_thresh035, pc_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh035 <- merge(behav_data_GNGr_thresh035, ndi_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh035 <- merge(behav_data_GNGr_thresh035, clustcoeff_netavg_035_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh035 <- behav_data_GNGr_thresh035 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh035, file = "brainbehav_data_GNGr_thresh035.RData")


# Threshold 0.40 #
dat_wholebrain_GNGr_040 <- filter(whole_brain_data_040, task == 'task-GNGr')
behav_data_GNGr_thresh040 <- merge(dat_wholebrain_GNGr_040, behav_data, by = 'sub')
behav_data_GNGr_thresh040 <- merge(behav_data_GNGr_thresh040, pc_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh040 <- merge(behav_data_GNGr_thresh040, ndi_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh040 <- merge(behav_data_GNGr_thresh040, clustcoeff_netavg_040_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh040 <- behav_data_GNGr_thresh040 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh040, file = "brainbehav_data_GNGr_thresh040.RData")


# Threshold 0.45 #
dat_wholebrain_GNGr_045 <- filter(whole_brain_data_045, task == 'task-GNGr')
behav_data_GNGr_thresh045 <- merge(dat_wholebrain_GNGr_045, behav_data, by = 'sub')
behav_data_GNGr_thresh045 <- merge(behav_data_GNGr_thresh045, pc_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh045 <- merge(behav_data_GNGr_thresh045, ndi_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender'))
behav_data_GNGr_thresh045 <- merge(behav_data_GNGr_thresh045, clustcoeff_netavg_045_Wdiag, by = c('sub', 'task', 'density', 'diag', 'age', 'gender', 'ses'))
brainbehav_data_GNGr_thresh045 <- behav_data_GNGr_thresh045 %>% dplyr::select(sub:ge, 'Overall_Repeat_Mean RT':'Overall_Repeat_Ommision Rate', avg_pc_net1:clustcoeff_visual)

save(brainbehav_data_GNGr_thresh045, file = "brainbehav_data_GNGr_thresh045.RData")


### All Data for Difference Score Analyses ###
dat_allGNG_thresh000 <- merge(behav_data_GNGs_thresh000, behav_data_GNGr_thresh000, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh005 <- merge(behav_data_GNGs_thresh005, behav_data_GNGr_thresh005, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh010 <- merge(behav_data_GNGs_thresh010, behav_data_GNGr_thresh010, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh015 <- merge(behav_data_GNGs_thresh015, behav_data_GNGr_thresh015, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh020 <- merge(behav_data_GNGs_thresh020, behav_data_GNGr_thresh020, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh025 <- merge(behav_data_GNGs_thresh025, behav_data_GNGr_thresh025, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh030 <- merge(behav_data_GNGs_thresh030, behav_data_GNGr_thresh030, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh035 <- merge(behav_data_GNGs_thresh035, behav_data_GNGr_thresh035, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh040 <- merge(behav_data_GNGs_thresh040, behav_data_GNGr_thresh040, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))

dat_allGNG_thresh045 <- merge(behav_data_GNGs_thresh045, behav_data_GNGr_thresh045, by = c('sub', 'diag', 'age', 'gender', 'ses', 'Overall_Simple_Mean RT', 'Overall_Simple_Std RT', 'Overall_Simple_CV RT', 'Overall_Simple_Mu', 'Overall_Simple_Sigma', 'Overall_Simple_Tau', 'Overall_Simple_Commision Rate', 'Overall_Simple_Ommision Rate', 'Overall_Repeat_Mean RT', 'Overall_Repeat_Std RT', 'Overall_Repeat_CV RT', 'Overall_Repeat_Mu', 'Overall_Repeat_Sigma', 'Overall_Repeat_Tau', 'Overall_Repeat_Commision Rate', 'Overall_Repeat_Ommision Rate'))


## Calculating Difference Scores - GNGr - GNGs ##
## These analyses are to look at whether changes in brain network organization from GNGs to GNGr are related to changes in behavioral performance from GNGs to GNGr ##
# Threshold 0.00 #
dat_allGNG_thresh000$mod_diff <- (dat_allGNG_thresh000$mod.y - dat_allGNG_thresh000$mod.x)
dat_allGNG_thresh000$ge_diff <- (dat_allGNG_thresh000$ge.y - dat_allGNG_thresh000$ge.x)
dat_allGNG_thresh000$density_diff <- (dat_allGNG_thresh000$density.y - dat_allGNG_thresh000$density.x)

dat_allGNG_thresh000$avg_pc_net1_diff <- (dat_allGNG_thresh000$avg_pc_net1.y - dat_allGNG_thresh000$avg_pc_net1.x)
dat_allGNG_thresh000$avg_pc_net2_diff <- (dat_allGNG_thresh000$avg_pc_net2.y - dat_allGNG_thresh000$avg_pc_net2.x)
dat_allGNG_thresh000$avg_pc_net3_diff <- (dat_allGNG_thresh000$avg_pc_net3.y - dat_allGNG_thresh000$avg_pc_net3.x)
dat_allGNG_thresh000$avg_pc_net4_diff <- (dat_allGNG_thresh000$avg_pc_net4.y - dat_allGNG_thresh000$avg_pc_net4.x)
dat_allGNG_thresh000$avg_pc_net5_diff <- (dat_allGNG_thresh000$avg_pc_net5.y - dat_allGNG_thresh000$avg_pc_net5.x)

dat_allGNG_thresh000$avg_ndi_net1_diff <- (dat_allGNG_thresh000$avg_ndi_net1.y - dat_allGNG_thresh000$avg_ndi_net1.x)
dat_allGNG_thresh000$avg_ndi_net2_diff <- (dat_allGNG_thresh000$avg_ndi_net2.y - dat_allGNG_thresh000$avg_ndi_net2.x)
dat_allGNG_thresh000$avg_ndi_net3_diff <- (dat_allGNG_thresh000$avg_ndi_net3.y - dat_allGNG_thresh000$avg_ndi_net3.x)
dat_allGNG_thresh000$avg_ndi_net4_diff <- (dat_allGNG_thresh000$avg_ndi_net4.y - dat_allGNG_thresh000$avg_ndi_net4.x)
dat_allGNG_thresh000$avg_ndi_net5_diff <- (dat_allGNG_thresh000$avg_ndi_net5.y - dat_allGNG_thresh000$avg_ndi_net5.x)

dat_allGNG_thresh000$clustcoeff_fpn_diff <- (dat_allGNG_thresh000$clustcoeff_fpn.y - dat_allGNG_thresh000$clustcoeff_fpn.x)
dat_allGNG_thresh000$clustcoeff_con_diff <- (dat_allGNG_thresh000$clustcoeff_con.y - dat_allGNG_thresh000$clustcoeff_con.x)
dat_allGNG_thresh000$clustcoeff_dmn_diff <- (dat_allGNG_thresh000$clustcoeff_dmn.y - dat_allGNG_thresh000$clustcoeff_dmn.x)
dat_allGNG_thresh000$clustcoeff_smd_diff <- (dat_allGNG_thresh000$clustcoeff_smd.y - dat_allGNG_thresh000$clustcoeff_smd.x)
dat_allGNG_thresh000$clustcoeff_visual_diff <- (dat_allGNG_thresh000$clustcoeff_visual.y - dat_allGNG_thresh000$clustcoeff_visual.x)

dat_allGNG_thresh000$`Overall_MeanRT_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Mean RT' - dat_allGNG_thresh000$'Overall_Simple_Mean RT'
dat_allGNG_thresh000$`Overall_SDRT_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Std RT' - dat_allGNG_thresh000$'Overall_Simple_Std RT'
dat_allGNG_thresh000$`Overall_CVRT_diff` <- dat_allGNG_thresh000$'Overall_Repeat_CV RT' - dat_allGNG_thresh000$'Overall_Simple_CV RT'
dat_allGNG_thresh000$`Overall_SDRT_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Std RT' - dat_allGNG_thresh000$'Overall_Simple_Std RT'
dat_allGNG_thresh000$`Overall_Mu_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Mu' - dat_allGNG_thresh000$'Overall_Simple_Mu'
dat_allGNG_thresh000$`Overall_Sigma_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Sigma' - dat_allGNG_thresh000$'Overall_Simple_Sigma'
dat_allGNG_thresh000$`Overall_Tau_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Tau' - dat_allGNG_thresh000$'Overall_Simple_Tau'
dat_allGNG_thresh000$`Overall_Commission_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh000$'Overall_Simple_Commision Rate'
dat_allGNG_thresh000$`Overall_Omission_diff` <- dat_allGNG_thresh000$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh000$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh000 <- dat_allGNG_thresh000 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh000, file = "brainbehav_data_diff_thresh000.RData")


# Threshold 0.05 #
dat_allGNG_thresh005$mod_diff <- (dat_allGNG_thresh005$mod.y - dat_allGNG_thresh005$mod.x)
dat_allGNG_thresh005$ge_diff <- (dat_allGNG_thresh005$ge.y - dat_allGNG_thresh005$ge.x)
dat_allGNG_thresh005$density_diff <- (dat_allGNG_thresh005$density.y - dat_allGNG_thresh005$density.x)

dat_allGNG_thresh005$avg_pc_net1_diff <- (dat_allGNG_thresh005$avg_pc_net1.y - dat_allGNG_thresh005$avg_pc_net1.x)
dat_allGNG_thresh005$avg_pc_net2_diff <- (dat_allGNG_thresh005$avg_pc_net2.y - dat_allGNG_thresh005$avg_pc_net2.x)
dat_allGNG_thresh005$avg_pc_net3_diff <- (dat_allGNG_thresh005$avg_pc_net3.y - dat_allGNG_thresh005$avg_pc_net3.x)
dat_allGNG_thresh005$avg_pc_net4_diff <- (dat_allGNG_thresh005$avg_pc_net4.y - dat_allGNG_thresh005$avg_pc_net4.x)
dat_allGNG_thresh005$avg_pc_net5_diff <- (dat_allGNG_thresh005$avg_pc_net5.y - dat_allGNG_thresh005$avg_pc_net5.x)

dat_allGNG_thresh005$avg_ndi_net1_diff <- (dat_allGNG_thresh005$avg_ndi_net1.y - dat_allGNG_thresh005$avg_ndi_net1.x)
dat_allGNG_thresh005$avg_ndi_net2_diff <- (dat_allGNG_thresh005$avg_ndi_net2.y - dat_allGNG_thresh005$avg_ndi_net2.x)
dat_allGNG_thresh005$avg_ndi_net3_diff <- (dat_allGNG_thresh005$avg_ndi_net3.y - dat_allGNG_thresh005$avg_ndi_net3.x)
dat_allGNG_thresh005$avg_ndi_net4_diff <- (dat_allGNG_thresh005$avg_ndi_net4.y - dat_allGNG_thresh005$avg_ndi_net4.x)
dat_allGNG_thresh005$avg_ndi_net5_diff <- (dat_allGNG_thresh005$avg_ndi_net5.y - dat_allGNG_thresh005$avg_ndi_net5.x)

dat_allGNG_thresh005$clustcoeff_fpn_diff <- (dat_allGNG_thresh005$clustcoeff_fpn.y - dat_allGNG_thresh005$clustcoeff_fpn.x)
dat_allGNG_thresh005$clustcoeff_con_diff <- (dat_allGNG_thresh005$clustcoeff_con.y - dat_allGNG_thresh005$clustcoeff_con.x)
dat_allGNG_thresh005$clustcoeff_dmn_diff <- (dat_allGNG_thresh005$clustcoeff_dmn.y - dat_allGNG_thresh005$clustcoeff_dmn.x)
dat_allGNG_thresh005$clustcoeff_smd_diff <- (dat_allGNG_thresh005$clustcoeff_smd.y - dat_allGNG_thresh005$clustcoeff_smd.x)
dat_allGNG_thresh005$clustcoeff_visual_diff <- (dat_allGNG_thresh005$clustcoeff_visual.y - dat_allGNG_thresh005$clustcoeff_visual.x)

dat_allGNG_thresh005$`Overall_MeanRT_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Mean RT' - dat_allGNG_thresh005$'Overall_Simple_Mean RT'
dat_allGNG_thresh005$`Overall_SDRT_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Std RT' - dat_allGNG_thresh005$'Overall_Simple_Std RT'
dat_allGNG_thresh005$`Overall_CVRT_diff` <- dat_allGNG_thresh005$'Overall_Repeat_CV RT' - dat_allGNG_thresh005$'Overall_Simple_CV RT'
dat_allGNG_thresh005$`Overall_SDRT_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Std RT' - dat_allGNG_thresh005$'Overall_Simple_Std RT'
dat_allGNG_thresh005$`Overall_Mu_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Mu' - dat_allGNG_thresh005$'Overall_Simple_Mu'
dat_allGNG_thresh005$`Overall_Sigma_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Sigma' - dat_allGNG_thresh005$'Overall_Simple_Sigma'
dat_allGNG_thresh005$`Overall_Tau_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Tau' - dat_allGNG_thresh005$'Overall_Simple_Tau'
dat_allGNG_thresh005$`Overall_Commission_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh005$'Overall_Simple_Commision Rate'
dat_allGNG_thresh005$`Overall_Omission_diff` <- dat_allGNG_thresh005$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh005$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh005 <- dat_allGNG_thresh005 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh005, file = "brainbehav_data_diff_thresh005.RData")


# Threshold 0.10 #
dat_allGNG_thresh010$mod_diff <- (dat_allGNG_thresh010$mod.y - dat_allGNG_thresh010$mod.x)
dat_allGNG_thresh010$ge_diff <- (dat_allGNG_thresh010$ge.y - dat_allGNG_thresh010$ge.x)
dat_allGNG_thresh010$density_diff <- (dat_allGNG_thresh010$density.y - dat_allGNG_thresh010$density.x)

dat_allGNG_thresh010$avg_pc_net1_diff <- (dat_allGNG_thresh010$avg_pc_net1.y - dat_allGNG_thresh010$avg_pc_net1.x)
dat_allGNG_thresh010$avg_pc_net2_diff <- (dat_allGNG_thresh010$avg_pc_net2.y - dat_allGNG_thresh010$avg_pc_net2.x)
dat_allGNG_thresh010$avg_pc_net3_diff <- (dat_allGNG_thresh010$avg_pc_net3.y - dat_allGNG_thresh010$avg_pc_net3.x)
dat_allGNG_thresh010$avg_pc_net4_diff <- (dat_allGNG_thresh010$avg_pc_net4.y - dat_allGNG_thresh010$avg_pc_net4.x)
dat_allGNG_thresh010$avg_pc_net5_diff <- (dat_allGNG_thresh010$avg_pc_net5.y - dat_allGNG_thresh010$avg_pc_net5.x)

dat_allGNG_thresh010$avg_ndi_net1_diff <- (dat_allGNG_thresh010$avg_ndi_net1.y - dat_allGNG_thresh010$avg_ndi_net1.x)
dat_allGNG_thresh010$avg_ndi_net2_diff <- (dat_allGNG_thresh010$avg_ndi_net2.y - dat_allGNG_thresh010$avg_ndi_net2.x)
dat_allGNG_thresh010$avg_ndi_net3_diff <- (dat_allGNG_thresh010$avg_ndi_net3.y - dat_allGNG_thresh010$avg_ndi_net3.x)
dat_allGNG_thresh010$avg_ndi_net4_diff <- (dat_allGNG_thresh010$avg_ndi_net4.y - dat_allGNG_thresh010$avg_ndi_net4.x)
dat_allGNG_thresh010$avg_ndi_net5_diff <- (dat_allGNG_thresh010$avg_ndi_net5.y - dat_allGNG_thresh010$avg_ndi_net5.x)

dat_allGNG_thresh010$clustcoeff_fpn_diff <- (dat_allGNG_thresh010$clustcoeff_fpn.y - dat_allGNG_thresh010$clustcoeff_fpn.x)
dat_allGNG_thresh010$clustcoeff_con_diff <- (dat_allGNG_thresh010$clustcoeff_con.y - dat_allGNG_thresh010$clustcoeff_con.x)
dat_allGNG_thresh010$clustcoeff_dmn_diff <- (dat_allGNG_thresh010$clustcoeff_dmn.y - dat_allGNG_thresh010$clustcoeff_dmn.x)
dat_allGNG_thresh010$clustcoeff_smd_diff <- (dat_allGNG_thresh010$clustcoeff_smd.y - dat_allGNG_thresh010$clustcoeff_smd.x)
dat_allGNG_thresh010$clustcoeff_visual_diff <- (dat_allGNG_thresh010$clustcoeff_visual.y - dat_allGNG_thresh010$clustcoeff_visual.x)

dat_allGNG_thresh010$`Overall_MeanRT_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Mean RT' - dat_allGNG_thresh010$'Overall_Simple_Mean RT'
dat_allGNG_thresh010$`Overall_SDRT_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Std RT' - dat_allGNG_thresh010$'Overall_Simple_Std RT'
dat_allGNG_thresh010$`Overall_CVRT_diff` <- dat_allGNG_thresh010$'Overall_Repeat_CV RT' - dat_allGNG_thresh010$'Overall_Simple_CV RT'
dat_allGNG_thresh010$`Overall_SDRT_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Std RT' - dat_allGNG_thresh010$'Overall_Simple_Std RT'
dat_allGNG_thresh010$`Overall_Mu_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Mu' - dat_allGNG_thresh010$'Overall_Simple_Mu'
dat_allGNG_thresh010$`Overall_Sigma_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Sigma' - dat_allGNG_thresh010$'Overall_Simple_Sigma'
dat_allGNG_thresh010$`Overall_Tau_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Tau' - dat_allGNG_thresh010$'Overall_Simple_Tau'
dat_allGNG_thresh010$`Overall_Commission_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh010$'Overall_Simple_Commision Rate'
dat_allGNG_thresh010$`Overall_Omission_diff` <- dat_allGNG_thresh010$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh010$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh010 <- dat_allGNG_thresh010 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh010, file = "brainbehav_data_diff_thresh010.RData")


# Threshold 0.15 #
dat_allGNG_thresh015$mod_diff <- (dat_allGNG_thresh015$mod.y - dat_allGNG_thresh015$mod.x)
dat_allGNG_thresh015$ge_diff <- (dat_allGNG_thresh015$ge.y - dat_allGNG_thresh015$ge.x)
dat_allGNG_thresh015$density_diff <- (dat_allGNG_thresh015$density.y - dat_allGNG_thresh015$density.x)

dat_allGNG_thresh015$avg_pc_net1_diff <- (dat_allGNG_thresh015$avg_pc_net1.y - dat_allGNG_thresh015$avg_pc_net1.x)
dat_allGNG_thresh015$avg_pc_net2_diff <- (dat_allGNG_thresh015$avg_pc_net2.y - dat_allGNG_thresh015$avg_pc_net2.x)
dat_allGNG_thresh015$avg_pc_net3_diff <- (dat_allGNG_thresh015$avg_pc_net3.y - dat_allGNG_thresh015$avg_pc_net3.x)
dat_allGNG_thresh015$avg_pc_net4_diff <- (dat_allGNG_thresh015$avg_pc_net4.y - dat_allGNG_thresh015$avg_pc_net4.x)
dat_allGNG_thresh015$avg_pc_net5_diff <- (dat_allGNG_thresh015$avg_pc_net5.y - dat_allGNG_thresh015$avg_pc_net5.x)

dat_allGNG_thresh015$avg_ndi_net1_diff <- (dat_allGNG_thresh015$avg_ndi_net1.y - dat_allGNG_thresh015$avg_ndi_net1.x)
dat_allGNG_thresh015$avg_ndi_net2_diff <- (dat_allGNG_thresh015$avg_ndi_net2.y - dat_allGNG_thresh015$avg_ndi_net2.x)
dat_allGNG_thresh015$avg_ndi_net3_diff <- (dat_allGNG_thresh015$avg_ndi_net3.y - dat_allGNG_thresh015$avg_ndi_net3.x)
dat_allGNG_thresh015$avg_ndi_net4_diff <- (dat_allGNG_thresh015$avg_ndi_net4.y - dat_allGNG_thresh015$avg_ndi_net4.x)
dat_allGNG_thresh015$avg_ndi_net5_diff <- (dat_allGNG_thresh015$avg_ndi_net5.y - dat_allGNG_thresh015$avg_ndi_net5.x)

dat_allGNG_thresh015$clustcoeff_fpn_diff <- (dat_allGNG_thresh015$clustcoeff_fpn.y - dat_allGNG_thresh015$clustcoeff_fpn.x)
dat_allGNG_thresh015$clustcoeff_con_diff <- (dat_allGNG_thresh015$clustcoeff_con.y - dat_allGNG_thresh015$clustcoeff_con.x)
dat_allGNG_thresh015$clustcoeff_dmn_diff <- (dat_allGNG_thresh015$clustcoeff_dmn.y - dat_allGNG_thresh015$clustcoeff_dmn.x)
dat_allGNG_thresh015$clustcoeff_smd_diff <- (dat_allGNG_thresh015$clustcoeff_smd.y - dat_allGNG_thresh015$clustcoeff_smd.x)
dat_allGNG_thresh015$clustcoeff_visual_diff <- (dat_allGNG_thresh015$clustcoeff_visual.y - dat_allGNG_thresh015$clustcoeff_visual.x)

dat_allGNG_thresh015$`Overall_MeanRT_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Mean RT' - dat_allGNG_thresh015$'Overall_Simple_Mean RT'
dat_allGNG_thresh015$`Overall_SDRT_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Std RT' - dat_allGNG_thresh015$'Overall_Simple_Std RT'
dat_allGNG_thresh015$`Overall_CVRT_diff` <- dat_allGNG_thresh015$'Overall_Repeat_CV RT' - dat_allGNG_thresh015$'Overall_Simple_CV RT'
dat_allGNG_thresh015$`Overall_SDRT_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Std RT' - dat_allGNG_thresh015$'Overall_Simple_Std RT'
dat_allGNG_thresh015$`Overall_Mu_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Mu' - dat_allGNG_thresh015$'Overall_Simple_Mu'
dat_allGNG_thresh015$`Overall_Sigma_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Sigma' - dat_allGNG_thresh015$'Overall_Simple_Sigma'
dat_allGNG_thresh015$`Overall_Tau_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Tau' - dat_allGNG_thresh015$'Overall_Simple_Tau'
dat_allGNG_thresh015$`Overall_Commission_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh015$'Overall_Simple_Commision Rate'
dat_allGNG_thresh015$`Overall_Omission_diff` <- dat_allGNG_thresh015$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh015$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh015 <- dat_allGNG_thresh015 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh015, file = "brainbehav_data_diff_thresh015.RData")


# Threshold 0.20 #
dat_allGNG_thresh020$mod_diff <- (dat_allGNG_thresh020$mod.y - dat_allGNG_thresh020$mod.x)
dat_allGNG_thresh020$ge_diff <- (dat_allGNG_thresh020$ge.y - dat_allGNG_thresh020$ge.x)
dat_allGNG_thresh020$density_diff <- (dat_allGNG_thresh020$density.y - dat_allGNG_thresh020$density.x)

dat_allGNG_thresh020$avg_pc_net1_diff <- (dat_allGNG_thresh020$avg_pc_net1.y - dat_allGNG_thresh020$avg_pc_net1.x)
dat_allGNG_thresh020$avg_pc_net2_diff <- (dat_allGNG_thresh020$avg_pc_net2.y - dat_allGNG_thresh020$avg_pc_net2.x)
dat_allGNG_thresh020$avg_pc_net3_diff <- (dat_allGNG_thresh020$avg_pc_net3.y - dat_allGNG_thresh020$avg_pc_net3.x)
dat_allGNG_thresh020$avg_pc_net4_diff <- (dat_allGNG_thresh020$avg_pc_net4.y - dat_allGNG_thresh020$avg_pc_net4.x)
dat_allGNG_thresh020$avg_pc_net5_diff <- (dat_allGNG_thresh020$avg_pc_net5.y - dat_allGNG_thresh020$avg_pc_net5.x)

dat_allGNG_thresh020$avg_ndi_net1_diff <- (dat_allGNG_thresh020$avg_ndi_net1.y - dat_allGNG_thresh020$avg_ndi_net1.x)
dat_allGNG_thresh020$avg_ndi_net2_diff <- (dat_allGNG_thresh020$avg_ndi_net2.y - dat_allGNG_thresh020$avg_ndi_net2.x)
dat_allGNG_thresh020$avg_ndi_net3_diff <- (dat_allGNG_thresh020$avg_ndi_net3.y - dat_allGNG_thresh020$avg_ndi_net3.x)
dat_allGNG_thresh020$avg_ndi_net4_diff <- (dat_allGNG_thresh020$avg_ndi_net4.y - dat_allGNG_thresh020$avg_ndi_net4.x)
dat_allGNG_thresh020$avg_ndi_net5_diff <- (dat_allGNG_thresh020$avg_ndi_net5.y - dat_allGNG_thresh020$avg_ndi_net5.x)

dat_allGNG_thresh020$clustcoeff_fpn_diff <- (dat_allGNG_thresh020$clustcoeff_fpn.y - dat_allGNG_thresh020$clustcoeff_fpn.x)
dat_allGNG_thresh020$clustcoeff_con_diff <- (dat_allGNG_thresh020$clustcoeff_con.y - dat_allGNG_thresh020$clustcoeff_con.x)
dat_allGNG_thresh020$clustcoeff_dmn_diff <- (dat_allGNG_thresh020$clustcoeff_dmn.y - dat_allGNG_thresh020$clustcoeff_dmn.x)
dat_allGNG_thresh020$clustcoeff_smd_diff <- (dat_allGNG_thresh020$clustcoeff_smd.y - dat_allGNG_thresh020$clustcoeff_smd.x)
dat_allGNG_thresh020$clustcoeff_visual_diff <- (dat_allGNG_thresh020$clustcoeff_visual.y - dat_allGNG_thresh020$clustcoeff_visual.x)

dat_allGNG_thresh020$`Overall_MeanRT_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Mean RT' - dat_allGNG_thresh020$'Overall_Simple_Mean RT'
dat_allGNG_thresh020$`Overall_SDRT_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Std RT' - dat_allGNG_thresh020$'Overall_Simple_Std RT'
dat_allGNG_thresh020$`Overall_CVRT_diff` <- dat_allGNG_thresh020$'Overall_Repeat_CV RT' - dat_allGNG_thresh020$'Overall_Simple_CV RT'
dat_allGNG_thresh020$`Overall_SDRT_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Std RT' - dat_allGNG_thresh020$'Overall_Simple_Std RT'
dat_allGNG_thresh020$`Overall_Mu_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Mu' - dat_allGNG_thresh020$'Overall_Simple_Mu'
dat_allGNG_thresh020$`Overall_Sigma_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Sigma' - dat_allGNG_thresh020$'Overall_Simple_Sigma'
dat_allGNG_thresh020$`Overall_Tau_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Tau' - dat_allGNG_thresh020$'Overall_Simple_Tau'
dat_allGNG_thresh020$`Overall_Commission_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh020$'Overall_Simple_Commision Rate'
dat_allGNG_thresh020$`Overall_Omission_diff` <- dat_allGNG_thresh020$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh020$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh020 <- dat_allGNG_thresh020 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh020, file = "brainbehav_data_diff_thresh020.RData")


# Threshold 0.25 #
dat_allGNG_thresh025$mod_diff <- (dat_allGNG_thresh025$mod.y - dat_allGNG_thresh025$mod.x)
dat_allGNG_thresh025$ge_diff <- (dat_allGNG_thresh025$ge.y - dat_allGNG_thresh025$ge.x)
dat_allGNG_thresh025$density_diff <- (dat_allGNG_thresh025$density.y - dat_allGNG_thresh025$density.x)

dat_allGNG_thresh025$avg_pc_net1_diff <- (dat_allGNG_thresh025$avg_pc_net1.y - dat_allGNG_thresh025$avg_pc_net1.x)
dat_allGNG_thresh025$avg_pc_net2_diff <- (dat_allGNG_thresh025$avg_pc_net2.y - dat_allGNG_thresh025$avg_pc_net2.x)
dat_allGNG_thresh025$avg_pc_net3_diff <- (dat_allGNG_thresh025$avg_pc_net3.y - dat_allGNG_thresh025$avg_pc_net3.x)
dat_allGNG_thresh025$avg_pc_net4_diff <- (dat_allGNG_thresh025$avg_pc_net4.y - dat_allGNG_thresh025$avg_pc_net4.x)
dat_allGNG_thresh025$avg_pc_net5_diff <- (dat_allGNG_thresh025$avg_pc_net5.y - dat_allGNG_thresh025$avg_pc_net5.x)

dat_allGNG_thresh025$avg_ndi_net1_diff <- (dat_allGNG_thresh025$avg_ndi_net1.y - dat_allGNG_thresh025$avg_ndi_net1.x)
dat_allGNG_thresh025$avg_ndi_net2_diff <- (dat_allGNG_thresh025$avg_ndi_net2.y - dat_allGNG_thresh025$avg_ndi_net2.x)
dat_allGNG_thresh025$avg_ndi_net3_diff <- (dat_allGNG_thresh025$avg_ndi_net3.y - dat_allGNG_thresh025$avg_ndi_net3.x)
dat_allGNG_thresh025$avg_ndi_net4_diff <- (dat_allGNG_thresh025$avg_ndi_net4.y - dat_allGNG_thresh025$avg_ndi_net4.x)
dat_allGNG_thresh025$avg_ndi_net5_diff <- (dat_allGNG_thresh025$avg_ndi_net5.y - dat_allGNG_thresh025$avg_ndi_net5.x)

dat_allGNG_thresh025$clustcoeff_fpn_diff <- (dat_allGNG_thresh025$clustcoeff_fpn.y - dat_allGNG_thresh025$clustcoeff_fpn.x)
dat_allGNG_thresh025$clustcoeff_con_diff <- (dat_allGNG_thresh025$clustcoeff_con.y - dat_allGNG_thresh025$clustcoeff_con.x)
dat_allGNG_thresh025$clustcoeff_dmn_diff <- (dat_allGNG_thresh025$clustcoeff_dmn.y - dat_allGNG_thresh025$clustcoeff_dmn.x)
dat_allGNG_thresh025$clustcoeff_smd_diff <- (dat_allGNG_thresh025$clustcoeff_smd.y - dat_allGNG_thresh025$clustcoeff_smd.x)
dat_allGNG_thresh025$clustcoeff_visual_diff <- (dat_allGNG_thresh025$clustcoeff_visual.y - dat_allGNG_thresh025$clustcoeff_visual.x)

dat_allGNG_thresh025$`Overall_MeanRT_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Mean RT' - dat_allGNG_thresh025$'Overall_Simple_Mean RT'
dat_allGNG_thresh025$`Overall_SDRT_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Std RT' - dat_allGNG_thresh025$'Overall_Simple_Std RT'
dat_allGNG_thresh025$`Overall_CVRT_diff` <- dat_allGNG_thresh025$'Overall_Repeat_CV RT' - dat_allGNG_thresh025$'Overall_Simple_CV RT'
dat_allGNG_thresh025$`Overall_SDRT_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Std RT' - dat_allGNG_thresh025$'Overall_Simple_Std RT'
dat_allGNG_thresh025$`Overall_Mu_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Mu' - dat_allGNG_thresh025$'Overall_Simple_Mu'
dat_allGNG_thresh025$`Overall_Sigma_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Sigma' - dat_allGNG_thresh025$'Overall_Simple_Sigma'
dat_allGNG_thresh025$`Overall_Tau_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Tau' - dat_allGNG_thresh025$'Overall_Simple_Tau'
dat_allGNG_thresh025$`Overall_Commission_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh025$'Overall_Simple_Commision Rate'
dat_allGNG_thresh025$`Overall_Omission_diff` <- dat_allGNG_thresh025$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh025$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh025 <- dat_allGNG_thresh025 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh025, file = "brainbehav_data_diff_thresh025.RData")


# Threshold 0.30 #
dat_allGNG_thresh030$mod_diff <- (dat_allGNG_thresh030$mod.y - dat_allGNG_thresh030$mod.x)
dat_allGNG_thresh030$ge_diff <- (dat_allGNG_thresh030$ge.y - dat_allGNG_thresh030$ge.x)
dat_allGNG_thresh030$density_diff <- (dat_allGNG_thresh030$density.y - dat_allGNG_thresh030$density.x)

dat_allGNG_thresh030$avg_pc_net1_diff <- (dat_allGNG_thresh030$avg_pc_net1.y - dat_allGNG_thresh030$avg_pc_net1.x)
dat_allGNG_thresh030$avg_pc_net2_diff <- (dat_allGNG_thresh030$avg_pc_net2.y - dat_allGNG_thresh030$avg_pc_net2.x)
dat_allGNG_thresh030$avg_pc_net3_diff <- (dat_allGNG_thresh030$avg_pc_net3.y - dat_allGNG_thresh030$avg_pc_net3.x)
dat_allGNG_thresh030$avg_pc_net4_diff <- (dat_allGNG_thresh030$avg_pc_net4.y - dat_allGNG_thresh030$avg_pc_net4.x)
dat_allGNG_thresh030$avg_pc_net5_diff <- (dat_allGNG_thresh030$avg_pc_net5.y - dat_allGNG_thresh030$avg_pc_net5.x)

dat_allGNG_thresh030$avg_ndi_net1_diff <- (dat_allGNG_thresh030$avg_ndi_net1.y - dat_allGNG_thresh030$avg_ndi_net1.x)
dat_allGNG_thresh030$avg_ndi_net2_diff <- (dat_allGNG_thresh030$avg_ndi_net2.y - dat_allGNG_thresh030$avg_ndi_net2.x)
dat_allGNG_thresh030$avg_ndi_net3_diff <- (dat_allGNG_thresh030$avg_ndi_net3.y - dat_allGNG_thresh030$avg_ndi_net3.x)
dat_allGNG_thresh030$avg_ndi_net4_diff <- (dat_allGNG_thresh030$avg_ndi_net4.y - dat_allGNG_thresh030$avg_ndi_net4.x)
dat_allGNG_thresh030$avg_ndi_net5_diff <- (dat_allGNG_thresh030$avg_ndi_net5.y - dat_allGNG_thresh030$avg_ndi_net5.x)

dat_allGNG_thresh030$clustcoeff_fpn_diff <- (dat_allGNG_thresh030$clustcoeff_fpn.y - dat_allGNG_thresh030$clustcoeff_fpn.x)
dat_allGNG_thresh030$clustcoeff_con_diff <- (dat_allGNG_thresh030$clustcoeff_con.y - dat_allGNG_thresh030$clustcoeff_con.x)
dat_allGNG_thresh030$clustcoeff_dmn_diff <- (dat_allGNG_thresh030$clustcoeff_dmn.y - dat_allGNG_thresh030$clustcoeff_dmn.x)
dat_allGNG_thresh030$clustcoeff_smd_diff <- (dat_allGNG_thresh030$clustcoeff_smd.y - dat_allGNG_thresh030$clustcoeff_smd.x)
dat_allGNG_thresh030$clustcoeff_visual_diff <- (dat_allGNG_thresh030$clustcoeff_visual.y - dat_allGNG_thresh030$clustcoeff_visual.x)

dat_allGNG_thresh030$`Overall_MeanRT_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Mean RT' - dat_allGNG_thresh030$'Overall_Simple_Mean RT'
dat_allGNG_thresh030$`Overall_SDRT_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Std RT' - dat_allGNG_thresh030$'Overall_Simple_Std RT'
dat_allGNG_thresh030$`Overall_CVRT_diff` <- dat_allGNG_thresh030$'Overall_Repeat_CV RT' - dat_allGNG_thresh030$'Overall_Simple_CV RT'
dat_allGNG_thresh030$`Overall_SDRT_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Std RT' - dat_allGNG_thresh030$'Overall_Simple_Std RT'
dat_allGNG_thresh030$`Overall_Mu_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Mu' - dat_allGNG_thresh030$'Overall_Simple_Mu'
dat_allGNG_thresh030$`Overall_Sigma_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Sigma' - dat_allGNG_thresh030$'Overall_Simple_Sigma'
dat_allGNG_thresh030$`Overall_Tau_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Tau' - dat_allGNG_thresh030$'Overall_Simple_Tau'
dat_allGNG_thresh030$`Overall_Commission_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh030$'Overall_Simple_Commision Rate'
dat_allGNG_thresh030$`Overall_Omission_diff` <- dat_allGNG_thresh030$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh030$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh030 <- dat_allGNG_thresh030 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh030, file = "brainbehav_data_diff_thresh030.RData")


# Threshold 0.35 #
dat_allGNG_thresh035$mod_diff <- (dat_allGNG_thresh035$mod.y - dat_allGNG_thresh035$mod.x)
dat_allGNG_thresh035$ge_diff <- (dat_allGNG_thresh035$ge.y - dat_allGNG_thresh035$ge.x)
dat_allGNG_thresh035$density_diff <- (dat_allGNG_thresh035$density.y - dat_allGNG_thresh035$density.x)

dat_allGNG_thresh035$avg_pc_net1_diff <- (dat_allGNG_thresh035$avg_pc_net1.y - dat_allGNG_thresh035$avg_pc_net1.x)
dat_allGNG_thresh035$avg_pc_net2_diff <- (dat_allGNG_thresh035$avg_pc_net2.y - dat_allGNG_thresh035$avg_pc_net2.x)
dat_allGNG_thresh035$avg_pc_net3_diff <- (dat_allGNG_thresh035$avg_pc_net3.y - dat_allGNG_thresh035$avg_pc_net3.x)
dat_allGNG_thresh035$avg_pc_net4_diff <- (dat_allGNG_thresh035$avg_pc_net4.y - dat_allGNG_thresh035$avg_pc_net4.x)
dat_allGNG_thresh035$avg_pc_net5_diff <- (dat_allGNG_thresh035$avg_pc_net5.y - dat_allGNG_thresh035$avg_pc_net5.x)

dat_allGNG_thresh035$avg_ndi_net1_diff <- (dat_allGNG_thresh035$avg_ndi_net1.y - dat_allGNG_thresh035$avg_ndi_net1.x)
dat_allGNG_thresh035$avg_ndi_net2_diff <- (dat_allGNG_thresh035$avg_ndi_net2.y - dat_allGNG_thresh035$avg_ndi_net2.x)
dat_allGNG_thresh035$avg_ndi_net3_diff <- (dat_allGNG_thresh035$avg_ndi_net3.y - dat_allGNG_thresh035$avg_ndi_net3.x)
dat_allGNG_thresh035$avg_ndi_net4_diff <- (dat_allGNG_thresh035$avg_ndi_net4.y - dat_allGNG_thresh035$avg_ndi_net4.x)
dat_allGNG_thresh035$avg_ndi_net5_diff <- (dat_allGNG_thresh035$avg_ndi_net5.y - dat_allGNG_thresh035$avg_ndi_net5.x)

dat_allGNG_thresh035$clustcoeff_fpn_diff <- (dat_allGNG_thresh035$clustcoeff_fpn.y - dat_allGNG_thresh035$clustcoeff_fpn.x)
dat_allGNG_thresh035$clustcoeff_con_diff <- (dat_allGNG_thresh035$clustcoeff_con.y - dat_allGNG_thresh035$clustcoeff_con.x)
dat_allGNG_thresh035$clustcoeff_dmn_diff <- (dat_allGNG_thresh035$clustcoeff_dmn.y - dat_allGNG_thresh035$clustcoeff_dmn.x)
dat_allGNG_thresh035$clustcoeff_smd_diff <- (dat_allGNG_thresh035$clustcoeff_smd.y - dat_allGNG_thresh035$clustcoeff_smd.x)
dat_allGNG_thresh035$clustcoeff_visual_diff <- (dat_allGNG_thresh035$clustcoeff_visual.y - dat_allGNG_thresh035$clustcoeff_visual.x)

dat_allGNG_thresh035$`Overall_MeanRT_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Mean RT' - dat_allGNG_thresh035$'Overall_Simple_Mean RT'
dat_allGNG_thresh035$`Overall_SDRT_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Std RT' - dat_allGNG_thresh035$'Overall_Simple_Std RT'
dat_allGNG_thresh035$`Overall_CVRT_diff` <- dat_allGNG_thresh035$'Overall_Repeat_CV RT' - dat_allGNG_thresh035$'Overall_Simple_CV RT'
dat_allGNG_thresh035$`Overall_SDRT_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Std RT' - dat_allGNG_thresh035$'Overall_Simple_Std RT'
dat_allGNG_thresh035$`Overall_Mu_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Mu' - dat_allGNG_thresh035$'Overall_Simple_Mu'
dat_allGNG_thresh035$`Overall_Sigma_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Sigma' - dat_allGNG_thresh035$'Overall_Simple_Sigma'
dat_allGNG_thresh035$`Overall_Tau_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Tau' - dat_allGNG_thresh035$'Overall_Simple_Tau'
dat_allGNG_thresh035$`Overall_Commission_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh035$'Overall_Simple_Commision Rate'
dat_allGNG_thresh035$`Overall_Omission_diff` <- dat_allGNG_thresh035$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh035$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh035 <- dat_allGNG_thresh035 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh035, file = "brainbehav_data_diff_thresh035.RData")


# Threshold 0.40 #
dat_allGNG_thresh040$mod_diff <- (dat_allGNG_thresh040$mod.y - dat_allGNG_thresh040$mod.x)
dat_allGNG_thresh040$ge_diff <- (dat_allGNG_thresh040$ge.y - dat_allGNG_thresh040$ge.x)
dat_allGNG_thresh040$density_diff <- (dat_allGNG_thresh040$density.y - dat_allGNG_thresh040$density.x)

dat_allGNG_thresh040$avg_pc_net1_diff <- (dat_allGNG_thresh040$avg_pc_net1.y - dat_allGNG_thresh040$avg_pc_net1.x)
dat_allGNG_thresh040$avg_pc_net2_diff <- (dat_allGNG_thresh040$avg_pc_net2.y - dat_allGNG_thresh040$avg_pc_net2.x)
dat_allGNG_thresh040$avg_pc_net3_diff <- (dat_allGNG_thresh040$avg_pc_net3.y - dat_allGNG_thresh040$avg_pc_net3.x)
dat_allGNG_thresh040$avg_pc_net4_diff <- (dat_allGNG_thresh040$avg_pc_net4.y - dat_allGNG_thresh040$avg_pc_net4.x)
dat_allGNG_thresh040$avg_pc_net5_diff <- (dat_allGNG_thresh040$avg_pc_net5.y - dat_allGNG_thresh040$avg_pc_net5.x)

dat_allGNG_thresh040$avg_ndi_net1_diff <- (dat_allGNG_thresh040$avg_ndi_net1.y - dat_allGNG_thresh040$avg_ndi_net1.x)
dat_allGNG_thresh040$avg_ndi_net2_diff <- (dat_allGNG_thresh040$avg_ndi_net2.y - dat_allGNG_thresh040$avg_ndi_net2.x)
dat_allGNG_thresh040$avg_ndi_net3_diff <- (dat_allGNG_thresh040$avg_ndi_net3.y - dat_allGNG_thresh040$avg_ndi_net3.x)
dat_allGNG_thresh040$avg_ndi_net4_diff <- (dat_allGNG_thresh040$avg_ndi_net4.y - dat_allGNG_thresh040$avg_ndi_net4.x)
dat_allGNG_thresh040$avg_ndi_net5_diff <- (dat_allGNG_thresh040$avg_ndi_net5.y - dat_allGNG_thresh040$avg_ndi_net5.x)

dat_allGNG_thresh040$clustcoeff_fpn_diff <- (dat_allGNG_thresh040$clustcoeff_fpn.y - dat_allGNG_thresh040$clustcoeff_fpn.x)
dat_allGNG_thresh040$clustcoeff_con_diff <- (dat_allGNG_thresh040$clustcoeff_con.y - dat_allGNG_thresh040$clustcoeff_con.x)
dat_allGNG_thresh040$clustcoeff_dmn_diff <- (dat_allGNG_thresh040$clustcoeff_dmn.y - dat_allGNG_thresh040$clustcoeff_dmn.x)
dat_allGNG_thresh040$clustcoeff_smd_diff <- (dat_allGNG_thresh040$clustcoeff_smd.y - dat_allGNG_thresh040$clustcoeff_smd.x)
dat_allGNG_thresh040$clustcoeff_visual_diff <- (dat_allGNG_thresh040$clustcoeff_visual.y - dat_allGNG_thresh040$clustcoeff_visual.x)

dat_allGNG_thresh040$`Overall_MeanRT_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Mean RT' - dat_allGNG_thresh040$'Overall_Simple_Mean RT'
dat_allGNG_thresh040$`Overall_SDRT_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Std RT' - dat_allGNG_thresh040$'Overall_Simple_Std RT'
dat_allGNG_thresh040$`Overall_CVRT_diff` <- dat_allGNG_thresh040$'Overall_Repeat_CV RT' - dat_allGNG_thresh040$'Overall_Simple_CV RT'
dat_allGNG_thresh040$`Overall_SDRT_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Std RT' - dat_allGNG_thresh040$'Overall_Simple_Std RT'
dat_allGNG_thresh040$`Overall_Mu_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Mu' - dat_allGNG_thresh040$'Overall_Simple_Mu'
dat_allGNG_thresh040$`Overall_Sigma_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Sigma' - dat_allGNG_thresh040$'Overall_Simple_Sigma'
dat_allGNG_thresh040$`Overall_Tau_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Tau' - dat_allGNG_thresh040$'Overall_Simple_Tau'
dat_allGNG_thresh040$`Overall_Commission_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh040$'Overall_Simple_Commision Rate'
dat_allGNG_thresh040$`Overall_Omission_diff` <- dat_allGNG_thresh040$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh040$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh040 <- dat_allGNG_thresh040 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh040, file = "brainbehav_data_diff_thresh040.RData")


# Threshold 0.45 #
dat_allGNG_thresh045$mod_diff <- (dat_allGNG_thresh045$mod.y - dat_allGNG_thresh045$mod.x)
dat_allGNG_thresh045$ge_diff <- (dat_allGNG_thresh045$ge.y - dat_allGNG_thresh045$ge.x)
dat_allGNG_thresh045$density_diff <- (dat_allGNG_thresh045$density.y - dat_allGNG_thresh045$density.x)

dat_allGNG_thresh045$avg_pc_net1_diff <- (dat_allGNG_thresh045$avg_pc_net1.y - dat_allGNG_thresh045$avg_pc_net1.x)
dat_allGNG_thresh045$avg_pc_net2_diff <- (dat_allGNG_thresh045$avg_pc_net2.y - dat_allGNG_thresh045$avg_pc_net2.x)
dat_allGNG_thresh045$avg_pc_net3_diff <- (dat_allGNG_thresh045$avg_pc_net3.y - dat_allGNG_thresh045$avg_pc_net3.x)
dat_allGNG_thresh045$avg_pc_net4_diff <- (dat_allGNG_thresh045$avg_pc_net4.y - dat_allGNG_thresh045$avg_pc_net4.x)
dat_allGNG_thresh045$avg_pc_net5_diff <- (dat_allGNG_thresh045$avg_pc_net5.y - dat_allGNG_thresh045$avg_pc_net5.x)

dat_allGNG_thresh045$avg_ndi_net1_diff <- (dat_allGNG_thresh045$avg_ndi_net1.y - dat_allGNG_thresh045$avg_ndi_net1.x)
dat_allGNG_thresh045$avg_ndi_net2_diff <- (dat_allGNG_thresh045$avg_ndi_net2.y - dat_allGNG_thresh045$avg_ndi_net2.x)
dat_allGNG_thresh045$avg_ndi_net3_diff <- (dat_allGNG_thresh045$avg_ndi_net3.y - dat_allGNG_thresh045$avg_ndi_net3.x)
dat_allGNG_thresh045$avg_ndi_net4_diff <- (dat_allGNG_thresh045$avg_ndi_net4.y - dat_allGNG_thresh045$avg_ndi_net4.x)
dat_allGNG_thresh045$avg_ndi_net5_diff <- (dat_allGNG_thresh045$avg_ndi_net5.y - dat_allGNG_thresh045$avg_ndi_net5.x)

dat_allGNG_thresh045$clustcoeff_fpn_diff <- (dat_allGNG_thresh045$clustcoeff_fpn.y - dat_allGNG_thresh045$clustcoeff_fpn.x)
dat_allGNG_thresh045$clustcoeff_con_diff <- (dat_allGNG_thresh045$clustcoeff_con.y - dat_allGNG_thresh045$clustcoeff_con.x)
dat_allGNG_thresh045$clustcoeff_dmn_diff <- (dat_allGNG_thresh045$clustcoeff_dmn.y - dat_allGNG_thresh045$clustcoeff_dmn.x)
dat_allGNG_thresh045$clustcoeff_smd_diff <- (dat_allGNG_thresh045$clustcoeff_smd.y - dat_allGNG_thresh045$clustcoeff_smd.x)
dat_allGNG_thresh045$clustcoeff_visual_diff <- (dat_allGNG_thresh045$clustcoeff_visual.y - dat_allGNG_thresh045$clustcoeff_visual.x)

dat_allGNG_thresh045$`Overall_MeanRT_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Mean RT' - dat_allGNG_thresh045$'Overall_Simple_Mean RT'
dat_allGNG_thresh045$`Overall_SDRT_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Std RT' - dat_allGNG_thresh045$'Overall_Simple_Std RT'
dat_allGNG_thresh045$`Overall_CVRT_diff` <- dat_allGNG_thresh045$'Overall_Repeat_CV RT' - dat_allGNG_thresh045$'Overall_Simple_CV RT'
dat_allGNG_thresh045$`Overall_SDRT_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Std RT' - dat_allGNG_thresh045$'Overall_Simple_Std RT'
dat_allGNG_thresh045$`Overall_Mu_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Mu' - dat_allGNG_thresh045$'Overall_Simple_Mu'
dat_allGNG_thresh045$`Overall_Sigma_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Sigma' - dat_allGNG_thresh045$'Overall_Simple_Sigma'
dat_allGNG_thresh045$`Overall_Tau_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Tau' - dat_allGNG_thresh045$'Overall_Simple_Tau'
dat_allGNG_thresh045$`Overall_Commission_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Commision Rate' - dat_allGNG_thresh045$'Overall_Simple_Commision Rate'
dat_allGNG_thresh045$`Overall_Omission_diff` <- dat_allGNG_thresh045$'Overall_Repeat_Ommision Rate' - dat_allGNG_thresh045$'Overall_Simple_Ommision Rate'

brainbehav_data_diff_thresh045 <- dat_allGNG_thresh045 %>% dplyr::select(sub:gender, mod_diff:Overall_Omission_diff)

save(brainbehav_data_diff_thresh045, file = "brainbehav_data_diff_thresh045.RData")



#### Loop through the GNGs RData files generated above, compute partial brain x behavior correlations controlling for density, and note them down in the all_corr_data data frame ####

data_list_gngs <- c("./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh000.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh005.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh010.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh015.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh020.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh025.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh030.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh035.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh040.RData", "./brainxbehav_data_GNGs/brainbehav_data_GNGs_thresh045.RData")


# This loop looks at GNGs brain network organization and GNGs task performance #
for (i in 1:length(data_list_gngs)) {
  
  df_temp <- get(load(data_list_gngs[i]))
  temp_corr_data <- data.frame(comparison = NA,
                               threshold = NA,
                               diagnosis = NA,
                               r.value = NA,
                               p.value = NA)
  
  thresh = sapply(data_list_gngs[i], function (s) strsplit(strsplit(s, "_")[[1]][6], ".R")[[1]][1])  
  
  
  # Make 5 lists, 1 for each variable of the data frame
  comparison_list = c()
  threshold_list = c()
  diagnosis_list = c()
  r_list = c()
  p_list = c()
  
  
  #### Run each correlation of interest and append the relevant values to each list ####
  ### Whole-Brain Metrics ###
  ## Modularity ##
  GNGs_mod_meanRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_mod_meanRT_TD$p.value)
  
  GNGs_mod_meanRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_meanRT_ADHD$p.value)
  
  GNGs_mod_SDRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_mod_SDRT_TD$p.value)
  
  GNGs_mod_SDRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_SDRT_ADHD$p.value)
  
  GNGs_mod_CVRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_mod_CVRT_TD$p.value)
  
  GNGs_mod_CVRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_CVRT_ADHD$p.value)
  
  GNGs_mod_mu_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_mu_TD$estimate)
  p_list <- c(p_list, GNGs_mod_mu_TD$p.value)
  
  GNGs_mod_mu_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_mu_ADHD$p.value)
  
  GNGs_mod_sigma_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_mod_sigma_TD$p.value)
  
  GNGs_mod_sigma_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_sigma_ADHD$p.value)
  
  GNGs_mod_tau_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_tau_TD$estimate)
  p_list <- c(p_list, GNGs_mod_tau_TD$p.value)
  
  GNGs_mod_tau_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_tau_ADHD$p.value)
  
  GNGs_mod_commission_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_commission_TD$estimate)
  p_list <- c(p_list, GNGs_mod_commission_TD$p.value)
  
  GNGs_mod_commission_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_commission_ADHD$p.value)
  
  GNGs_mod_omission_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_mod_omission_TD$estimate)
  p_list <- c(p_list, GNGs_mod_omission_TD$p.value)
  
  GNGs_mod_omission_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_mod_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_mod_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_mod_omission_ADHD$p.value)
  
  
  ## Global Efficiency ##
  GNGs_ge_meanRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ge_meanRT_TD$p.value)
  
  GNGs_ge_meanRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_meanRT_ADHD$p.value)
  
  GNGs_ge_SDRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ge_SDRT_TD$p.value)
  
  GNGs_ge_SDRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_SDRT_ADHD$p.value)
  
  GNGs_ge_CVRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ge_CVRT_TD$p.value)
  
  GNGs_ge_CVRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_CVRT_ADHD$p.value)
  
  GNGs_ge_mu_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ge_mu_TD$p.value)
  
  GNGs_ge_mu_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_mu_ADHD$p.value)
  
  GNGs_ge_sigma_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ge_sigma_TD$p.value)
  
  GNGs_ge_sigma_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_sigma_ADHD$p.value)
  
  GNGs_ge_tau_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ge_tau_TD$p.value)
  
  GNGs_ge_tau_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_tau_ADHD$p.value)
  
  GNGs_ge_commission_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ge_commission_TD$p.value)
  
  GNGs_ge_commission_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_commission_ADHD$p.value)
  
  GNGs_ge_omission_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ge_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ge_omission_TD$p.value)
  
  GNGs_ge_omission_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGs_ge_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ge_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ge_omission_ADHD$p.value)
  
  
  ### Nodal Metrics ###
  ## Participation Coefficient ##
  # Frontoparietal Network #
  GNGs_pc_fpn_meanRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_meanRT_TD$p.value)
  
  GNGs_pc_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_meanRT_ADHD$p.value)
  
  GNGs_pc_fpn_SDRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_SDRT_TD$p.value)
  
  GNGs_pc_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_SDRT_ADHD$p.value)
  
  GNGs_pc_fpn_CVRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_CVRT_TD$p.value)
  
  GNGs_pc_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_CVRT_ADHD$p.value)
  
  GNGs_pc_fpn_mu_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_mu_TD$p.value)
  
  GNGs_pc_fpn_mu_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_mu_ADHD$p.value)
  
  GNGs_pc_fpn_sigma_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_sigma_TD$p.value)
  
  GNGs_pc_fpn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_sigma_ADHD$p.value)
  
  GNGs_pc_fpn_tau_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_tau_TD$p.value)
  
  GNGs_pc_fpn_tau_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_tau_ADHD$p.value)
  
  GNGs_pc_fpn_commission_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_commission_TD$p.value)
  
  GNGs_pc_fpn_commission_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_commission_ADHD$p.value)
  
  GNGs_pc_fpn_omission_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_omission_TD$p.value)
  
  GNGs_pc_fpn_omission_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGs_pc_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGs_pc_con_meanRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_meanRT_TD$p.value)
  
  GNGs_pc_con_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_meanRT_ADHD$p.value)
  
  GNGs_pc_con_SDRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_SDRT_TD$p.value)
  
  GNGs_pc_con_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_SDRT_ADHD$p.value)
  
  GNGs_pc_con_CVRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_CVRT_TD$p.value)
  
  GNGs_pc_con_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_CVRT_ADHD$p.value)
  
  GNGs_pc_con_mu_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_mu_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_mu_TD$p.value)
  
  GNGs_pc_con_mu_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_mu_ADHD$p.value)
  
  GNGs_pc_con_sigma_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_sigma_TD$p.value)
  
  GNGs_pc_con_sigma_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_sigma_ADHD$p.value)
  
  GNGs_pc_con_tau_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_tau_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_tau_TD$p.value)
  
  GNGs_pc_con_tau_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_tau_ADHD$p.value)
  
  GNGs_pc_con_commission_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_commission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_commission_TD$p.value)
  
  GNGs_pc_con_commission_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_commission_ADHD$p.value)
  
  GNGs_pc_con_omission_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_con_omission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_con_omission_TD$p.value)
  
  GNGs_pc_con_omission_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGs_pc_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGs_pc_dmn_meanRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_meanRT_TD$p.value)
  
  GNGs_pc_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_meanRT_ADHD$p.value)
  
  GNGs_pc_dmn_SDRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_SDRT_TD$p.value)
  
  GNGs_pc_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_SDRT_ADHD$p.value)
  
  GNGs_pc_dmn_CVRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_CVRT_TD$p.value)
  
  GNGs_pc_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_CVRT_ADHD$p.value)
  
  GNGs_pc_dmn_mu_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_mu_TD$p.value)
  
  GNGs_pc_dmn_mu_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_mu_ADHD$p.value)
  
  GNGs_pc_dmn_sigma_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_sigma_TD$p.value)
  
  GNGs_pc_dmn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_sigma_ADHD$p.value)
  
  GNGs_pc_dmn_tau_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_tau_TD$p.value)
  
  GNGs_pc_dmn_tau_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_tau_ADHD$p.value)
  
  GNGs_pc_dmn_commission_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_commission_TD$p.value)
  
  GNGs_pc_dmn_commission_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_commission_ADHD$p.value)
  
  GNGs_pc_dmn_omission_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_omission_TD$p.value)
  
  GNGs_pc_dmn_omission_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGs_pc_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGs_pc_smd_meanRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_meanRT_TD$p.value)
  
  GNGs_pc_smd_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_meanRT_ADHD$p.value)
  
  GNGs_pc_smd_SDRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_SDRT_TD$p.value)
  
  GNGs_pc_smd_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_SDRT_ADHD$p.value)
  
  GNGs_pc_smd_CVRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_CVRT_TD$p.value)
  
  GNGs_pc_smd_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_CVRT_ADHD$p.value)
  
  GNGs_pc_smd_mu_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_mu_TD$p.value)
  
  GNGs_pc_smd_mu_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_mu_ADHD$p.value)
  
  GNGs_pc_smd_sigma_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_sigma_TD$p.value)
  
  GNGs_pc_smd_sigma_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_sigma_ADHD$p.value)
  
  GNGs_pc_smd_tau_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_tau_TD$p.value)
  
  GNGs_pc_smd_tau_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_tau_ADHD$p.value)
  
  GNGs_pc_smd_commission_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_commission_TD$p.value)
  
  GNGs_pc_smd_commission_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_commission_ADHD$p.value)
  
  GNGs_pc_smd_omission_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_omission_TD$p.value)
  
  GNGs_pc_smd_omission_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGs_pc_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGs_pc_visual_meanRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_meanRT_TD$p.value)
  
  GNGs_pc_visual_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_meanRT_ADHD$p.value)
  
  GNGs_pc_visual_SDRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_SDRT_TD$p.value)
  
  GNGs_pc_visual_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_SDRT_ADHD$p.value)
  
  GNGs_pc_visual_CVRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_CVRT_TD$p.value)
  
  GNGs_pc_visual_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_CVRT_ADHD$p.value)
  
  GNGs_pc_visual_mu_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_mu_TD$p.value)
  
  GNGs_pc_visual_mu_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_mu_ADHD$p.value)
  
  GNGs_pc_visual_sigma_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_sigma_TD$p.value)
  
  GNGs_pc_visual_sigma_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_sigma_ADHD$p.value)
  
  GNGs_pc_visual_tau_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_tau_TD$p.value)
  
  GNGs_pc_visual_tau_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_tau_ADHD$p.value)
  
  GNGs_pc_visual_commission_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_commission_TD$p.value)
  
  GNGs_pc_visual_commission_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_commission_ADHD$p.value)
  
  GNGs_pc_visual_omission_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_pc_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_omission_TD$p.value)
  
  GNGs_pc_visual_omission_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGs_pc_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_pc_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_pc_visual_omission_ADHD$p.value)
  
  
  ## Node Dissociation Index ##
  # Frontoparietal Network #
  GNGs_ndi_fpn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_meanRT_TD$p.value)
  
  GNGs_ndi_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_meanRT_ADHD$p.value)
  
  GNGs_ndi_fpn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_SDRT_TD$p.value)
  
  GNGs_ndi_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_SDRT_ADHD$p.value)
  
  GNGs_ndi_fpn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_CVRT_TD$p.value)
  
  GNGs_ndi_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_CVRT_ADHD$p.value)
  
  GNGs_ndi_fpn_mu_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_mu_TD$p.value)
  
  GNGs_ndi_fpn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_mu_ADHD$p.value)
  
  GNGs_ndi_fpn_sigma_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_sigma_TD$p.value)
  
  GNGs_ndi_fpn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_sigma_ADHD$p.value)
  
  GNGs_ndi_fpn_tau_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_tau_TD$p.value)
  
  GNGs_ndi_fpn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_tau_ADHD$p.value)
  
  GNGs_ndi_fpn_commission_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_commission_TD$p.value)
  
  GNGs_ndi_fpn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_commission_ADHD$p.value)
  
  GNGs_ndi_fpn_omission_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_omission_TD$p.value)
  
  GNGs_ndi_fpn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGs_ndi_con_meanRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_meanRT_TD$p.value)
  
  GNGs_ndi_con_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_meanRT_ADHD$p.value)
  
  GNGs_ndi_con_SDRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_SDRT_TD$p.value)
  
  GNGs_ndi_con_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_SDRT_ADHD$p.value)
  
  GNGs_ndi_con_CVRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_CVRT_TD$p.value)
  
  GNGs_ndi_con_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_CVRT_ADHD$p.value)
  
  GNGs_ndi_con_mu_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_mu_TD$p.value)
  
  GNGs_ndi_con_mu_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_mu_ADHD$p.value)
  
  GNGs_ndi_con_sigma_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_sigma_TD$p.value)
  
  GNGs_ndi_con_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_sigma_ADHD$p.value)
  
  GNGs_ndi_con_tau_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_tau_TD$p.value)
  
  GNGs_ndi_con_tau_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_tau_ADHD$p.value)
  
  GNGs_ndi_con_commission_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_commission_TD$p.value)
  
  GNGs_ndi_con_commission_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_commission_ADHD$p.value)
  
  GNGs_ndi_con_omission_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_con_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_omission_TD$p.value)
  
  GNGs_ndi_con_omission_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGs_ndi_dmn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_meanRT_TD$p.value)
  
  GNGs_ndi_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_meanRT_ADHD$p.value)
  
  GNGs_ndi_dmn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_SDRT_TD$p.value)
  
  GNGs_ndi_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_SDRT_ADHD$p.value)
  
  GNGs_ndi_dmn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_CVRT_TD$p.value)
  
  GNGs_ndi_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_CVRT_ADHD$p.value)
  
  GNGs_ndi_dmn_mu_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_mu_TD$p.value)
  
  GNGs_ndi_dmn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_mu_ADHD$p.value)
  
  GNGs_ndi_dmn_sigma_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_sigma_TD$p.value)
  
  GNGs_ndi_dmn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_sigma_ADHD$p.value)
  
  GNGs_ndi_dmn_tau_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_tau_TD$p.value)
  
  GNGs_ndi_dmn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_tau_ADHD$p.value)
  
  GNGs_ndi_dmn_commission_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_commission_TD$p.value)
  
  GNGs_ndi_dmn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_commission_ADHD$p.value)
  
  GNGs_ndi_dmn_omission_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_omission_TD$p.value)
  
  GNGs_ndi_dmn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGs_ndi_smd_meanRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_meanRT_TD$p.value)
  
  GNGs_ndi_smd_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_meanRT_ADHD$p.value)
  
  GNGs_ndi_smd_SDRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_SDRT_TD$p.value)
  
  GNGs_ndi_smd_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_SDRT_ADHD$p.value)
  
  GNGs_ndi_smd_CVRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_CVRT_TD$p.value)
  
  GNGs_ndi_smd_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_CVRT_ADHD$p.value)
  
  GNGs_ndi_smd_mu_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_mu_TD$p.value)
  
  GNGs_ndi_smd_mu_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_mu_ADHD$p.value)
  
  GNGs_ndi_smd_sigma_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_sigma_TD$p.value)
  
  GNGs_ndi_smd_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_sigma_ADHD$p.value)
  
  GNGs_ndi_smd_tau_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_tau_TD$p.value)
  
  GNGs_ndi_smd_tau_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_tau_ADHD$p.value)
  
  GNGs_ndi_smd_commission_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_commission_TD$p.value)
  
  GNGs_ndi_smd_commission_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_commission_ADHD$p.value)
  
  GNGs_ndi_smd_omission_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_omission_TD$p.value)
  
  GNGs_ndi_smd_omission_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGs_ndi_visual_meanRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_meanRT_TD$p.value)
  
  GNGs_ndi_visual_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_meanRT_ADHD$p.value)
  
  GNGs_ndi_visual_SDRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_SDRT_TD$p.value)
  
  GNGs_ndi_visual_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_SDRT_ADHD$p.value)
  
  GNGs_ndi_visual_CVRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_CVRT_TD$p.value)
  
  GNGs_ndi_visual_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_CVRT_ADHD$p.value)
  
  GNGs_ndi_visual_mu_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_mu_TD$p.value)
  
  GNGs_ndi_visual_mu_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_mu_ADHD$p.value)
  
  GNGs_ndi_visual_sigma_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_sigma_TD$p.value)
  
  GNGs_ndi_visual_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_sigma_ADHD$p.value)
  
  GNGs_ndi_visual_tau_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_tau_TD$p.value)
  
  GNGs_ndi_visual_tau_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_tau_ADHD$p.value)
  
  GNGs_ndi_visual_commission_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_commission_TD$p.value)
  
  GNGs_ndi_visual_commission_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_commission_ADHD$p.value)
  
  GNGs_ndi_visual_omission_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_ndi_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_omission_TD$p.value)
  
  GNGs_ndi_visual_omission_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGs_ndi_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_ndi_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_ndi_visual_omission_ADHD$p.value)
  
  
  ## Clustering Coefficient ##
  # Frontoparietal Network #
  GNGs_clustcoeff_fpn_meanRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_meanRT_TD$p.value)
  
  GNGs_clustcoeff_fpn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_meanRT_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_SDRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_SDRT_TD$p.value)
  
  GNGs_clustcoeff_fpn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_SDRT_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_CVRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_CVRT_TD$p.value)
  
  GNGs_clustcoeff_fpn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_CVRT_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_mu_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_mu_TD$p.value)
  
  GNGs_clustcoeff_fpn_mu_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_mu_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_sigma_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_sigma_TD$p.value)
  
  GNGs_clustcoeff_fpn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_sigma_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_tau_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_tau_TD$p.value)
  
  GNGs_clustcoeff_fpn_tau_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_tau_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_commission_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_commission_TD$p.value)
  
  GNGs_clustcoeff_fpn_commission_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_commission_ADHD$p.value)
  
  GNGs_clustcoeff_fpn_omission_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_omission_TD$p.value)
  
  GNGs_clustcoeff_fpn_omission_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGs_clustcoeff_con_meanRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_meanRT_TD$p.value)
  
  GNGs_clustcoeff_con_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_meanRT_ADHD$p.value)
  
  GNGs_clustcoeff_con_SDRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_SDRT_TD$p.value)
  
  GNGs_clustcoeff_con_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_SDRT_ADHD$p.value)
  
  GNGs_clustcoeff_con_CVRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_CVRT_TD$p.value)
  
  GNGs_clustcoeff_con_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_CVRT_ADHD$p.value)
  
  GNGs_clustcoeff_con_mu_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_mu_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_mu_TD$p.value)
  
  GNGs_clustcoeff_con_mu_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_mu_ADHD$p.value)
  
  GNGs_clustcoeff_con_sigma_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_sigma_TD$p.value)
  
  GNGs_clustcoeff_con_sigma_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_sigma_ADHD$p.value)
  
  GNGs_clustcoeff_con_tau_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_tau_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_tau_TD$p.value)
  
  GNGs_clustcoeff_con_tau_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_tau_ADHD$p.value)
  
  GNGs_clustcoeff_con_commission_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_commission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_commission_TD$p.value)
  
  GNGs_clustcoeff_con_commission_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_commission_ADHD$p.value)
  
  GNGs_clustcoeff_con_omission_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_con_omission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_omission_TD$p.value)
  
  GNGs_clustcoeff_con_omission_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGs_clustcoeff_dmn_meanRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_meanRT_TD$p.value)
  
  GNGs_clustcoeff_dmn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_meanRT_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_SDRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_SDRT_TD$p.value)
  
  GNGs_clustcoeff_dmn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_SDRT_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_CVRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_CVRT_TD$p.value)
  
  GNGs_clustcoeff_dmn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_CVRT_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_mu_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_mu_TD$p.value)
  
  GNGs_clustcoeff_dmn_mu_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_mu_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_sigma_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_sigma_TD$p.value)
  
  GNGs_clustcoeff_dmn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_sigma_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_tau_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_tau_TD$p.value)
  
  GNGs_clustcoeff_dmn_tau_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_tau_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_commission_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_commission_TD$p.value)
  
  GNGs_clustcoeff_dmn_commission_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_commission_ADHD$p.value)
  
  GNGs_clustcoeff_dmn_omission_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_omission_TD$p.value)
  
  GNGs_clustcoeff_dmn_omission_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGs_clustcoeff_smd_meanRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_meanRT_TD$p.value)
  
  GNGs_clustcoeff_smd_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_meanRT_ADHD$p.value)
  
  GNGs_clustcoeff_smd_SDRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_SDRT_TD$p.value)
  
  GNGs_clustcoeff_smd_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_SDRT_ADHD$p.value)
  
  GNGs_clustcoeff_smd_CVRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_CVRT_TD$p.value)
  
  GNGs_clustcoeff_smd_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_CVRT_ADHD$p.value)
  
  GNGs_clustcoeff_smd_mu_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_mu_TD$p.value)
  
  GNGs_clustcoeff_smd_mu_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_mu_ADHD$p.value)
  
  GNGs_clustcoeff_smd_sigma_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_sigma_TD$p.value)
  
  GNGs_clustcoeff_smd_sigma_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_sigma_ADHD$p.value)
  
  GNGs_clustcoeff_smd_tau_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_tau_TD$p.value)
  
  GNGs_clustcoeff_smd_tau_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_tau_ADHD$p.value)
  
  GNGs_clustcoeff_smd_commission_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_commission_TD$p.value)
  
  GNGs_clustcoeff_smd_commission_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_commission_ADHD$p.value)
  
  GNGs_clustcoeff_smd_omission_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_omission_TD$p.value)
  
  GNGs_clustcoeff_smd_omission_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGs_clustcoeff_visual_meanRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_meanRT_TD$p.value)
  
  GNGs_clustcoeff_visual_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_meanRT_ADHD$p.value)
  
  GNGs_clustcoeff_visual_SDRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_SDRT_TD$p.value)
  
  GNGs_clustcoeff_visual_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_SDRT_ADHD$p.value)
  
  GNGs_clustcoeff_visual_CVRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_CVRT_TD$p.value)
  
  GNGs_clustcoeff_visual_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_CVRT_ADHD$p.value)
  
  GNGs_clustcoeff_visual_mu_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_mu_TD$p.value)
  
  GNGs_clustcoeff_visual_mu_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_mu_ADHD$p.value)
  
  GNGs_clustcoeff_visual_sigma_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_sigma_TD$p.value)
  
  GNGs_clustcoeff_visual_sigma_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_sigma_ADHD$p.value)
  
  GNGs_clustcoeff_visual_tau_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_tau_TD$p.value)
  
  GNGs_clustcoeff_visual_tau_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_tau_ADHD$p.value)
  
  GNGs_clustcoeff_visual_commission_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_commission_TD$p.value)
  
  GNGs_clustcoeff_visual_commission_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_commission_ADHD$p.value)
  
  GNGs_clustcoeff_visual_omission_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_omission_TD$p.value)
  
  GNGs_clustcoeff_visual_omission_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Simple_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGs_clustcoeff_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGs_clustcoeff_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGs_clustcoeff_visual_omission_ADHD$p.value)
  
  
  #### Now that the comparison/threshold/diagnosis/r/p lists are made, add them to the temp_corr_data data frame and then bind this to the all_corr_data data frame ####
  temp_corr_data = data.frame(comparison = comparison_list, threshold = threshold_list, diagnosis = diagnosis_list, r.value = r_list, p.value = p_list)
  
  all_corr_data = rbind(all_corr_data, temp_corr_data)
  
}




#### Loop through the GNGr RData files generated above, compute partial brain x behavior correlations controlling for density, and note them down in the all_corr_data data frame ####

data_list_gngr <- c("./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh000.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh005.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh010.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh015.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh020.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh025.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh030.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh035.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh040.RData", "./brainxbehav_data_GNGr/brainbehav_data_GNGr_thresh045.RData")

# This loop looks at GNGr brain network organization and GNGr task performance #
for (i in 1:length(data_list_gngr)) {
  
  df_temp <- get(load(data_list_gngr[i]))
  temp_corr_data <- data.frame(comparison = NA,
                               threshold = NA,
                               diagnosis = NA,
                               r.value = NA,
                               p.value = NA)
  
  thresh = sapply(data_list_gngr[i], function (s) strsplit(strsplit(s, "_")[[1]][6], ".R")[[1]][1])  
  
  
  # Make 5 lists, 1 for each variable of the data frame
  comparison_list = c()
  threshold_list = c()
  diagnosis_list = c()
  r_list = c()
  p_list = c()
  
  
  #### Run each correlation of interest and append the relevant values to each list ####
  ### Whole-Brain Metrics ###
  ## Modularity ##
  GNGr_mod_meanRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_mod_meanRT_TD$p.value)
  
  GNGr_mod_meanRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_meanRT_ADHD$p.value)
  
  GNGr_mod_SDRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_mod_SDRT_TD$p.value)
  
  GNGr_mod_SDRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_SDRT_ADHD$p.value)
  
  GNGr_mod_CVRT_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_mod_CVRT_TD$p.value)
  
  GNGr_mod_CVRT_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_CVRT_ADHD$p.value)
  
  GNGr_mod_mu_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_mu_TD$estimate)
  p_list <- c(p_list, GNGr_mod_mu_TD$p.value)
  
  GNGr_mod_mu_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_mu_ADHD$p.value)
  
  GNGr_mod_sigma_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_mod_sigma_TD$p.value)
  
  GNGr_mod_sigma_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_sigma_ADHD$p.value)
  
  GNGr_mod_tau_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_tau_TD$estimate)
  p_list <- c(p_list, GNGr_mod_tau_TD$p.value)
  
  GNGr_mod_tau_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_tau_ADHD$p.value)
  
  GNGr_mod_commission_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_commission_TD$estimate)
  p_list <- c(p_list, GNGr_mod_commission_TD$p.value)
  
  GNGr_mod_commission_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_commission_ADHD$p.value)
  
  GNGr_mod_omission_TD <- pcor.test(df_temp$mod[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_mod_omission_TD$estimate)
  p_list <- c(p_list, GNGr_mod_omission_TD$p.value)
  
  GNGr_mod_omission_ADHD <- pcor.test(df_temp$mod[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_mod_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_mod_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_mod_omission_ADHD$p.value)
  
  
  ## Global Efficiency ##
  GNGr_ge_meanRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ge_meanRT_TD$p.value)
  
  GNGr_ge_meanRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_meanRT_ADHD$p.value)
  
  GNGr_ge_SDRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ge_SDRT_TD$p.value)
  
  GNGr_ge_SDRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_SDRT_ADHD$p.value)
  
  GNGr_ge_CVRT_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ge_CVRT_TD$p.value)
  
  GNGr_ge_CVRT_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_CVRT_ADHD$p.value)
  
  GNGr_ge_mu_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ge_mu_TD$p.value)
  
  GNGr_ge_mu_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_mu_ADHD$p.value)
  
  GNGr_ge_sigma_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ge_sigma_TD$p.value)
  
  GNGr_ge_sigma_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_sigma_ADHD$p.value)
  
  GNGr_ge_tau_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ge_tau_TD$p.value)
  
  GNGr_ge_tau_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_tau_ADHD$p.value)
  
  GNGr_ge_commission_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ge_commission_TD$p.value)
  
  GNGr_ge_commission_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_commission_ADHD$p.value)
  
  GNGr_ge_omission_TD <- pcor.test(df_temp$ge[df_temp$diag == 'TD'], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD'], df_temp$density[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ge_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ge_omission_TD$p.value)
  
  GNGr_ge_omission_ADHD <- pcor.test(df_temp$ge[df_temp$diag == 'ADHD'], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD'], df_temp$density[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'GNGr_ge_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ge_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ge_omission_ADHD$p.value)
  
  
  ### Nodal Metrics ###
  ## Participation Coefficient ##
  # Frontoparietal Network #
  GNGr_pc_fpn_meanRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_meanRT_TD$p.value)
  
  GNGr_pc_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_meanRT_ADHD$p.value)
  
  GNGr_pc_fpn_SDRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_SDRT_TD$p.value)
  
  GNGr_pc_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_SDRT_ADHD$p.value)
  
  GNGr_pc_fpn_CVRT_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_CVRT_TD$p.value)
  
  GNGr_pc_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_CVRT_ADHD$p.value)
  
  GNGr_pc_fpn_mu_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_mu_TD$p.value)
  
  GNGr_pc_fpn_mu_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_mu_ADHD$p.value)
  
  GNGr_pc_fpn_sigma_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_sigma_TD$p.value)
  
  GNGr_pc_fpn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_sigma_ADHD$p.value)
  
  GNGr_pc_fpn_tau_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_tau_TD$p.value)
  
  GNGr_pc_fpn_tau_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_tau_ADHD$p.value)
  
  GNGr_pc_fpn_commission_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_commission_TD$p.value)
  
  GNGr_pc_fpn_commission_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_commission_ADHD$p.value)
  
  GNGr_pc_fpn_omission_TD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_omission_TD$p.value)
  
  GNGr_pc_fpn_omission_ADHD <- pcor.test(df_temp$avg_pc_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3)])
  comparison_list <- c(comparison_list, 'GNGr_pc_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGr_pc_con_meanRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_meanRT_TD$p.value)
  
  GNGr_pc_con_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_meanRT_ADHD$p.value)
  
  GNGr_pc_con_SDRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_SDRT_TD$p.value)
  
  GNGr_pc_con_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_SDRT_ADHD$p.value)
  
  GNGr_pc_con_CVRT_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_CVRT_TD$p.value)
  
  GNGr_pc_con_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_CVRT_ADHD$p.value)
  
  GNGr_pc_con_mu_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_mu_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_mu_TD$p.value)
  
  GNGr_pc_con_mu_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_mu_ADHD$p.value)
  
  GNGr_pc_con_sigma_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_sigma_TD$p.value)
  
  GNGr_pc_con_sigma_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_sigma_ADHD$p.value)
  
  GNGr_pc_con_tau_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_tau_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_tau_TD$p.value)
  
  GNGr_pc_con_tau_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_tau_ADHD$p.value)
  
  GNGr_pc_con_commission_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_commission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_commission_TD$p.value)
  
  GNGr_pc_con_commission_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_commission_ADHD$p.value)
  
  GNGr_pc_con_omission_TD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_con_omission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_con_omission_TD$p.value)
  
  GNGr_pc_con_omission_ADHD <- pcor.test(df_temp$avg_pc_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1)])
  comparison_list <- c(comparison_list, 'GNGr_pc_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGr_pc_dmn_meanRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_meanRT_TD$p.value)
  
  GNGr_pc_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_meanRT_ADHD$p.value)
  
  GNGr_pc_dmn_SDRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_SDRT_TD$p.value)
  
  GNGr_pc_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_SDRT_ADHD$p.value)
  
  GNGr_pc_dmn_CVRT_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_CVRT_TD$p.value)
  
  GNGr_pc_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_CVRT_ADHD$p.value)
  
  GNGr_pc_dmn_mu_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_mu_TD$p.value)
  
  GNGr_pc_dmn_mu_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_mu_ADHD$p.value)
  
  GNGr_pc_dmn_sigma_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_sigma_TD$p.value)
  
  GNGr_pc_dmn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_sigma_ADHD$p.value)
  
  GNGr_pc_dmn_tau_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_tau_TD$p.value)
  
  GNGr_pc_dmn_tau_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_tau_ADHD$p.value)
  
  GNGr_pc_dmn_commission_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_commission_TD$p.value)
  
  GNGr_pc_dmn_commission_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_commission_ADHD$p.value)
  
  GNGr_pc_dmn_omission_TD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_omission_TD$p.value)
  
  GNGr_pc_dmn_omission_ADHD <- pcor.test(df_temp$avg_pc_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2)])
  comparison_list <- c(comparison_list, 'GNGr_pc_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGr_pc_smd_meanRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_meanRT_TD$p.value)
  
  GNGr_pc_smd_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_meanRT_ADHD$p.value)
  
  GNGr_pc_smd_SDRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_SDRT_TD$p.value)
  
  GNGr_pc_smd_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_SDRT_ADHD$p.value)
  
  GNGr_pc_smd_CVRT_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_CVRT_TD$p.value)
  
  GNGr_pc_smd_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_CVRT_ADHD$p.value)
  
  GNGr_pc_smd_mu_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_mu_TD$p.value)
  
  GNGr_pc_smd_mu_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_mu_ADHD$p.value)
  
  GNGr_pc_smd_sigma_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_sigma_TD$p.value)
  
  GNGr_pc_smd_sigma_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_sigma_ADHD$p.value)
  
  GNGr_pc_smd_tau_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_tau_TD$p.value)
  
  GNGr_pc_smd_tau_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_tau_ADHD$p.value)
  
  GNGr_pc_smd_commission_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_commission_TD$p.value)
  
  GNGr_pc_smd_commission_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_commission_ADHD$p.value)
  
  GNGr_pc_smd_omission_TD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_omission_TD$p.value)
  
  GNGr_pc_smd_omission_ADHD <- pcor.test(df_temp$avg_pc_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4)])
  comparison_list <- c(comparison_list, 'GNGr_pc_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGr_pc_visual_meanRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_meanRT_TD$p.value)
  
  GNGr_pc_visual_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_meanRT_ADHD$p.value)
  
  GNGr_pc_visual_SDRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_SDRT_TD$p.value)
  
  GNGr_pc_visual_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_SDRT_ADHD$p.value)
  
  GNGr_pc_visual_CVRT_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_CVRT_TD$p.value)
  
  GNGr_pc_visual_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_CVRT_ADHD$p.value)
  
  GNGr_pc_visual_mu_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_mu_TD$p.value)
  
  GNGr_pc_visual_mu_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_mu_ADHD$p.value)
  
  GNGr_pc_visual_sigma_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_sigma_TD$p.value)
  
  GNGr_pc_visual_sigma_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_sigma_ADHD$p.value)
  
  GNGr_pc_visual_tau_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_tau_TD$p.value)
  
  GNGr_pc_visual_tau_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_tau_ADHD$p.value)
  
  GNGr_pc_visual_commission_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_commission_TD$p.value)
  
  GNGr_pc_visual_commission_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_commission_ADHD$p.value)
  
  GNGr_pc_visual_omission_TD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_pc_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_omission_TD$p.value)
  
  GNGr_pc_visual_omission_ADHD <- pcor.test(df_temp$avg_pc_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5)])
  comparison_list <- c(comparison_list, 'GNGr_pc_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_pc_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_pc_visual_omission_ADHD$p.value)
  
  
  ## Node Dissociation Index ##
  # Frontoparietal Network #
  GNGr_ndi_fpn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_meanRT_TD$p.value)
  
  GNGr_ndi_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_meanRT_ADHD$p.value)
  
  GNGr_ndi_fpn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_SDRT_TD$p.value)
  
  GNGr_ndi_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_SDRT_ADHD$p.value)
  
  GNGr_ndi_fpn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_CVRT_TD$p.value)
  
  GNGr_ndi_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_CVRT_ADHD$p.value)
  
  GNGr_ndi_fpn_mu_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_mu_TD$p.value)
  
  GNGr_ndi_fpn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_mu_ADHD$p.value)
  
  GNGr_ndi_fpn_sigma_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_sigma_TD$p.value)
  
  GNGr_ndi_fpn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_sigma_ADHD$p.value)
  
  GNGr_ndi_fpn_tau_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_tau_TD$p.value)
  
  GNGr_ndi_fpn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_tau_ADHD$p.value)
  
  GNGr_ndi_fpn_commission_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_commission_TD$p.value)
  
  GNGr_ndi_fpn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_commission_ADHD$p.value)
  
  GNGr_ndi_fpn_omission_TD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_omission_TD$p.value)
  
  GNGr_ndi_fpn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net3[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGr_ndi_con_meanRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_meanRT_TD$p.value)
  
  GNGr_ndi_con_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_meanRT_ADHD$p.value)
  
  GNGr_ndi_con_SDRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_SDRT_TD$p.value)
  
  GNGr_ndi_con_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_SDRT_ADHD$p.value)
  
  GNGr_ndi_con_CVRT_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_CVRT_TD$p.value)
  
  GNGr_ndi_con_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_CVRT_ADHD$p.value)
  
  GNGr_ndi_con_mu_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_mu_TD$p.value)
  
  GNGr_ndi_con_mu_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_mu_ADHD$p.value)
  
  GNGr_ndi_con_sigma_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_sigma_TD$p.value)
  
  GNGr_ndi_con_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_sigma_ADHD$p.value)
  
  GNGr_ndi_con_tau_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_tau_TD$p.value)
  
  GNGr_ndi_con_tau_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_tau_ADHD$p.value)
  
  GNGr_ndi_con_commission_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_commission_TD$p.value)
  
  GNGr_ndi_con_commission_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_commission_ADHD$p.value)
  
  GNGr_ndi_con_omission_TD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_con_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_omission_TD$p.value)
  
  GNGr_ndi_con_omission_ADHD <- pcor.test(df_temp$avg_ndi_net1[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGr_ndi_dmn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_meanRT_TD$p.value)
  
  GNGr_ndi_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_meanRT_ADHD$p.value)
  
  GNGr_ndi_dmn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_SDRT_TD$p.value)
  
  GNGr_ndi_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_SDRT_ADHD$p.value)
  
  GNGr_ndi_dmn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_CVRT_TD$p.value)
  
  GNGr_ndi_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_CVRT_ADHD$p.value)
  
  GNGr_ndi_dmn_mu_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_mu_TD$p.value)
  
  GNGr_ndi_dmn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_mu_ADHD$p.value)
  
  GNGr_ndi_dmn_sigma_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_sigma_TD$p.value)
  
  GNGr_ndi_dmn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_sigma_ADHD$p.value)
  
  GNGr_ndi_dmn_tau_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_tau_TD$p.value)
  
  GNGr_ndi_dmn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_tau_ADHD$p.value)
  
  GNGr_ndi_dmn_commission_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_commission_TD$p.value)
  
  GNGr_ndi_dmn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_commission_ADHD$p.value)
  
  GNGr_ndi_dmn_omission_TD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_omission_TD$p.value)
  
  GNGr_ndi_dmn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net2[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGr_ndi_smd_meanRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_meanRT_TD$p.value)
  
  GNGr_ndi_smd_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_meanRT_ADHD$p.value)
  
  GNGr_ndi_smd_SDRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_SDRT_TD$p.value)
  
  GNGr_ndi_smd_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_SDRT_ADHD$p.value)
  
  GNGr_ndi_smd_CVRT_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_CVRT_TD$p.value)
  
  GNGr_ndi_smd_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_CVRT_ADHD$p.value)
  
  GNGr_ndi_smd_mu_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_mu_TD$p.value)
  
  GNGr_ndi_smd_mu_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_mu_ADHD$p.value)
  
  GNGr_ndi_smd_sigma_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_sigma_TD$p.value)
  
  GNGr_ndi_smd_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_sigma_ADHD$p.value)
  
  GNGr_ndi_smd_tau_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_tau_TD$p.value)
  
  GNGr_ndi_smd_tau_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_tau_ADHD$p.value)
  
  GNGr_ndi_smd_commission_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_commission_TD$p.value)
  
  GNGr_ndi_smd_commission_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_commission_ADHD$p.value)
  
  GNGr_ndi_smd_omission_TD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_omission_TD$p.value)
  
  GNGr_ndi_smd_omission_ADHD <- pcor.test(df_temp$avg_ndi_net4[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGr_ndi_visual_meanRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_meanRT_TD$p.value)
  
  GNGr_ndi_visual_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_meanRT_ADHD$p.value)
  
  GNGr_ndi_visual_SDRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_SDRT_TD$p.value)
  
  GNGr_ndi_visual_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_SDRT_ADHD$p.value)
  
  GNGr_ndi_visual_CVRT_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_CVRT_TD$p.value)
  
  GNGr_ndi_visual_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_CVRT_ADHD$p.value)
  
  GNGr_ndi_visual_mu_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_mu_TD$p.value)
  
  GNGr_ndi_visual_mu_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_mu_ADHD$p.value)
  
  GNGr_ndi_visual_sigma_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_sigma_TD$p.value)
  
  GNGr_ndi_visual_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_sigma_ADHD$p.value)
  
  GNGr_ndi_visual_tau_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_tau_TD$p.value)
  
  GNGr_ndi_visual_tau_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_tau_ADHD$p.value)
  
  GNGr_ndi_visual_commission_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_commission_TD$p.value)
  
  GNGr_ndi_visual_commission_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_commission_ADHD$p.value)
  
  GNGr_ndi_visual_omission_TD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_ndi_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_omission_TD$p.value)
  
  GNGr_ndi_visual_omission_ADHD <- pcor.test(df_temp$avg_ndi_net5[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5)])
  comparison_list <- c(comparison_list, 'GNGr_ndi_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_ndi_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_ndi_visual_omission_ADHD$p.value)
  
  
  ## Clustering Coefficient ##
  # Frontoparietal Network #
  GNGr_clustcoeff_fpn_meanRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_meanRT_TD$p.value)
  
  GNGr_clustcoeff_fpn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_meanRT_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_SDRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_SDRT_TD$p.value)
  
  GNGr_clustcoeff_fpn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_SDRT_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_CVRT_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_CVRT_TD$p.value)
  
  GNGr_clustcoeff_fpn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_CVRT_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_mu_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_mu_TD$p.value)
  
  GNGr_clustcoeff_fpn_mu_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_mu_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_sigma_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_sigma_TD$p.value)
  
  GNGr_clustcoeff_fpn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_sigma_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_tau_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_tau_TD$p.value)
  
  GNGr_clustcoeff_fpn_tau_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_tau_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_commission_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_commission_TD$p.value)
  
  GNGr_clustcoeff_fpn_commission_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_commission_ADHD$p.value)
  
  GNGr_clustcoeff_fpn_omission_TD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_omission_TD$p.value)
  
  GNGr_clustcoeff_fpn_omission_ADHD <- pcor.test(df_temp$clustcoeff_fpn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  GNGr_clustcoeff_con_meanRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_meanRT_TD$p.value)
  
  GNGr_clustcoeff_con_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_meanRT_ADHD$p.value)
  
  GNGr_clustcoeff_con_SDRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_SDRT_TD$p.value)
  
  GNGr_clustcoeff_con_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_SDRT_ADHD$p.value)
  
  GNGr_clustcoeff_con_CVRT_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_CVRT_TD$p.value)
  
  GNGr_clustcoeff_con_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_CVRT_ADHD$p.value)
  
  GNGr_clustcoeff_con_mu_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_mu_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_mu_TD$p.value)
  
  GNGr_clustcoeff_con_mu_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_mu_ADHD$p.value)
  
  GNGr_clustcoeff_con_sigma_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_sigma_TD$p.value)
  
  GNGr_clustcoeff_con_sigma_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_sigma_ADHD$p.value)
  
  GNGr_clustcoeff_con_tau_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_tau_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_tau_TD$p.value)
  
  GNGr_clustcoeff_con_tau_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_tau_ADHD$p.value)
  
  GNGr_clustcoeff_con_commission_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_commission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_commission_TD$p.value)
  
  GNGr_clustcoeff_con_commission_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_commission_ADHD$p.value)
  
  GNGr_clustcoeff_con_omission_TD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_con_omission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_omission_TD$p.value)
  
  GNGr_clustcoeff_con_omission_ADHD <- pcor.test(df_temp$clustcoeff_con[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_con_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  GNGr_clustcoeff_dmn_meanRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_meanRT_TD$p.value)
  
  GNGr_clustcoeff_dmn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_meanRT_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_SDRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_SDRT_TD$p.value)
  
  GNGr_clustcoeff_dmn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_SDRT_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_CVRT_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_CVRT_TD$p.value)
  
  GNGr_clustcoeff_dmn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_CVRT_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_mu_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_mu_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_mu_TD$p.value)
  
  GNGr_clustcoeff_dmn_mu_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_mu_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_sigma_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_sigma_TD$p.value)
  
  GNGr_clustcoeff_dmn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_sigma_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_tau_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_tau_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_tau_TD$p.value)
  
  GNGr_clustcoeff_dmn_tau_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_tau_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_commission_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_commission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_commission_TD$p.value)
  
  GNGr_clustcoeff_dmn_commission_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_commission_ADHD$p.value)
  
  GNGr_clustcoeff_dmn_omission_TD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_omission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_omission_TD$p.value)
  
  GNGr_clustcoeff_dmn_omission_ADHD <- pcor.test(df_temp$clustcoeff_dmn[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  GNGr_clustcoeff_smd_meanRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_meanRT_TD$p.value)
  
  GNGr_clustcoeff_smd_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_meanRT_ADHD$p.value)
  
  GNGr_clustcoeff_smd_SDRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_SDRT_TD$p.value)
  
  GNGr_clustcoeff_smd_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_SDRT_ADHD$p.value)
  
  GNGr_clustcoeff_smd_CVRT_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_CVRT_TD$p.value)
  
  GNGr_clustcoeff_smd_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_CVRT_ADHD$p.value)
  
  GNGr_clustcoeff_smd_mu_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_mu_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_mu_TD$p.value)
  
  GNGr_clustcoeff_smd_mu_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_mu_ADHD$p.value)
  
  GNGr_clustcoeff_smd_sigma_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_sigma_TD$p.value)
  
  GNGr_clustcoeff_smd_sigma_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_sigma_ADHD$p.value)
  
  GNGr_clustcoeff_smd_tau_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_tau_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_tau_TD$p.value)
  
  GNGr_clustcoeff_smd_tau_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_tau_ADHD$p.value)
  
  GNGr_clustcoeff_smd_commission_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_commission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_commission_TD$p.value)
  
  GNGr_clustcoeff_smd_commission_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_commission_ADHD$p.value)
  
  GNGr_clustcoeff_smd_omission_TD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_omission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_omission_TD$p.value)
  
  GNGr_clustcoeff_smd_omission_ADHD <- pcor.test(df_temp$clustcoeff_smd[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_smd_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  GNGr_clustcoeff_visual_meanRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_meanRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_meanRT_TD$p.value)
  
  GNGr_clustcoeff_visual_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Mean RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_meanRT_ADHD$p.value)
  
  GNGr_clustcoeff_visual_SDRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_SDRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_SDRT_TD$p.value)
  
  GNGr_clustcoeff_visual_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Std RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_SDRT_ADHD$p.value)
  
  GNGr_clustcoeff_visual_CVRT_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_CVRT_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_CVRT_TD$p.value)
  
  GNGr_clustcoeff_visual_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_CV RT'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_CVRT_ADHD$p.value)
  
  GNGr_clustcoeff_visual_mu_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_mu_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_mu_TD$p.value)
  
  GNGr_clustcoeff_visual_mu_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Mu'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_mu_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_mu_ADHD$p.value)
  
  GNGr_clustcoeff_visual_sigma_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_sigma_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_sigma_TD$p.value)
  
  GNGr_clustcoeff_visual_sigma_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Sigma'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_sigma_ADHD$p.value)
  
  GNGr_clustcoeff_visual_tau_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_tau_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_tau_TD$p.value)
  
  GNGr_clustcoeff_visual_tau_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Tau'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_tau_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_tau_ADHD$p.value)
  
  GNGr_clustcoeff_visual_commission_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_commission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_commission_TD$p.value)
  
  GNGr_clustcoeff_visual_commission_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Commision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_commission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_commission_ADHD$p.value)
  
  GNGr_clustcoeff_visual_omission_TD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_omission_TD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_omission_TD$p.value)
  
  GNGr_clustcoeff_visual_omission_ADHD <- pcor.test(df_temp$clustcoeff_visual[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$'Overall_Repeat_Ommision Rate'[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)], df_temp$density[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual)])
  comparison_list <- c(comparison_list, 'GNGr_clustcoeff_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, GNGr_clustcoeff_visual_omission_ADHD$estimate)
  p_list <- c(p_list, GNGr_clustcoeff_visual_omission_ADHD$p.value)
  
  
  #### Now that the comparison/threshold/diagnosis/r/p lists are made, add them to the temp_corr_data data frame and then bind this to the all_corr_data data frame ####
  temp_corr_data = data.frame(comparison = comparison_list, threshold = threshold_list, diagnosis = diagnosis_list, r.value = r_list, p.value = p_list)
  
  all_corr_data = rbind(all_corr_data, temp_corr_data)
  
}



#### Loop through the GNGs RData files generated above, compute partial brain x behavior correlations controlling for density, and note them down in the all_corr_data data frame ####

data_list_diff <- c("./brainxbehav_data_diff/brainbehav_data_diff_thresh000.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh005.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh010.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh015.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh020.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh025.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh030.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh035.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh040.RData", "./brainxbehav_data_diff/brainbehav_data_diff_thresh045.RData")


# This loop looks at how changes in brain network organization from GNGs to GNGr are associated with changes in behavioral performance from GNGs to GNGr #
for (i in 1:length(data_list_diff)) {
  
  df_temp <- get(load(data_list_diff[i]))
  temp_corr_data <- data.frame(comparison = NA,
                               threshold = NA,
                               diagnosis = NA,
                               r.value = NA,
                               p.value = NA)
  
  thresh = sapply(data_list_diff[i], function (s) strsplit(strsplit(s, "_")[[1]][6], ".R")[[1]][1])  
  
  
  # Make 5 lists, 1 for each variable of the data frame
  comparison_list = c()
  threshold_list = c()
  diagnosis_list = c()
  r_list = c()
  p_list = c()
  
  
  #### Run each correlation of interest and append the relevant values to each list ####
  ### Whole-Brain Metrics ###
  ## Modularity ##
  diff_mod_meanRT_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_MeanRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_meanRT_TD$estimate)
  p_list <- c(p_list, diff_mod_meanRT_TD$p.value)
  
  diff_mod_meanRT_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_MeanRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_mod_meanRT_ADHD$p.value)
  
  diff_mod_SDRT_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_SDRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_SDRT_TD$estimate)
  p_list <- c(p_list, diff_mod_SDRT_TD$p.value)
  
  diff_mod_SDRT_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_SDRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_mod_SDRT_ADHD$p.value)
  
  diff_mod_CVRT_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_CVRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_CVRT_TD$estimate)
  p_list <- c(p_list, diff_mod_CVRT_TD$p.value)
  
  diff_mod_CVRT_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_CVRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_mod_CVRT_ADHD$p.value)
  
  diff_mod_mu_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_Mu_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_mu_TD$estimate)
  p_list <- c(p_list, diff_mod_mu_TD$p.value)
  
  diff_mod_mu_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Mu_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_mu_ADHD$estimate)
  p_list <- c(p_list, diff_mod_mu_ADHD$p.value)
  
  diff_mod_sigma_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_Sigma_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_sigma_TD$estimate)
  p_list <- c(p_list, diff_mod_sigma_TD$p.value)
  
  diff_mod_sigma_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Sigma_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_mod_sigma_ADHD$p.value)
  
  diff_mod_tau_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_Tau_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_tau_TD$estimate)
  p_list <- c(p_list, diff_mod_tau_TD$p.value)
  
  diff_mod_tau_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Tau_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_tau_ADHD$estimate)
  p_list <- c(p_list, diff_mod_tau_ADHD$p.value)
  
  diff_mod_commission_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_Commission_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_commission_TD$estimate)
  p_list <- c(p_list, diff_mod_commission_TD$p.value)
  
  diff_mod_commission_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Commission_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_commission_ADHD$estimate)
  p_list <- c(p_list, diff_mod_commission_ADHD$p.value)
  
  diff_mod_omission_TD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'TD'], df_temp$'Overall_Omission_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_mod_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_mod_omission_TD$estimate)
  p_list <- c(p_list, diff_mod_omission_TD$p.value)
  
  diff_mod_omission_ADHD <- pcor.test(df_temp$mod_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Omission_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_mod_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_mod_omission_ADHD$estimate)
  p_list <- c(p_list, diff_mod_omission_ADHD$p.value)
  
  
  ## Global Efficiency ##
  diff_ge_meanRT_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_MeanRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ge_meanRT_TD$p.value)
  
  diff_ge_meanRT_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_MeanRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ge_meanRT_ADHD$p.value)
  
  diff_ge_SDRT_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_SDRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ge_SDRT_TD$p.value)
  
  diff_ge_SDRT_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_SDRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ge_SDRT_ADHD$p.value)
  
  diff_ge_CVRT_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_CVRT_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ge_CVRT_TD$p.value)
  
  diff_ge_CVRT_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_CVRT_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ge_CVRT_ADHD$p.value)
  
  diff_ge_mu_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_Mu_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_mu_TD$estimate)
  p_list <- c(p_list, diff_ge_mu_TD$p.value)
  
  diff_ge_mu_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Mu_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ge_mu_ADHD$p.value)
  
  diff_ge_sigma_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_Sigma_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_sigma_TD$estimate)
  p_list <- c(p_list, diff_ge_sigma_TD$p.value)
  
  diff_ge_sigma_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Sigma_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ge_sigma_ADHD$p.value)
  
  diff_ge_tau_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_Tau_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_tau_TD$estimate)
  p_list <- c(p_list, diff_ge_tau_TD$p.value)
  
  diff_ge_tau_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Tau_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ge_tau_ADHD$p.value)
  
  diff_ge_commission_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_Commission_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_commission_TD$estimate)
  p_list <- c(p_list, diff_ge_commission_TD$p.value)
  
  diff_ge_commission_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Commission_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ge_commission_ADHD$p.value)
  
  diff_ge_omission_TD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'TD'], df_temp$'Overall_Omission_diff'[df_temp$diag == 'TD'], df_temp$density_diff[df_temp$diag == 'TD'])
  comparison_list <- c(comparison_list, 'diff_ge_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ge_omission_TD$estimate)
  p_list <- c(p_list, diff_ge_omission_TD$p.value)
  
  diff_ge_omission_ADHD <- pcor.test(df_temp$ge_diff[df_temp$diag == 'ADHD'], df_temp$'Overall_Omission_diff'[df_temp$diag == 'ADHD'], df_temp$density_diff[df_temp$diag == 'ADHD'])
  comparison_list <- c(comparison_list, 'diff_ge_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ge_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ge_omission_ADHD$p.value)
  
  
  ### Nodal Metrics ###
  ## Participation Coefficient ##
  # Frontoparietal Network #
  diff_pc_fpn_meanRT_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_meanRT_TD$p.value)
  
  diff_pc_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_meanRT_ADHD$p.value)
  
  diff_pc_fpn_SDRT_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_SDRT_TD$p.value)
  
  diff_pc_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_SDRT_ADHD$p.value)
  
  diff_pc_fpn_CVRT_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_CVRT_TD$p.value)
  
  diff_pc_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_CVRT_ADHD$p.value)
  
  diff_pc_fpn_mu_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_mu_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_mu_TD$p.value)
  
  diff_pc_fpn_mu_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_mu_ADHD$p.value)
  
  diff_pc_fpn_sigma_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_sigma_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_sigma_TD$p.value)
  
  diff_pc_fpn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_sigma_ADHD$p.value)
  
  diff_pc_fpn_tau_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_tau_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_tau_TD$p.value)
  
  diff_pc_fpn_tau_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_tau_ADHD$p.value)
  
  diff_pc_fpn_commission_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_commission_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_commission_TD$p.value)
  
  diff_pc_fpn_commission_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_commission_ADHD$p.value)
  
  diff_pc_fpn_omission_TD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_fpn_omission_TD$estimate)
  p_list <- c(p_list, diff_pc_fpn_omission_TD$p.value)
  
  diff_pc_fpn_omission_ADHD <- pcor.test(df_temp$avg_pc_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  diff_pc_con_meanRT_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_meanRT_TD$estimate)
  p_list <- c(p_list, diff_pc_con_meanRT_TD$p.value)
  
  diff_pc_con_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_meanRT_ADHD$p.value)
  
  diff_pc_con_SDRT_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_SDRT_TD$estimate)
  p_list <- c(p_list, diff_pc_con_SDRT_TD$p.value)
  
  diff_pc_con_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_SDRT_ADHD$p.value)
  
  diff_pc_con_CVRT_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_CVRT_TD$estimate)
  p_list <- c(p_list, diff_pc_con_CVRT_TD$p.value)
  
  diff_pc_con_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_CVRT_ADHD$p.value)
  
  diff_pc_con_mu_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_mu_TD$estimate)
  p_list <- c(p_list, diff_pc_con_mu_TD$p.value)
  
  diff_pc_con_mu_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_mu_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_mu_ADHD$p.value)
  
  diff_pc_con_sigma_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_sigma_TD$estimate)
  p_list <- c(p_list, diff_pc_con_sigma_TD$p.value)
  
  diff_pc_con_sigma_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_sigma_ADHD$p.value)
  
  diff_pc_con_tau_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_tau_TD$estimate)
  p_list <- c(p_list, diff_pc_con_tau_TD$p.value)
  
  diff_pc_con_tau_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_tau_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_tau_ADHD$p.value)
  
  diff_pc_con_commission_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_commission_TD$estimate)
  p_list <- c(p_list, diff_pc_con_commission_TD$p.value)
  
  diff_pc_con_commission_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_commission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_commission_ADHD$p.value)
  
  diff_pc_con_omission_TD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_con_omission_TD$estimate)
  p_list <- c(p_list, diff_pc_con_omission_TD$p.value)
  
  diff_pc_con_omission_ADHD <- pcor.test(df_temp$avg_pc_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_con_omission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  diff_pc_dmn_meanRT_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_meanRT_TD$p.value)
  
  diff_pc_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_meanRT_ADHD$p.value)
  
  diff_pc_dmn_SDRT_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_SDRT_TD$p.value)
  
  diff_pc_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_SDRT_ADHD$p.value)
  
  diff_pc_dmn_CVRT_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_CVRT_TD$p.value)
  
  diff_pc_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_CVRT_ADHD$p.value)
  
  diff_pc_dmn_mu_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_mu_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_mu_TD$p.value)
  
  diff_pc_dmn_mu_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_mu_ADHD$p.value)
  
  diff_pc_dmn_sigma_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_sigma_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_sigma_TD$p.value)
  
  diff_pc_dmn_sigma_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_sigma_ADHD$p.value)
  
  diff_pc_dmn_tau_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_tau_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_tau_TD$p.value)
  
  diff_pc_dmn_tau_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_tau_ADHD$p.value)
  
  diff_pc_dmn_commission_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_commission_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_commission_TD$p.value)
  
  diff_pc_dmn_commission_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_commission_ADHD$p.value)
  
  diff_pc_dmn_omission_TD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_dmn_omission_TD$estimate)
  p_list <- c(p_list, diff_pc_dmn_omission_TD$p.value)
  
  diff_pc_dmn_omission_ADHD <- pcor.test(df_temp$avg_pc_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  diff_pc_smd_meanRT_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_meanRT_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_meanRT_TD$p.value)
  
  diff_pc_smd_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_meanRT_ADHD$p.value)
  
  diff_pc_smd_SDRT_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_SDRT_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_SDRT_TD$p.value)
  
  diff_pc_smd_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_SDRT_ADHD$p.value)
  
  diff_pc_smd_CVRT_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_CVRT_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_CVRT_TD$p.value)
  
  diff_pc_smd_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_CVRT_ADHD$p.value)
  
  diff_pc_smd_mu_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_mu_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_mu_TD$p.value)
  
  diff_pc_smd_mu_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_mu_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_mu_ADHD$p.value)
  
  diff_pc_smd_sigma_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_sigma_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_sigma_TD$p.value)
  
  diff_pc_smd_sigma_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_sigma_ADHD$p.value)
  
  diff_pc_smd_tau_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_tau_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_tau_TD$p.value)
  
  diff_pc_smd_tau_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_tau_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_tau_ADHD$p.value)
  
  diff_pc_smd_commission_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_commission_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_commission_TD$p.value)
  
  diff_pc_smd_commission_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_commission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_commission_ADHD$p.value)
  
  diff_pc_smd_omission_TD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_smd_omission_TD$estimate)
  p_list <- c(p_list, diff_pc_smd_omission_TD$p.value)
  
  diff_pc_smd_omission_ADHD <- pcor.test(df_temp$avg_pc_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_smd_omission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  diff_pc_visual_meanRT_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_meanRT_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_meanRT_TD$p.value)
  
  diff_pc_visual_meanRT_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_meanRT_ADHD$p.value)
  
  diff_pc_visual_SDRT_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_SDRT_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_SDRT_TD$p.value)
  
  diff_pc_visual_SDRT_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_SDRT_ADHD$p.value)
  
  diff_pc_visual_CVRT_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_CVRT_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_CVRT_TD$p.value)
  
  diff_pc_visual_CVRT_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_CVRT_ADHD$p.value)
  
  diff_pc_visual_mu_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_mu_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_mu_TD$p.value)
  
  diff_pc_visual_mu_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_mu_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_mu_ADHD$p.value)
  
  diff_pc_visual_sigma_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_sigma_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_sigma_TD$p.value)
  
  diff_pc_visual_sigma_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_sigma_ADHD$p.value)
  
  diff_pc_visual_tau_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_tau_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_tau_TD$p.value)
  
  diff_pc_visual_tau_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_tau_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_tau_ADHD$p.value)
  
  diff_pc_visual_commission_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_commission_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_commission_TD$p.value)
  
  diff_pc_visual_commission_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_commission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_commission_ADHD$p.value)
  
  diff_pc_visual_omission_TD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_pc_visual_omission_TD$estimate)
  p_list <- c(p_list, diff_pc_visual_omission_TD$p.value)
  
  diff_pc_visual_omission_ADHD <- pcor.test(df_temp$avg_pc_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_pc_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_pc_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_pc_visual_omission_ADHD$estimate)
  p_list <- c(p_list, diff_pc_visual_omission_ADHD$p.value)
  
  
  ## Node Dissociation Index ##
  # Frontoparietal Network #
  diff_ndi_fpn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_meanRT_TD$p.value)
  
  diff_ndi_fpn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_meanRT_ADHD$p.value)
  
  diff_ndi_fpn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_SDRT_TD$p.value)
  
  diff_ndi_fpn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_SDRT_ADHD$p.value)
  
  diff_ndi_fpn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_CVRT_TD$p.value)
  
  diff_ndi_fpn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_CVRT_ADHD$p.value)
  
  diff_ndi_fpn_mu_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_mu_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_mu_TD$p.value)
  
  diff_ndi_fpn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_mu_ADHD$p.value)
  
  diff_ndi_fpn_sigma_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_sigma_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_sigma_TD$p.value)
  
  diff_ndi_fpn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_sigma_ADHD$p.value)
  
  diff_ndi_fpn_tau_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_tau_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_tau_TD$p.value)
  
  diff_ndi_fpn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_tau_ADHD$p.value)
  
  diff_ndi_fpn_commission_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_commission_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_commission_TD$p.value)
  
  diff_ndi_fpn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_commission_ADHD$p.value)
  
  diff_ndi_fpn_omission_TD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_fpn_omission_TD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_omission_TD$p.value)
  
  diff_ndi_fpn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net3_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net3_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  diff_ndi_con_meanRT_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_meanRT_TD$p.value)
  
  diff_ndi_con_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_meanRT_ADHD$p.value)
  
  diff_ndi_con_SDRT_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_SDRT_TD$p.value)
  
  diff_ndi_con_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_SDRT_ADHD$p.value)
  
  diff_ndi_con_CVRT_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_CVRT_TD$p.value)
  
  diff_ndi_con_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_CVRT_ADHD$p.value)
  
  diff_ndi_con_mu_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_mu_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_mu_TD$p.value)
  
  diff_ndi_con_mu_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_mu_ADHD$p.value)
  
  diff_ndi_con_sigma_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_sigma_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_sigma_TD$p.value)
  
  diff_ndi_con_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_sigma_ADHD$p.value)
  
  diff_ndi_con_tau_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_tau_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_tau_TD$p.value)
  
  diff_ndi_con_tau_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_tau_ADHD$p.value)
  
  diff_ndi_con_commission_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_commission_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_commission_TD$p.value)
  
  diff_ndi_con_commission_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_commission_ADHD$p.value)
  
  diff_ndi_con_omission_TD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_con_omission_TD$estimate)
  p_list <- c(p_list, diff_ndi_con_omission_TD$p.value)
  
  diff_ndi_con_omission_ADHD <- pcor.test(df_temp$avg_ndi_net1_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net1_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_con_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  diff_ndi_dmn_meanRT_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_meanRT_TD$p.value)
  
  diff_ndi_dmn_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_meanRT_ADHD$p.value)
  
  diff_ndi_dmn_SDRT_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_SDRT_TD$p.value)
  
  diff_ndi_dmn_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_SDRT_ADHD$p.value)
  
  diff_ndi_dmn_CVRT_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_CVRT_TD$p.value)
  
  diff_ndi_dmn_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_CVRT_ADHD$p.value)
  
  diff_ndi_dmn_mu_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_mu_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_mu_TD$p.value)
  
  diff_ndi_dmn_mu_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_mu_ADHD$p.value)
  
  diff_ndi_dmn_sigma_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_sigma_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_sigma_TD$p.value)
  
  diff_ndi_dmn_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_sigma_ADHD$p.value)
  
  diff_ndi_dmn_tau_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_tau_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_tau_TD$p.value)
  
  diff_ndi_dmn_tau_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_tau_ADHD$p.value)
  
  diff_ndi_dmn_commission_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_commission_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_commission_TD$p.value)
  
  diff_ndi_dmn_commission_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_commission_ADHD$p.value)
  
  diff_ndi_dmn_omission_TD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_dmn_omission_TD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_omission_TD$p.value)
  
  diff_ndi_dmn_omission_ADHD <- pcor.test(df_temp$avg_ndi_net2_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net2_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  diff_ndi_smd_meanRT_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_meanRT_TD$p.value)
  
  diff_ndi_smd_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_meanRT_ADHD$p.value)
  
  diff_ndi_smd_SDRT_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_SDRT_TD$p.value)
  
  diff_ndi_smd_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_SDRT_ADHD$p.value)
  
  diff_ndi_smd_CVRT_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_CVRT_TD$p.value)
  
  diff_ndi_smd_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_CVRT_ADHD$p.value)
  
  diff_ndi_smd_mu_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_mu_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_mu_TD$p.value)
  
  diff_ndi_smd_mu_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_mu_ADHD$p.value)
  
  diff_ndi_smd_sigma_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_sigma_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_sigma_TD$p.value)
  
  diff_ndi_smd_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_sigma_ADHD$p.value)
  
  diff_ndi_smd_tau_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_tau_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_tau_TD$p.value)
  
  diff_ndi_smd_tau_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_tau_ADHD$p.value)
  
  diff_ndi_smd_commission_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_commission_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_commission_TD$p.value)
  
  diff_ndi_smd_commission_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_commission_ADHD$p.value)
  
  diff_ndi_smd_omission_TD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_smd_omission_TD$estimate)
  p_list <- c(p_list, diff_ndi_smd_omission_TD$p.value)
  
  diff_ndi_smd_omission_ADHD <- pcor.test(df_temp$avg_ndi_net4_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net4_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_smd_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  diff_ndi_visual_meanRT_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_meanRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_meanRT_TD$p.value)
  
  diff_ndi_visual_meanRT_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_meanRT_ADHD$p.value)
  
  diff_ndi_visual_SDRT_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_SDRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_SDRT_TD$p.value)
  
  diff_ndi_visual_SDRT_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_SDRT_ADHD$p.value)
  
  diff_ndi_visual_CVRT_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_CVRT_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_CVRT_TD$p.value)
  
  diff_ndi_visual_CVRT_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_CVRT_ADHD$p.value)
  
  diff_ndi_visual_mu_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_mu_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_mu_TD$p.value)
  
  diff_ndi_visual_mu_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_mu_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_mu_ADHD$p.value)
  
  diff_ndi_visual_sigma_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_sigma_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_sigma_TD$p.value)
  
  diff_ndi_visual_sigma_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_sigma_ADHD$p.value)
  
  diff_ndi_visual_tau_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_tau_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_tau_TD$p.value)
  
  diff_ndi_visual_tau_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_tau_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_tau_ADHD$p.value)
  
  diff_ndi_visual_commission_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_commission_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_commission_TD$p.value)
  
  diff_ndi_visual_commission_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_commission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_commission_ADHD$p.value)
  
  diff_ndi_visual_omission_TD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_ndi_visual_omission_TD$estimate)
  p_list <- c(p_list, diff_ndi_visual_omission_TD$p.value)
  
  diff_ndi_visual_omission_ADHD <- pcor.test(df_temp$avg_ndi_net5_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$avg_ndi_net5_diff)])
  comparison_list <- c(comparison_list, 'diff_ndi_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_ndi_visual_omission_ADHD$estimate)
  p_list <- c(p_list, diff_ndi_visual_omission_ADHD$p.value)
  
  
  ## Clustering Coefficient ##
  # Frontoparietal Network #
  diff_clustcoeff_fpn_meanRT_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_meanRT_TD$p.value)
  
  diff_clustcoeff_fpn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_meanRT_ADHD$p.value)
  
  diff_clustcoeff_fpn_SDRT_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_SDRT_TD$p.value)
  
  diff_clustcoeff_fpn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_SDRT_ADHD$p.value)
  
  diff_clustcoeff_fpn_CVRT_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_CVRT_TD$p.value)
  
  diff_clustcoeff_fpn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_CVRT_ADHD$p.value)
  
  diff_clustcoeff_fpn_mu_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_mu_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_mu_TD$p.value)
  
  diff_clustcoeff_fpn_mu_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_mu_ADHD$p.value)
  
  diff_clustcoeff_fpn_sigma_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_sigma_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_sigma_TD$p.value)
  
  diff_clustcoeff_fpn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_sigma_ADHD$p.value)
  
  diff_clustcoeff_fpn_tau_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_tau_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_tau_TD$p.value)
  
  diff_clustcoeff_fpn_tau_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_tau_ADHD$p.value)
  
  diff_clustcoeff_fpn_commission_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_commission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_commission_TD$p.value)
  
  diff_clustcoeff_fpn_commission_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_commission_ADHD$p.value)
  
  diff_clustcoeff_fpn_omission_TD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_fpn_omission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_omission_TD$p.value)
  
  diff_clustcoeff_fpn_omission_ADHD <- pcor.test(df_temp$clustcoeff_fpn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_fpn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_fpn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_fpn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_fpn_omission_ADHD$p.value)
  
  
  # Cingulo-Opercular Network #
  diff_clustcoeff_con_meanRT_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_meanRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_meanRT_TD$p.value)
  
  diff_clustcoeff_con_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_meanRT_ADHD$p.value)
  
  diff_clustcoeff_con_SDRT_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_SDRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_SDRT_TD$p.value)
  
  diff_clustcoeff_con_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_SDRT_ADHD$p.value)
  
  diff_clustcoeff_con_CVRT_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_CVRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_CVRT_TD$p.value)
  
  diff_clustcoeff_con_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_CVRT_ADHD$p.value)
  
  diff_clustcoeff_con_mu_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_mu_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_mu_TD$p.value)
  
  diff_clustcoeff_con_mu_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_mu_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_mu_ADHD$p.value)
  
  diff_clustcoeff_con_sigma_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_sigma_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_sigma_TD$p.value)
  
  diff_clustcoeff_con_sigma_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_sigma_ADHD$p.value)
  
  diff_clustcoeff_con_tau_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_tau_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_tau_TD$p.value)
  
  diff_clustcoeff_con_tau_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_tau_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_tau_ADHD$p.value)
  
  diff_clustcoeff_con_commission_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_commission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_commission_TD$p.value)
  
  diff_clustcoeff_con_commission_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_commission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_commission_ADHD$p.value)
  
  diff_clustcoeff_con_omission_TD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_con_omission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_omission_TD$p.value)
  
  diff_clustcoeff_con_omission_ADHD <- pcor.test(df_temp$clustcoeff_con_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_con_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_con_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_con_omission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_con_omission_ADHD$p.value)
  
  
  # Default Mode Network #
  diff_clustcoeff_dmn_meanRT_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_meanRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_meanRT_TD$p.value)
  
  diff_clustcoeff_dmn_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_meanRT_ADHD$p.value)
  
  diff_clustcoeff_dmn_SDRT_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_SDRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_SDRT_TD$p.value)
  
  diff_clustcoeff_dmn_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_SDRT_ADHD$p.value)
  
  diff_clustcoeff_dmn_CVRT_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_CVRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_CVRT_TD$p.value)
  
  diff_clustcoeff_dmn_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_CVRT_ADHD$p.value)
  
  diff_clustcoeff_dmn_mu_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_mu_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_mu_TD$p.value)
  
  diff_clustcoeff_dmn_mu_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_mu_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_mu_ADHD$p.value)
  
  diff_clustcoeff_dmn_sigma_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_sigma_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_sigma_TD$p.value)
  
  diff_clustcoeff_dmn_sigma_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_sigma_ADHD$p.value)
  
  diff_clustcoeff_dmn_tau_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_tau_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_tau_TD$p.value)
  
  diff_clustcoeff_dmn_tau_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_tau_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_tau_ADHD$p.value)
  
  diff_clustcoeff_dmn_commission_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_commission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_commission_TD$p.value)
  
  diff_clustcoeff_dmn_commission_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_commission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_commission_ADHD$p.value)
  
  diff_clustcoeff_dmn_omission_TD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_dmn_omission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_omission_TD$p.value)
  
  diff_clustcoeff_dmn_omission_ADHD <- pcor.test(df_temp$clustcoeff_dmn_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_dmn_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_dmn_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_dmn_omission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_dmn_omission_ADHD$p.value)
  
  
  # Somatomotor Dorsal Network #
  diff_clustcoeff_smd_meanRT_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_meanRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_meanRT_TD$p.value)
  
  diff_clustcoeff_smd_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_meanRT_ADHD$p.value)
  
  diff_clustcoeff_smd_SDRT_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_SDRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_SDRT_TD$p.value)
  
  diff_clustcoeff_smd_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_SDRT_ADHD$p.value)
  
  diff_clustcoeff_smd_CVRT_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_CVRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_CVRT_TD$p.value)
  
  diff_clustcoeff_smd_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_CVRT_ADHD$p.value)
  
  diff_clustcoeff_smd_mu_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_mu_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_mu_TD$p.value)
  
  diff_clustcoeff_smd_mu_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_mu_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_mu_ADHD$p.value)
  
  diff_clustcoeff_smd_sigma_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_sigma_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_sigma_TD$p.value)
  
  diff_clustcoeff_smd_sigma_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_sigma_ADHD$p.value)
  
  diff_clustcoeff_smd_tau_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_tau_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_tau_TD$p.value)
  
  diff_clustcoeff_smd_tau_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_tau_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_tau_ADHD$p.value)
  
  diff_clustcoeff_smd_commission_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_commission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_commission_TD$p.value)
  
  diff_clustcoeff_smd_commission_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_commission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_commission_ADHD$p.value)
  
  diff_clustcoeff_smd_omission_TD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_smd_omission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_omission_TD$p.value)
  
  diff_clustcoeff_smd_omission_ADHD <- pcor.test(df_temp$clustcoeff_smd_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_smd_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_smd_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_smd_omission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_smd_omission_ADHD$p.value)
  
  
  # Visual Network #
  diff_clustcoeff_visual_meanRT_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_meanRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_meanRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_meanRT_TD$p.value)
  
  diff_clustcoeff_visual_meanRT_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_MeanRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_meanRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_meanRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_meanRT_ADHD$p.value)
  
  diff_clustcoeff_visual_SDRT_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_SDRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_SDRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_SDRT_TD$p.value)
  
  diff_clustcoeff_visual_SDRT_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_SDRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_SDRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_SDRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_SDRT_ADHD$p.value)
  
  diff_clustcoeff_visual_CVRT_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_CVRT_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_CVRT_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_CVRT_TD$p.value)
  
  diff_clustcoeff_visual_CVRT_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_CVRT_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_CVRT_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_CVRT_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_CVRT_ADHD$p.value)
  
  diff_clustcoeff_visual_mu_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_mu_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_mu_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_mu_TD$p.value)
  
  diff_clustcoeff_visual_mu_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Mu_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_mu_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_mu_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_mu_ADHD$p.value)
  
  diff_clustcoeff_visual_sigma_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_sigma_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_sigma_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_sigma_TD$p.value)
  
  diff_clustcoeff_visual_sigma_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Sigma_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_sigma_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_sigma_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_sigma_ADHD$p.value)
  
  diff_clustcoeff_visual_tau_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_tau_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_tau_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_tau_TD$p.value)
  
  diff_clustcoeff_visual_tau_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Tau_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_tau_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_tau_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_tau_ADHD$p.value)
  
  diff_clustcoeff_visual_commission_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_commission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_commission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_commission_TD$p.value)
  
  diff_clustcoeff_visual_commission_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Commission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_commission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_commission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_commission_ADHD$p.value)
  
  diff_clustcoeff_visual_omission_TD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'TD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_omission_TD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'TD')
  r_list <- c(r_list, diff_clustcoeff_visual_omission_TD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_omission_TD$p.value)
  
  diff_clustcoeff_visual_omission_ADHD <- pcor.test(df_temp$clustcoeff_visual_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$Overall_Omission_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)], df_temp$density_diff[df_temp$diag == 'ADHD' & !is.na(df_temp$clustcoeff_visual_diff)])
  comparison_list <- c(comparison_list, 'diff_clustcoeff_visual_omission_ADHD')
  threshold_list <- c(threshold_list, thresh)
  diagnosis_list <- c(diagnosis_list, 'ADHD')
  r_list <- c(r_list, diff_clustcoeff_visual_omission_ADHD$estimate)
  p_list <- c(p_list, diff_clustcoeff_visual_omission_ADHD$p.value)
 
  
  #### Now that the comparison/threshold/diagnosis/r/p lists are made, add them to the temp_corr_data data frame and then bind this to the all_corr_data data frame ####
  temp_corr_data = data.frame(comparison = comparison_list, threshold = threshold_list, diagnosis = diagnosis_list, r.value = r_list, p.value = p_list)
  
  all_corr_data = rbind(all_corr_data, temp_corr_data)
  
}


write.csv(all_corr_data, file = 'all_corr_data.csv', row.names = FALSE)



##### Brain x Behavior Correlation Plots Across Thresholds #####
all_corr_data <- read.csv('all_corr_data.csv', sep = ',')
all_corr_data$threshold <- gsub("[a-z]", "", all_corr_data$threshold)
all_corr_data$threshold[all_corr_data$threshold == '000'] <- '0.00'
all_corr_data$threshold[all_corr_data$threshold == '005'] <- '0.05'
all_corr_data$threshold[all_corr_data$threshold == '010'] <- '0.10'
all_corr_data$threshold[all_corr_data$threshold == '015'] <- '0.15'
all_corr_data$threshold[all_corr_data$threshold == '020'] <- '0.20'
all_corr_data$threshold[all_corr_data$threshold == '025'] <- '0.25'
all_corr_data$threshold[all_corr_data$threshold == '030'] <- '0.30'
all_corr_data$threshold[all_corr_data$threshold == '035'] <- '0.35'
all_corr_data$threshold[all_corr_data$threshold == '040'] <- '0.40'
all_corr_data$threshold[all_corr_data$threshold == '045'] <- '0.45'
all_corr_data$threshold <- as.numeric(as.character(all_corr_data$threshold))


#### GNGs ####
### Modularity Plots ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))
    
# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_mod_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Global Efficiency Plots ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ge_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Participation Coefficient Plots ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Node Dissociation Index ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Clustering Coefficient ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGs_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGs_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))



#### GNGr ####
### Modularity Plots ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_mod_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Global Efficiency Plots ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ge_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Participation Coefficient Plots ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Node Dissociation Index ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Clustering Coefficient ###
## Frontoparietal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


# Somatomotor Dorsal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'GNGr_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'GNGr_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


#### Difference Scores ####
### Modularity Plots ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_mod_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_mu_TD', threshold, NA), ifelse(comparison == 'diff_mod_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_mu_TD', threshold, NA), ifelse(comparison == 'diff_mod_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_sigma_TD', threshold, NA), ifelse(comparison == 'diff_mod_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_sigma_TD', threshold, NA), ifelse(comparison == 'diff_mod_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_tau_TD', threshold, NA), ifelse(comparison == 'diff_mod_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_tau_TD', threshold, NA), ifelse(comparison == 'diff_mod_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_commission_TD', threshold, NA), ifelse(comparison == 'diff_mod_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_commission_TD', threshold, NA), ifelse(comparison == 'diff_mod_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_mod_omission_TD', threshold, NA), ifelse(comparison == 'diff_mod_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_omission_TD', threshold, NA), ifelse(comparison == 'diff_mod_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_mod_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_mod_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Global Efficiency ###
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ge_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_mu_TD', threshold, NA), ifelse(comparison == 'diff_ge_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_mu_TD', threshold, NA), ifelse(comparison == 'diff_ge_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ge_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ge_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_tau_TD', threshold, NA), ifelse(comparison == 'diff_ge_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_tau_TD', threshold, NA), ifelse(comparison == 'diff_ge_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_commission_TD', threshold, NA), ifelse(comparison == 'diff_ge_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_commission_TD', threshold, NA), ifelse(comparison == 'diff_ge_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ge_omission_TD', threshold, NA), ifelse(comparison == 'diff_ge_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_omission_TD', threshold, NA), ifelse(comparison == 'diff_ge_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ge_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ge_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Participation Coefficient ###
## Frontoparietal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Somatomotor Dorsal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_pc_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_pc_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_pc_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Node Dissociation Index ###
## Frontoparietal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Somatomotor Dorsal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_ndi_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_ndi_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


### Clustering Coefficient ###
## Frontoparietal Network #
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_fpn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_fpn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Cingulo-Opercular Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_con_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_con_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Default Mode Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_dmn_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_dmn_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Somatomotor Dorsal Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_smd_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_smd_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))


## Visual Network ##
# Mean RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_meanRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_meanRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_meanRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_meanRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# SD RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_SDRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_SDRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_SDRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_SDRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# CV RT #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_CVRT_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_CVRT_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_CVRT_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_CVRT_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Mu #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_mu_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_mu_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_mu_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_mu_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Sigma #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_sigma_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_sigma_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_sigma_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_sigma_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Tau #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_tau_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_tau_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_tau_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_tau_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Commission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_commission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_commission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_commission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_commission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

# Omission Error Rate #
ggplot(all_corr_data) +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_omission_TD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_omission_TD', r.value, NA), colour = 'TD'), method = 'lm') +
  geom_point(aes(ifelse(comparison == 'diff_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), size = 0.8) +
  geom_smooth(aes(ifelse(comparison == 'diff_clustcoeff_visual_omission_ADHD', threshold, NA), ifelse(comparison == 'diff_clustcoeff_visual_omission_ADHD', r.value, NA), colour = 'ADHD'), method = 'lm') +
  scale_colour_manual(name = 'Diagnosis', values = c(ADHD = 'red', TD = 'blue')) +
  xlab('Threshold') + ylab('Partial Correlation') +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45))

