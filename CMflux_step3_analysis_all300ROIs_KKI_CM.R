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
#install.packages("lmerTest")

library(rstudioapi)
#current_path <- getActiveDocumentContext()$path 
#setwd(dirname(current_path))
print( getwd() )

## NOTE: INCLUDE DENSITY IN THE MODELS, COMPARE NDI AND PC

######## look at differences in whole brain metrics by threshold ########

## load in graph metrics modules

datthreshfiles = list.files("./thresh_testing_all300ROIs/graph_metrics_thresh015/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE)
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh020/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh025/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh030/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh035/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh040/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))
datthreshfiles = append(datthreshfiles, list.files("./thresh_testing_all300ROIs/graph_metrics_thresh045/", pattern=glob2rx("CMflux_wholebrain_metrics_03.16.2020.csv"), full.names = TRUE))

names = sapply(datthreshfiles, function (s) strsplit(s, '//')[[1]][2])
#sess = sapply(datthreshfiles, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][2])
thresh = sapply(datthreshfiles, function (s) strsplit(strsplit(s, '//')[[1]][1], '_')[[1]][5])

#### look at modularity across different thresholds
dat_thresh = data.frame(thresh = NA, 
                 task = NA,
                 N_mod = NA,
                 mean_mod = NA,
                 sd_mod = NA,
                 se_mod = NA,
                 N_density = NA,
                 mean_density = NA,
                 sd_density = NA,
                 se_density = NA,
                 N_ge = NA,
                 mean_ge = NA,
                 sd_ge = NA,
                 se_ge = NA)

# loop through all thresh csv files; pull out mean modularity, density, and global efficiency and paste them into a new dataframe with info from each threshold
for (i in 1:length(datthreshfiles)){
  print(paste0(i,' ', names[i]))
  threshmetrics=read.csv(datthreshfiles[i], sep=",", stringsAsFactors = F)
  
  temp_datthresh <- ddply(threshmetrics, c("task"), summarise,
                 N_mod    = length(mod),
                 mean_mod = mean(mod, na.rm=T),
                 sd_mod   = sd(mod, na.rm=T),
                 se_mod   = sd_mod / sqrt(N_mod),
                 N_density    = length(density),
                 mean_density = mean(density, na.rm=T),
                 sd_density  = sd(density, na.rm=T),
                 se_density   = sd_density / sqrt(N_density),
                 N_ge    = length(ge),
                 mean_ge = mean(ge, na.rm=T),
                 sd_ge  = sd(ge, na.rm=T),
                 se_ge   = sd_ge / sqrt(N_ge))
  
  temp_datthresh$thresh=as.character(thresh[i])
  dat_thresh = rbind(dat_thresh,temp_datthresh)
}

save(dat_thresh, file = "CMFlux_thresh_comparison_all300ROIs.RData")

ggplot(dat_thresh[!is.na(dat_thresh$task),], aes(x = thresh, y = mean_mod)) + geom_point(aes(color = task)) + geom_line(aes(group = task, color = task)) + geom_errorbar(aes(ymin=mean_mod-se_mod, ymax=mean_mod+se_mod, color = task), width=.1) + theme_classic()+ theme(axis.line.x = element_line(), axis.line.y = element_line(), axis.title.y=element_text(size = rel(1.0), vjust=0.7), axis.text.y = element_text(size = rel(0.75), vjust=0.4), axis.title.x = element_text(size = rel(1.0), vjust=0.3), axis.text.x = element_text(angle = 30, size = rel(0.7), vjust=1, hjust=1), title= element_text(size = rel(0.75), vjust=0.4), legend.position= "right") + xlab("Threshold") + ylab("Average Modularity") + ggtitle("") 

ggplot(dat_thresh[!is.na(dat_thresh$task),], aes(x = thresh, y = mean_density)) + geom_point(aes(color = task)) + geom_line(aes(group = task, color = task)) + geom_errorbar(aes(ymin=mean_density-se_density, ymax=mean_density+se_density, color = task), width=.1) + theme_classic()+ theme(axis.line.x = element_line(), axis.line.y = element_line(), axis.title.y=element_text(size = rel(1.0), vjust=0.7), axis.text.y = element_text(size = rel(0.75), vjust=0.4), axis.title.x = element_text(size = rel(1.0), vjust=0.3), axis.text.x = element_text(angle = 30, size = rel(0.7), vjust=1, hjust=1), title= element_text(size = rel(0.75), vjust=0.4), legend.position= "right") + xlab("Threshold") + ylab("Average Density") + ggtitle("") 

ggplot(dat_thresh[!is.na(dat_thresh$task),], aes(x = thresh, y = mean_ge)) + geom_point(aes(color = task)) + geom_line(aes(group = task, color = task)) + geom_errorbar(aes(ymin=mean_ge-se_ge, ymax=mean_ge+se_ge, color = task), width=.1) + theme_classic() + theme(axis.line.x = element_line(), axis.line.y = element_line(), axis.title.y=element_text(size = rel(1.0), vjust=0.7), axis.text.y = element_text(size = rel(0.75), vjust=0.4), axis.title.x = element_text(size = rel(1.0), vjust=0.3), axis.text.x = element_text(angle = 30, size = rel(0.7), vjust=1, hjust=1), title= element_text(size = rel(0.75), vjust=0.4), legend.position= "right") + xlab("Threshold") + ylab("Average Global Efficiency") + ggtitle("") 





####### decide what threshold to use, then run the rest of this #######

## using 0.45 for now 
dat=read.csv("./thresh_testing_all300ROIs/graph_metrics_thresh045/CMflux_wholebrain_metrics_03.16.2020.csv", sep=",", stringsAsFactors = F)
# change the names of the tasks so that it uses GNGs as the baseline comparison in the lmer models below
dat$task[dat$task == 'task-GNGr'] = 'task3-GNGr'
dat$task[dat$task == 'task-GNGs'] = 'task1-GNGs'
dat$task[dat$task == 'task-rest'] = 'task2-rest'

diag_df <- read.csv('subject_demographics.csv', sep = ',', stringsAsFactors = FALSE)
dat <- merge(dat, diag_df, by = ('sub'))

  
###### Whole-Brain Metrics ######
## this multilevel model (more like an ANOVA (than a standard linear regression) but accounting for variance across cell means and individual level variance (within-subject; RE = random effects))
## checking for differences in conditions, not predicting a continuous relationship
## when t value is neg, the condition you are comparing to the baseline condition (GNGs) is less than the baseline (i.e., GNGr less than GNGs)
#### whenever you run a t-test, whatever you put first the sign on the t-statistic tells you the direction of the relationship (i.e. if t value is positive, the first item is greater than the second item) // this would make it seem that in this model the item put first is the comparison group and the baseline group (e.g., GNGs in these models) is put second so you are saying is this one greater than or less than GNGs and if the t value is positive it is greater
##t.test(dat$ge[dat$task == "task-GNGs"],dat$ge[dat$task == "task-GNGs"])
## report in text: beta, p
## report in table: beta, se, t, p


### H1: higher global efficiency and lower modularity from Rest -> GNGs -> GNGr
### H1a: lower global efficiency and higher modularity in ADHD relative to TD across task states (e.g., Lin et al., 2014) (and lower relative changes across tasks)

# Global Efficiency
# asking: are there differences in global efficiency (ge) between tasks and ADHD status controlling for density and the random effect of within-subject variance?
summary(lmer(ge ~ task + density + (1|sub), data = dat))
  # rest marginally higher GE than GNGs
  # no significant difference between GNGs and GNGr

summary(lmer(ge ~ task + diag + density + (1|sub) + task*diag, data = dat))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction


# Modularity
summary(lmer(mod ~ task + density + (1|sub), data = dat))
  # no significant difference between rest and GNGs, or between GNGs and GNGr

summary(lmer(mod ~ task + diag + density + (1|sub) + task*diag, data = dat))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction


# Density - this isn't a variable of most interest; included to be able to control for differences between thresholding
#summary(lmer(density ~ task +(1|sub), data = dat))





###### Nodal Metrics within Networks ######

# pull density column out of the whole-brain dataframe
require(dplyr)
dat_dens <- dat %>% select(sub, task, density)
# merge density column into each nodal metric dataframe (in each section below)


#### Local Efficiency ####
le_netavg = read.csv("./thresh_testing_all300ROIs/graph_metrics_thresh030/CMflux_localeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
le_netavg$task[le_netavg$task == 'task-GNGr'] = 'task3-GNGr'
le_netavg$task[le_netavg$task == 'task-GNGs'] = 'task1-GNGs'
le_netavg$task[le_netavg$task == 'task-rest'] = 'task2-rest'
le_netavg_Wdens = merge(le_netavg, dat_dens, by = c("sub", "task"))
le_netavg_Wdiag = merge(le_netavg_Wdens, diag_df, by = ('sub'))

# Overall Local Efficiency (Whole-Brain Metric)
summary(lmer(overall_le ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in rest than GNGs
  # no significant difference in local efficiency between GNGs and GNGr

summary(lmer(overall_le ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # marginally higher local efficiency in rest than GNGs
  # no significant difference in local efficiency between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Frontoparietal Network
summary(lmer(avg_le_net5 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # no significant difference between rest and GNGs, or between GNGs and GNGr

summary(lmer(avg_le_net5 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Dorsal Attention Network
summary(lmer(avg_le_net4 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in DAN in rest than GNGs, but no difference between GNGs and GNGr

summary(lmer(avg_le_net4 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Cingulo-Opercular Network
summary(lmer(avg_le_net2 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # no significant difference between rest and GNGs, or between GNGs and GNGr

summary(lmer(avg_le_net2 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Default Mode Network
summary(lmer(avg_le_net3 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in DMN in rest than GNGs, but no difference between GNGs and GNGr

summary(lmer(avg_le_net3 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # no significant difference between rest and GNGs, or between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Somatomotor Dorsal Network
summary(lmer(avg_le_net10 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in SMDN in rest than GNGs, but no difference between GNGs and GNGr

summary(lmer(avg_le_net10 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # higher local efficiency in SMDN in rest than GNGs, but no difference between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Somatomotor Ventral Network
summary(lmer(avg_le_net11 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in SMVN in rest than GNGs, but no difference between GNGs and GNGr

summary(lmer(avg_le_net11 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # higher local efficiency in SMVN in rest than GNGs, but no difference between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction

# Visual Network
summary(lmer(avg_le_net14 ~ task + density +(1|sub), data = le_netavg_Wdens))
  # higher local efficiency in VN in rest than GNGs, but no difference between GNGs and GNGr

summary(lmer(avg_le_net14 ~ task + diag + density +(1|sub) + task*diag, data = le_netavg_Wdiag))
  # marginally higher local efficiency in VN in rest than GNGs, but no difference between GNGs and GNGr
  # no significant effect of diagnosis, or task x diagnosis interaction




#### Participation Coefficient ####
pc_netavg = read.csv("./thresh_testing_all300ROIs/graph_metrics_thresh030/CMflux_particcoeffsNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
pc_netavg$task[pc_netavg$task == 'task-GNGr'] = 'task3-GNGr'
pc_netavg$task[pc_netavg$task == 'task-GNGs'] = 'task1-GNGs'
pc_netavg$task[pc_netavg$task == 'task-rest'] = 'task2-rest'
pc_netavg_Wdens = merge(pc_netavg, dat_dens, by = c("sub", "task"))
pc_netavg_Wdiag = merge(pc_netavg_Wdens, diag_df, by = ('sub'))

# Frontoparietal Network
summary(lmer(avg_pc_net5 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net5 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Dorsal Attention Network
summary(lmer(avg_pc_net4 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net4 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Cingulo-Opercular Network
summary(lmer(avg_pc_net2 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net2 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Default Mode Network
summary(lmer(avg_pc_net3 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net3 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Somatomotor Dorsal Network
summary(lmer(avg_pc_net10 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net10 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Somatomotor Ventral Network
summary(lmer(avg_pc_net11 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net11 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))

# Visual Network
summary(lmer(avg_pc_net14 ~ task + density +(1|sub), data = pc_netavg_Wdens))
summary(lmer(avg_pc_net14 ~ task + diag + density +(1|sub) + task*diag, data = pc_netavg_Wdiag))




#### Node Dissociation Index ####
ndi_netavg = read.csv("./thresh_testing_all300ROIs/graph_metrics_thresh030/CMflux_ndisNETAVG_03.16.2020.csv", sep=",", stringsAsFactors = F)
ndi_netavg$task[ndi_netavg$task == 'task-GNGr'] = 'task3-GNGr'
ndi_netavg$task[ndi_netavg$task == 'task-GNGs'] = 'task1-GNGs'
ndi_netavg$task[ndi_netavg$task == 'task-rest'] = 'task2-rest'
ndi_netavg_Wdens = merge(ndi_netavg, dat_dens, by = c("sub", "task"))
ndi_netavg_Wdiag = merge(ndi_netavg_Wdens, diag_df, by = ('sub'))

# Frontoparietal Network
summary(lmer(avg_ndi_net5 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net5 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Dorsal Attention Network
summary(lmer(avg_ndi_net4 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net4 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Cingulo-Opercular Network
summary(lmer(avg_ndi_net2 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net2 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Default Mode Network
summary(lmer(avg_ndi_net3 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net3 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Somatomotor Dorsal Network
summary(lmer(avg_ndi_net10 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net10 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Somatomotor Ventral Network
summary(lmer(avg_ndi_net11 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net11 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Visual Network
summary(lmer(avg_ndi_net14 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net14 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))




#### within module degree - 02/25/2020, cannot look at network differences because z-scored within each network ####

## could also do this without z-scoring; would need to write the function







#### Clustering Coefficient (AKA Transitivity) ####
clustcoeff_netavg = read.csv("./thresh_testing_all300ROIs/graph_metrics_thresh030/CMflux_netclustcoeff_03.16.2020.csv", sep = ',', stringsAsFactors = FALSE)
clustcoeff_netavg$task[clustcoeff_netavg$task == 'task-GNGr'] = 'task3-GNGr'
clustcoeff_netavg$task[clustcoeff_netavg$task == 'task-GNGs'] = 'task1-GNGs'
clustcoeff_netavg$task[clustcoeff_netavg$task == 'task-rest'] = 'task2-rest'
clustcoeff_netavg_Wdens = merge(clustcoeff_netavg, dat_dens, by = c("sub", "task"))
clustcoeff_netavg_Wdiag = merge(clustcoeff_netavg_Wdens, diag_df, by = ('sub'))

##HELP: DO WE CONTROL FOR DENSITY OR NOT?? i am thinking not because clustcoeff is a percentage
##HELP: DO WE NEED TO DO THIS REGRESSION MODEL OR CAN WE DO T.TESTS?

# Frontoparietal Network
summary(lmer(clustcoeff_fpn ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_fpn ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))

# Dorsal Attention Network - CALCULATE THEM
summary(lmer(avg_ndi_net4 ~ task + density +(1|sub), data = ndi_netavg_Wdens))
summary(lmer(avg_ndi_net4 ~ task + diag + density +(1|sub) + task*diag, data = ndi_netavg_Wdiag))

# Cingulo-Opercular Network
summary(lmer(clustcoeff_con ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_con ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))

# Default Mode Network
summary(lmer(clustcoeff_dmn ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_dmn ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))

# Somatomotor Dorsal Network
summary(lmer(clustcoeff_smd ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_smd ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))

# Somatomotor Ventral Network
summary(lmer(clustcoeff_smv ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_smv ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))

# Visual Network
summary(lmer(clustcoeff_visual ~ task + density +(1|sub), data = clustcoeff_netavg_Wdens))
summary(lmer(clustcoeff_visual ~ task + diag + density +(1|sub) + task*diag, data = clustcoeff_netavg_Wdiag))






#### FDR adjust for multiple comparisons ####
## THERE'S GOTTA BE A BETTER WAYYYY - but hardcoded p values from lmer analyses
# make vector of all p-values ##NOTE: from talking with TRH on 3/4/2020, do this in chunks by question and graph metric and do not correct the whole brain metrics
aim1pc_ps = list(0.03889,
                 0.26763,
                 1.85E-06,
                 0.0845,
                 0.013,
                 0.04107)
aim1ndi_ps = list(0.0516,
                  0.22882,
                  1.57E-06,
                  0.536,
                  0.0232,
                  0.677)
aim1cc_ps = list(0.0814,
                 0.136,
                 0.1755,
                 0.932,
                 0.0465,
                 0.339)

aim2pc_ps = list(0.00048,
                 0.00526,
                 0.106,
                 0.022,
                 8.97E-05,
                 0.00872)
aim2ndi_ps = list(0.0378,
                  0.09849,
                  0.1526,
                  0.295,
                  0.035,
                  0.11)
aim2cc_ps = list(0.0065,
                 0.156,
                 0.0676,
                 0.492,
                 0.7712,
                 0.618)
# create vector of adjusted p values
adj_ps_aim1pc = p.adjust(aim1pc_ps, method = "fdr", n = length(aim1pc_ps))
adj_ps_aim1ndi = p.adjust(aim1ndi_ps, method = "fdr", n = length(aim1ndi_ps))
adj_ps_aim1cc = p.adjust(aim1cc_ps, method = "fdr", n = length(aim1cc_ps))

adj_ps_aim2pc = p.adjust(aim2pc_ps, method = "fdr", n = length(aim2pc_ps))
adj_ps_aim2ndi = p.adjust(aim2ndi_ps, method = "fdr", n = length(aim2ndi_ps))
adj_ps_aim2cc = p.adjust(aim2cc_ps, method = "fdr", n = length(aim2cc_ps))
## manually enter adj_ps into google sheet tracking results





######## calculate provincial and connector hubs ########
## ref: Guimera & Nunes Amaral_2005Nature
#### Provincial hubs: PC < 0.3 + WMD > 2.5
#### Connector hubs: PC [0.3 - 0.7] + WMD > 2.5
## rich club coefficients?









