#### CM Flux Analysis ####
## April 2020, CM ##
## This script takes correlation matrices indicating the functional connectivity between pairs of ROIs in the 5 core networks (FPN, CON, DMN, Somat Dorsal, Visual) to generate correlation plots ##
# This is because the analyses of the clustering coefficient for the CON are underpowered, so we want to visualize the FC plots to check whether there is anything weird in the CON ROIs

#network 1 = CON
#network 2 = DMN
#network 3 = FPN
#network 4 = SMD
#network 5 = VIS


rm(list=ls())

#install.packages("corrplot")
#install.packages("ggcorrplot")
#install.packages("popbio")
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(ggcorrplot)
library(popbio) # need this library to average across a list of matrices


#### A. Generate functional connectivity matrices for each task (rest, GNGs, GNGr) averaging across all subjects ####

## A1. Load in all the correlation matrices created in "Recentering & Rescaling ROI Time Series Across Runs" from ROI timeseries ##
# These will be used to plot correlation matrices for each task for each subject
corrmat_files <- list.files("./correlation_matrices_5CoreNetworks/", pattern = glob2rx("sub-*.csv"), full.names = TRUE)
#names <- sapply(corrmat_files, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
#sess = sapply(corrmat_files, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][2])
#tasks = sapply(corrmat_files, function (s) strsplit(strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][3], '.csv' [[1]][1]))


## A2. Load in all the correlation matrices created in "Recentering & Rescaling ROI Time Series Across Runs" from ROI timeseries for the resting-state condition, and construct a single correlation matrix averaging across subjects for the resting-state condition ##

corrmat_files_rest <- list.files("./correlation_matrices_5CoreNetworks/", pattern = glob2rx("sub-*rest.csv"), full.names = TRUE)  # read in only the GNGs correlation matrices
#list_matrices_rest <- lapply(lapply(lapply(corrmat_files_rest, read.csv), as.matrix), function(x){x[which(x < 0.00, arr.ind = TRUE)] = 0
#return(x)})  # read in individual subject's csv files, make them matrices, and then remove negative values from the individual subject matrices
list_matrices_rest <- lapply(lapply(corrmat_files_rest, read.csv), as.matrix)

avg_corrmat_rest <- mean.list(list_matrices_rest, na.rm = TRUE)
#image(avg_corrmat_rest)
ggcorrplot(avg_corrmat_rest, outline.color = 'white')
#avg_corrmat_rest[which(avg_corrmat_rest < 0.10, arr.ind = TRUE)] = 0  # absolute thresholding of the group average matrix
#image(avg_corrmat_rest)


## A3. Load in all the correlation matrices created in "Recentering & Rescaling ROI Time Series Across Runs" from ROI timeseries for the GNGs task, and construct a single correlation matrix averaging across subjects for the GNGs task ##

corrmat_files_gngs <- list.files("./correlation_matrices_5CoreNetworks/", pattern = glob2rx("sub-*GNGs.csv"), full.names = TRUE)  # read in only the GNGs correlation matrices
#list_matrices_gngs <- lapply(lapply(lapply(corrmat_files_gngs, read.csv), as.matrix), function(x){x[which(x < 0, arr.ind = TRUE)] = 0
#return(x)})  # read in individual subject csv files, make them matrices, and then remove negative values from the individual subject matrices
list_matrices_gngs <- lapply(lapply(corrmat_files_gngs, read.csv), as.matrix)

avg_corrmat_gngs <- mean.list(list_matrices_gngs, na.rm = TRUE)
#image(avg_corrmat_gngs)
ggcorrplot(avg_corrmat_gngs, outline.color = 'white')
#avg_corrmat_gngs[which(avg_corrmat_gngs < 0.1, arr.ind = TRUE)] = 0  # absolute thresholding of the group average matrix
#image(avg_corrmat_gngs)


## A4. Load in all the correlation matrices created in "Recentering & Rescaling ROI Time Series Across Runs" from ROI timeseries for the GNGr task, and construct a single correlation matrix averaging across subjects for the GNGr task ##

corrmat_files_gngr <- list.files("./correlation_matrices_5CoreNetworks/", pattern = glob2rx("sub-*GNGr.csv"), full.names = TRUE)  # read in only the GNGr correlation matrices
#list_matrices_gngr <- lapply(lapply(lapply(corrmat_files_gngr, read.csv), as.matrix), function(x){x[which(x < 0, arr.ind = TRUE)] = 0
#return(x)})  # read in individual subject csv files, make them matrices, and then remove negative values from the individual subject matrices
list_matrices_gngr <- lapply(lapply(corrmat_files_gngr, read.csv), as.matrix)

avg_corrmat_gngr <- mean.list(list_matrices_gngr, na.rm = TRUE)
#image(avg_corrmat_gngr)
ggcorrplot(avg_corrmat_gngr, outline.color = 'white')
#avg_corrmat_gngr[which(avg_corrmat_gngr < 0.1, arr.ind = TRUE)] = 0  # absolute thresholding of the group average matrix
#image(avg_corrmat_gngr)



#### B. Generate functional connectivity matrices for each task (rest, GNGs, GNGr) averaging across all subjects - within diagnosis (i.e., separate FC matrices for ADHD and TD subjects for each task) ####
sub_diag <- read.csv('subject_demographics.csv', sep = ',')

## B1. Resting-State Condition ##
corrmat_files_rest_adhd <- c()
corrmat_files_rest_td <- c()

for (i in 1:length(corrmat_files_rest)) {
  
  corrmat_file <- corrmat_files_rest[i]
  sub_id <- sapply(corrmat_file, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
  
  if (sub_diag$diag[sub_diag$sub == sub_id] == 'ADHD') {
    corrmat_files_rest_adhd <- c(corrmat_files_rest_adhd, corrmat_file)
  } else {
    corrmat_files_rest_td <- c(corrmat_files_rest_td, corrmat_file)
  }
  
} 


list_matrices_rest_adhd <- lapply(lapply(corrmat_files_rest_adhd, read.csv), as.matrix)
avg_corrmat_rest_adhd <- mean.list(list_matrices_rest_adhd, na.rm = TRUE)
ggcorrplot(avg_corrmat_rest_adhd, outline.color = 'white')

list_matrices_rest_td <- lapply(lapply(corrmat_files_rest_td, read.csv), as.matrix)
avg_corrmat_rest_td <- mean.list(list_matrices_rest_td, na.rm = TRUE)
ggcorrplot(avg_corrmat_rest_td, outline.color = 'white')


## B2. GNGs Condition #
corrmat_files_gngs_adhd <- c()
corrmat_files_gngs_td <- c()

for (i in 1:length(corrmat_files_gngs)) {
  
  corrmat_file <- corrmat_files_gngs[i]
  sub_id <- sapply(corrmat_file, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
  
  if (sub_diag$diag[sub_diag$sub == sub_id] == 'ADHD') {
    corrmat_files_gngs_adhd <- c(corrmat_files_gngs_adhd, corrmat_file)
  } else {
    corrmat_files_gngs_td <- c(corrmat_files_gngs_td, corrmat_file)
  }
  
} 


list_matrices_gngs_adhd <- lapply(lapply(corrmat_files_gngs_adhd, read.csv), as.matrix)
avg_corrmat_gngs_adhd <- mean.list(list_matrices_gngs_adhd, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngs_adhd, outline.color = 'white')

list_matrices_gngs_td <- lapply(lapply(corrmat_files_gngs_td, read.csv), as.matrix)
avg_corrmat_gngs_td <- mean.list(list_matrices_gngs_td, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngs_td, outline.color = 'white')


## B3. GNGr Condition #
corrmat_files_gngr_adhd <- c()
corrmat_files_gngr_td <- c()

for (i in 1:length(corrmat_files_gngr)) {
  
  corrmat_file <- corrmat_files_gngr[i]
  sub_id <- sapply(corrmat_file, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
  
  if (sub_diag$diag[sub_diag$sub == sub_id] == 'ADHD') {
    corrmat_files_gngr_adhd <- c(corrmat_files_gngr_adhd, corrmat_file)
  } else {
    corrmat_files_gngr_td <- c(corrmat_files_gngr_td, corrmat_file)
  }
  
} 


list_matrices_gngr_adhd <- lapply(lapply(corrmat_files_gngr_adhd, read.csv), as.matrix)
avg_corrmat_gngr_adhd <- mean.list(list_matrices_gngr_adhd, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngr_adhd, outline.color = 'white')

list_matrices_gngr_td <- lapply(lapply(corrmat_files_gngr_td, read.csv), as.matrix)
avg_corrmat_gngr_td <- mean.list(list_matrices_gngr_td, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngr_td, outline.color = 'white')



#### C. Generate functional connectivity matrices of difference scores (i. GNGs - rest - ii. GNGr - GNGs) averaging across all subjects ####
corrmat_files_rest_subs <- sapply(corrmat_files_rest, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
corrmat_files_gngs_subs <- sapply(corrmat_files_gngs, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
corrmat_files_gngr_subs <- sapply(corrmat_files_gngr, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])


## C1. GNGs - Rest ##
corrmat_files_gngsrest_rest <- c()
corrmat_files_gngsrest_gngs <- c()

# Get the GNGs files for the subjects who have both GNGs and resting-state files
for (i in 1:length(corrmat_files_gngs)) {
  if (corrmat_files_gngs_subs[i] %in% corrmat_files_rest_subs) {
    corrmat_files_gngsrest_gngs <- c(corrmat_files_gngsrest_gngs, corrmat_files_gngs[i])
  }
}

corrmat_files_gngsrest_gngs_subs <- sapply(corrmat_files_gngsrest_gngs, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])


# Get the resting-state files for the subjects who have both GNGs and resting-state files
for (i in 1:length(corrmat_files_rest)) {
  if (corrmat_files_rest_subs[i] %in% corrmat_files_gngsrest_gngs_subs) {
    corrmat_files_gngsrest_rest <- c(corrmat_files_gngsrest_rest, corrmat_files_rest[i])
  }
}


# Create a correlation matrix of GNGs - rest difference scores for each subject with both GNGs and resting-state data #
for (i in 1:length(corrmat_files_gngsrest_gngs)) {
  
  temp_corrmat_file_rest <- read.csv(corrmat_files_gngsrest_rest[i], sep = ',')
  temp_corrmat_file_gngs <- read.csv(corrmat_files_gngsrest_gngs[i], sep = ',')
  
  temp_corrmat_file_restgngs_diff <- temp_corrmat_file_gngs - temp_corrmat_file_rest
  temp_sub_id <- sapply(corrmat_files_gngsrest_gngs[i], function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
  temp_ses <- 'ses-1'
  temp_condition <- 'GNGs-rest'
  file_name <- paste(temp_sub_id, "_", temp_ses, "_", temp_condition, ".csv", sep = "")
  write.csv(temp_corrmat_file_restgngs_diff, file_name, row.names = FALSE)
}


# Generate a functional connectivity matrix for the GNGs-rest difference averaging across all subjects #
corrmat_files_gngsrest <- list.files("./correlation_matrices_GNGs-rest/", full.names = TRUE)
list_matrices_gngsrest <- lapply(lapply(corrmat_files_gngsrest, read.csv), as.matrix)
avg_corrmat_gngsrest <- mean.list(list_matrices_gngsrest, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngsrest, outline.color = 'white')



## C2. GNGr - GNGs ##
corrmat_files_gngrgngs_gngs <- c()
corrmat_files_gngrgngs_gngr <- c()

# Get the GNGs files for the subjects who have both GNGs and GNGr files
for (i in 1:length(corrmat_files_gngs)) {
  if (corrmat_files_gngs_subs[i] %in% corrmat_files_gngr_subs) {
    corrmat_files_gngrgngs_gngs <- c(corrmat_files_gngrgngs_gngs, corrmat_files_gngs[i])
  }
}

corrmat_files_gngrgngs_gngs_subs <- sapply(corrmat_files_gngrgngs_gngs, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])


# Get the GNGs files for the subjects who have both GNGs and GNGr files
for (i in 1:length(corrmat_files_gngr)) {
  if (corrmat_files_gngr_subs[i] %in% corrmat_files_gngrgngs_gngs_subs) {
    corrmat_files_gngrgngs_gngr <- c(corrmat_files_gngrgngs_gngr, corrmat_files_gngr[i])
  }
}


# Create a correlation matrix of GNGr - GNGs difference scores for each subject with both GNGs and resting-state data #
for (i in 1:length(corrmat_files_gngrgngs_gngs)) {
  
  temp_corrmat_file_gngs <- read.csv(corrmat_files_gngrgngs_gngs[i], sep = ',')
  temp_corrmat_file_gngr <- read.csv(corrmat_files_gngrgngs_gngr[i], sep = ',')
  
  temp_corrmat_file_gngrgngs_diff <- temp_corrmat_file_gngr - temp_corrmat_file_gngs
  temp_sub_id <- sapply(corrmat_files_gngrgngs_gngs[i], function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
  temp_ses <- 'ses-1'
  temp_condition <- 'GNGr-GNGs'
  file_name <- paste(temp_sub_id, "_", temp_ses, "_", temp_condition, ".csv", sep = "")
  write.csv(temp_corrmat_file_gngrgngs_diff, file_name, row.names = FALSE)
}


# Generate a functional connectivity matrix for the GNGr-GNGs difference averaging across all subjects #
corrmat_files_gngrgngs <- list.files("./correlation_matrices_GNGr-GNGs/", full.names = TRUE)
list_matrices_gngrgngs <- lapply(lapply(corrmat_files_gngrgngs, read.csv), as.matrix)
avg_corrmat_gngrgngs <- mean.list(list_matrices_gngrgngs, na.rm = TRUE)
ggcorrplot(avg_corrmat_gngrgngs, outline.color = 'white')

