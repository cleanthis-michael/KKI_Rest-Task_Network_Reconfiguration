#### Network Cleaning and Net Stats Calculation ####
## written by TRHenry, adapted for First Year Project 2020 by MEMitchell (February 2020)
## adapted for KKI task analyses for Flux Conference by Cleanthis Michael (March 2020)

## Step 2 in CM Flux analysis - to be completed after the 'ROI Time Series Across Runs.R' script that creates the correlation matrices from each subject scaled and concatenated across task runs
## purpose - to calculate graph metrics (whole brain and nodal) from each correlation matrix and output data.frames (whole brain metrics, nodal pc, nodal wmd, nodal le, nodal ndi, network pc, network ndi, network clust coeff)


# NOTE: 2/25/2020, teague confirmed my thoughts on WMD - averaging within module is useless bc it z-scores within modules. so instead we are doing to make subgraphs (an iGraph function) of each module/network separately and then run global clustering coeff on each (transitivity, an iGraph function). this will give a measure of the percentage of closed triads which functions in this analysis as a measure of integration within each network.

## NOTE: 3/12/2020, this script is step 1 in constructing graphs from correlation matrices with just the five core networks of interest (frontoparietal, cingulo-opercular, default mode, somatomotor dorsal, visual, 214 ROIs total); additionally, the networks are assigned new numbers for this analysis bc there are only 5:
#network 1 = CON
#network 2 = DMN
#network 3 = FPN
#network 4 = SMD
#network 5 = VIS


rm(list=ls())

library(dplyr)
library(plyr)
library(tidyr)
library(igraph)
library(brainGraph)
#library(netcontrol) # not available yet bc teague is building it (will be available on CRAN soon, allows for control theory analyses of brain data)
library(Matrix)
library(progress)

# set wd
library(rstudioapi)
#current_path <- getActiveDocumentContext()$path 
#setwd(dirname(current_path))
print( getwd() )


#### A. prep for loop that makes graphs and pulls graph metrics ####

# A1. specify the network membership of each ROI (**only using 214ROIs in select networks from the bigbrain (AKA extended Power atlas -- Seitzman et al. 2019) for this analysis)
    ## networks of interest: frontoparietal, cingulo-opercular, somatomotor dorsal, visual, and default mode
bigbrainatlas = read.csv("./bigbrain_rois/BigBrain300_ROIinfo.csv", sep=",") # should be created when you run the ROI timeseries extraction script  - kept only the 2 first columns from downloaded csv file
bigbrain_5corenets = bigbrainatlas[bigbrainatlas$Network.Assignment == 'Fronto-Parietal' | bigbrainatlas$Network.Assignment == 'Cingulo-Opercular' | bigbrainatlas$Network.Assignment == 'Somatomotor Dorsal' | bigbrainatlas$Network.Assignment == 'Visual' | bigbrainatlas$Network.Assignment == 'Default Mode',]
roi_membership_5corenets <- as.vector(bigbrain_5corenets[,2]); roi_membership_5corenets
roi_comms_5corenets <- as.numeric(as.factor(roi_membership_5corenets)); roi_comms_5corenets # creates a numeric vector of network membership for each roi




# A2. read in all correlation matrices created in the 'DataProcNetCon_MEMFYP_trhfix.R' script from ROI timecourses
corrmatfiles= list.files("./correlation_matrices_5CoreNetworks/", pattern=glob2rx("sub-*.csv"), full.names = TRUE)
names = sapply(corrmatfiles, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][1])
sess = sapply(corrmatfiles, function (s) strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][2])
tasks = sapply(corrmatfiles, function (s) strsplit(strsplit(strsplit(s, '//')[[1]][2], '_')[[1]][3], '.csv' [[1]][1]))

# A3. make empty data.frame that will be filled in with whole brain stats in the loop below
dat = data.frame(sub = NA, 
                ses = NA, 
                task = NA,
                density = NA, 
                strength = NA, 
                mod = NA, 
                ge = NA)

# A4. make empty data.frames that will be filled in with nodal metrics
roinums <- as.vector(bigbrain_5corenets[,1])
rois <- list()
for (r in 1:length(roinums)){
  roi = paste("V", roinums[r], sep = "")
  print(roi)
  rois <-append(rois,roi)
}
# make df to build PC values into
dat_pc <- data.frame(matrix(nrow=1, ncol=214))
colnames(dat_pc) = rois
# make df to build WMD values into
dat_wmd <- data.frame(matrix(nrow=1, ncol=214))
colnames(dat_wmd) = rois
# make df to build local efficiency (LE) values into
dat_le <- data.frame(matrix(nrow=1, ncol=214))
colnames(dat_le) = rois
# make df to build NDI values into
dat_ndi <- data.frame(matrix(nrow=1, ncol=214))
colnames(dat_ndi) = rois

# A5. load NDI function (TRH built this, February 2020)
NDI <- function(thresh_corr_mat, comm_vec){
  NDI_vec = vector()
  diag(thresh_corr_mat) = 0
  for(i in 1:dim(thresh_corr_mat)[[1]]){
    if(!all(is.na(thresh_corr_mat[i,]))){
      self_comm = comm_vec[i]
      not_self = -which(comm_vec == self_comm)
      deg = sum(thresh_corr_mat[i,], na.rm = T)
      not_self_deg = sum(thresh_corr_mat[i,][not_self], na.rm = T)
      NDI_vec[i] = not_self_deg/deg
    }else{
      NDI_vec[i] = NA
    }
  }
  return(NDI_vec)
}

# A6. make empty data.frame for clustering coefficient + make subgraph ROI vectors
## note: calculating global clustering coefficient of each network by making subgraphs of each network separately
# make empty data.frame to populate with network clustering coefficients
dat_clustcoeff = data.frame(sub = NA, 
                            ses = NA, 
                            task = NA,
                            clustcoeff_fpn = NA, 
                            clustcoeff_con = NA, 
                            clustcoeff_dmn = NA, 
                            clustcoeff_smd = NA,
                            clustcoeff_visual = NA)

# make ROI vectors that you will need in the subgraph function later
# look at options for network assignments in your preassigned atlas
unique(bigbrain_5corenets$Network.Assignment) 
# add a V in front of each ROI so that we can merge smoothly later
bigbrain_5corenets$ROI.Number2 <- paste("V", bigbrain_5corenets$ROI.Number, sep = "")

# make ROI vectors for each network (MEM commented out the networks that she is not interested in)
rois.frontoparietal <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Fronto-Parietal"]
rois.cinguloopercular <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Cingulo-Opercular"]
rois.dmn <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Default Mode"]
rois.somatdorsal <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Somatomotor Dorsal"]
rois.visual <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Visual"]

#rois.dorsalattn <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Dorsal Attention"]
#rois.somatventral <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Somatomotor Ventral"]
#rois.salience <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Salience"]
#rois.ventralattn <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Ventral Attention"]
#rois.reward <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Reward"]
#rois.unassigned <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Unassigned"]
#rois.medtemp <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Medial Temporal"]
#rois.auditory <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Auditory"]
#rois.parietalmemory <- bigbrain_5corenets$ROI.Number2[bigbrain_5corenets$Network.Assignment == "Parietal Memory"]







#### B. loop through each subject task file, create graphs, calculate graph metrics ####

# B1. loop through each subject task file, load in the correlation matrices, threshold, create igraph object, calculate whole brain graph metrics from it, and add a line to each of the data.frames created above (dat for whole brain metrics, dat_pc/dat_wmd/dat_ndi/dat_le for nodal metrics)
for (i in 1:length(corrmatfiles)){
  print(paste0(i,' ', names[i]))
  
  ###read in correlation matrix
  subcorr=read.csv(corrmatfiles[i], sep=",", stringsAsFactors = F)
  subcorr_mat <- as.matrix(subcorr)

  ###remove negative correlations + threshold **THE UNIVERSAL QUESTION: WHAT THRESHOLD SHOULD I USE?** 
  #CM: run at 0.15 - 0.45 in steps of 0.05 to graph
  subcorr_mat[which(subcorr_mat < 0.10, arr.ind = T)] = 0 #absolute thresholding
  subcorr_mat[which(is.na(subcorr_mat), arr.ind = T)] = 0 #treating ROIs that are outside the field of view (NAs) as an isolate
  hist(subcorr_mat)
  #subcorr_mat[which(subcorr_mat > 0.0,arr.ind = T)] = 1 #binarize the graph (if you want, not binarizing for MEMFYP)
  #subcorr_mat_prop = 1*(subcorr_mat>quantile(subcorr_mat[upper.tri(subcorr_mat)],probs=1-0.1)) #cost thresholding/proportional thresholding, where: probs=1-(percentage you want as your density)
  
  ###create igraph object (aka the graph) **note: setting 'diag = F' is crucial and not the default setting
  corr_nets <- graph.adjacency(subcorr_mat, mode = "undirected", diag = F, weighted = T) ## want a weighted graph (not a binarized one); this line makes an igraph object out of the already thresholded subject correlation matrix (makes igraphs from adjacency matrices = square matrices used to represent a finite graph)
  
  ###first, whole brain metrics
  dat.loop = data.frame(sub = NA, 
                   ses = NA, 
                   task = NA,
                   density = NA, 
                   strength = NA, 
                   mod = NA, 
                   ge = NA)
  dat.loop$sub=as.character(names[i])
  dat.loop$ses=as.character(sess[i])
  dat.loop$task=as.character(tasks[i])
  dat.loop$density = edge_density(corr_nets) # the ratio of the number of edges to the number of possible edges
  dat.loop$strength = mean(strength(graph = corr_nets), na.rm=T) # average sum up edge weights for each ROI
  dat.loop$mod = modularity(corr_nets, roi_comms_5corenets) # looking at modularity when utilizing the big brain atlas network membership assignments (preassigned network membership from adult atlas Power et al. 2011 & extended with additional subcortical and cerebellar coverage in Seitzman et al. 2019)
  dat.loop$ge = efficiency(g = corr_nets, type = 'global', weights = E(corr_nets)$weight) # this is a braingraph function; using a weighted graph ((if using a weighted graph: weights = E(corr_nets)$weight))
  dat=rbind(dat,dat.loop)
  

  ###now, nodal metrics
  ###calculate local efficiency for each node
  localeff = efficiency(g = corr_nets, type = 'local') #check this line, weights made it fail last time; 	"Numeric vector of edge weights; if NULL (the default), and if the graph has edge attribute weight, then that will be used."
  temp_le <- as.data.frame(localeff)
  temp_le_wide <- as.data.frame(t(temp_le))
  colnames(temp_le_wide) = rois 
  temp_le_wide$sub=as.character(names[i])
  temp_le_wide$ses=as.character(sess[i])
  temp_le_wide$task=as.character(tasks[i])
  dat_le=merge(dat_le,temp_le_wide, all=TRUE)
  
  ###calculate participation coefficient (PC) for each node
    # PC iterates through each network and asks each node "how much do you participate in this network? and what about this next one? i dont care which one you think you belong to."
    # PC doesnt care about which network each node belongs to; just how it participates in each network of the system (in this case, 14 networks)
  pc = part_coeff(corr_nets, roi_comms_5corenets) # feed it the igraph object and the vector of roi network assignments
  temp_pc <- as.data.frame(pc)
  temp_pc_wide <- as.data.frame(t(temp_pc))
  colnames(temp_pc_wide) = rois 
  temp_pc_wide$sub=as.character(names[i])
  temp_pc_wide$ses=as.character(sess[i])
  temp_pc_wide$task=as.character(tasks[i])
  dat_pc=merge(dat_pc,temp_pc_wide, all=TRUE)
  
  ###calculate node dissociation index (NDI) for each node
    ## NDI: an adaptation of PC; considers the node's network membership (Cary et al., 2017)
    ## NDI is better for preassigned atlases than PC (TRH)
  ndi = NDI(subcorr_mat, roi_comms_5corenets) #feed it the correlation matrix, not the graph; TRH wrote the function for this
  temp_ndi <- as.data.frame(ndi)
  temp_ndi_wide <- as.data.frame(t(temp_ndi))
  colnames(temp_ndi_wide) = rois 
  temp_ndi_wide$sub=as.character(names[i])
  temp_ndi_wide$ses=as.character(sess[i])
  temp_ndi_wide$task=as.character(tasks[i])
  dat_ndi=merge(dat_ndi,temp_ndi_wide, all=TRUE)
  
  ###calculate within-module degree (WMD) for each node
  wmd = within_module_deg_z_score(corr_nets, roi_comms_5corenets) # feed it the igraph object and the vector of roi network assignments
  temp_wmd <- as.data.frame(wmd)
  temp_wmd_wide <- as.data.frame(t(temp_wmd))
  colnames(temp_wmd_wide) = rois 
  temp_wmd_wide$sub=as.character(names[i])
  temp_wmd_wide$ses=as.character(sess[i])
  temp_wmd_wide$task=as.character(tasks[i])
  dat_wmd=merge(dat_wmd,temp_wmd_wide, all=TRUE)
  
  
  ###calculate clustering coefficients for each network by (1) making network specific subgraph, then (2) calculating global clustering coefficient with the transitivity() function
  # default mode clustering coefficient
  dmn_subgraph = induced_subgraph(graph = corr_nets, vids = rois.dmn, impl = c("auto"))
  transit_dmn <- transitivity(g = dmn_subgraph, type = c("global"), vids = NULL, weights = NULL, isolates = c("NaN", "zero"))
  
  # frontoparietal clustering coefficient
  fpn_subgraph = induced_subgraph(graph = corr_nets, vids = rois.frontoparietal, impl = c("auto"))
  transit_fpn <- transitivity(g = fpn_subgraph, type = c("global"), vids = NULL, weights = NULL, isolates = c("NaN", "zero"))
  
  # cingulo-opercular clustering coefficient
  con_subgraph = induced_subgraph(graph = corr_nets, vids = rois.cinguloopercular, impl = c("auto"))
  transit_con <- transitivity(g = con_subgraph, type = c("global"), vids = NULL, weights = NULL, isolates = c("NaN", "zero"))
  
  # somatomotor dorsal clustering coefficient
  smd_subgraph = induced_subgraph(graph = corr_nets, vids = rois.somatdorsal, impl = c("auto"))
  transit_smd <- transitivity(g = smd_subgraph, type = c("global"), vids = NULL, weights = NULL, isolates = c("NaN", "zero"))
  
  #  visual clustering coefficient
  visual_subgraph = induced_subgraph(graph = corr_nets, vids = rois.visual, impl = c("auto"))
  transit_visual <- transitivity(g = visual_subgraph, type = c("global"), vids = NULL, weights = NULL, isolates = c("NaN", "zero"))
  
  # now, build network clustering coefficients into a data.frame
  temp_clustcoeff = data.frame(sub = NA, 
                               ses = NA, 
                               task = NA,
                               clustcoeff_fpn = NA, 
                               clustcoeff_con = NA, 
                               clustcoeff_dmn = NA, 
                               clustcoeff_smd = NA,
                               clustcoeff_visual = NA)
  temp_clustcoeff$sub=as.character(names[i])
  temp_clustcoeff$ses=as.character(sess[i])
  temp_clustcoeff$task=as.character(tasks[i])
  temp_clustcoeff$clustcoeff_fpn =  transit_fpn
  temp_clustcoeff$clustcoeff_con =  transit_con
  temp_clustcoeff$clustcoeff_dmn =  transit_dmn
  temp_clustcoeff$clustcoeff_smd = transit_smd
  temp_clustcoeff$clustcoeff_visual = transit_visual
  
  dat_clustcoeff=rbind(dat_clustcoeff,temp_clustcoeff)
}


## B2. clean up the nodal metric data.frames

require(dplyr)
dat_pc <- dat_pc %>% select(sub, ses, task, V13:V300)
dat_pc[dat_pc == 0.0000000] <- NaN # set zeros to NA

dat_wmd <- dat_wmd %>% select(sub, ses, task, V13:V300)
dat_wmd[dat_wmd == 0.0000000] <- NaN # set zeros to NA

dat_le <- dat_le %>% select(sub, ses, task, V13:V300)

dat_ndi <- dat_ndi %>% select(sub, ses, task, V13:V300)
dat_ndi[dat_ndi == 0.0000000] <- NA

## B3. and now we save everything, but you have two options to save

# option 1: as csvs
### update names with (i) binarized or weighted, (ii) the threshold type and level
date=format(Sys.time(), "%m.%d.%Y")
write.csv(dat, sprintf("CMflux_wholebrain_metrics_%s.csv", date), row.names=FALSE)
write.csv(dat_pc, sprintf("CMflux_particcoeffs_%s.csv", date), row.names=FALSE)
write.csv(dat_wmd, sprintf("CMflux_withinmoddegs_%s.csv", date), row.names=FALSE)
write.csv(dat_le, sprintf("CMflux_localeffs_%s.csv", date), row.names=FALSE)
write.csv(dat_ndi, sprintf("CMflux_nodedissocindices_%s.csv", date), row.names=FALSE)
write.csv(dat_clustcoeff, sprintf("CMflux_netclustcoeff_%s.csv", date), row.names=FALSE)

# option 2: save these data.frames (and whatever else you want) as an R object
##load items saved in an R object (.RData files) into another R script with load("filename.RData")
date=format(Sys.time(), "%m.%d.%Y")
save(dat, dat_pc, dat_wmd, dat_le, dat_ndi, dat_clustcoeff, file = sprintf("CMflux_graphoutputs_weighted_absthresh010_%s.RData", date)) #update name with (i) binarized or weighted, (ii) the threshold type and level
#save(dat_clustcoeff, file = sprintf("CMflux_netclustcoeff_%s.RData", date))






#### C. construct network averages for WMD -- 2/25/2020, no longer doing this because you get 1e-17 for basically every module (because WMD is z-scored within each module) so there are no difference between tasks. will only use WMD to calculate provincial/connector hubs. ####
## NOTE: should write a function to calculate WMD without scaling (JRC + TRH suggestion)

# within_module_deg_z_score calculattion from brainGraph:
# function (g, memb) 
  # {
  #   stopifnot(is_igraph(g))
  #   N <- max(memb)
  #   A <- as_adj(g, sparse = FALSE, names = FALSE)
  #   z <- Ki <- rep(0, nrow(A))
  #   Ksi <- sigKsi <- rep(0, N)
  #   for (S in seq_len(N)) {
    #     x <- A[memb == S, ] %*% (memb == S)
    #     Ki[memb == S] <- x
    #     Ksi[S] <- mean(x)
    #     sigKsi[S] <- sd(x)
    #   }
  #   z <- (Ki - Ksi[memb])/sigKsi[memb]
  #   z <- ifelse(!is.finite(z), 0, z)
  #   return(z)
  # }
# <bytecode: 0x7fc9d37a83c8>
# <environment: namespace:brainGraph>







#### D. calculate average PC for each network ####

# D1. prep for reshape
dat_pc_byroi <- dat_pc
col_for_reshape_pc <- paste(colnames(dat_pc_byroi)[4:217], "_", roi_comms_5corenets, sep = "")
colnames(dat_pc_byroi)[4:217] <- paste(colnames(dat_pc_byroi)[4:217], "_", roi_comms_5corenets, sep = "") # adds network affiliation into the column name

# D2. reshape from wide to long keeping sub, ses, and task as columns
dat_pc_byroi
dat_pc_byroi_long <- reshape(dat_pc_byroi, 
                              direction = "long",
                              varying = list(names(dat_pc_byroi)[4:217]),
                              v.names = "pc",
                              idvar = c("sub", "ses", "task"),
                              timevar = "roi",
                              times = col_for_reshape_pc)

# D3. calculate averages for each network
pc_netavg <- ddply(dat_pc_byroi_long, .(sub, task), summarize, 
                    avg_pc_net1=mean(pc[grepl("_1", roi)], na.rm=TRUE), #network 1 = CON
                    avg_pc_net2=mean(pc[grepl("_2", roi)], na.rm=TRUE), #network 2 = DMN
                    avg_pc_net3=mean(pc[grepl("_3", roi)], na.rm=TRUE), #network 3 = FPN
                    avg_pc_net4=mean(pc[grepl("_4", roi)], na.rm=TRUE), #network 4 = SMD
                    avg_pc_net5=mean(pc[grepl("_5", roi)], na.rm=TRUE)) #network 5 = VIS

# D4. save
write.csv(pc_netavg, sprintf("CMflux_particcoeffsNETAVG_%s.csv", date), row.names=FALSE)
#save(pc_netavg, file = sprintf("MEMFYP_particcoeffsNEGTAVG_%s.RData", date))








#### E. calculate average local efficiency for each network ####

##NOTE FOR CLEANTHIS: calcualte local eff for each network and look for outliers

# E1. prep for reshape
dat_le_byroi <- dat_le
col_for_reshape_le <- paste(colnames(dat_le_byroi)[4:217], "_", roi_comms_5corenets, sep = "")
colnames(dat_le_byroi)[4:217] <- paste(colnames(dat_le_byroi)[4:217], "_", roi_comms_5corenets, sep = "") # adds network affiliation into the column name

# E2. reshape from wide to long keeping sub, ses, and task as columns
dat_le_byroi
dat_le_byroi_long <- reshape(dat_le_byroi, 
                             direction = "long",
                             varying = list(names(dat_le_byroi)[4:217]),
                             v.names = "le",
                             idvar = c("sub", "ses", "task"),
                             timevar = "roi",
                             times = col_for_reshape_le)

# E3. calculate averages for each network
le_netavg <- ddply(dat_le_byroi_long, .(sub, task), summarize, 
                   avg_le_net1=mean(le[grepl("_1", roi)], na.rm=TRUE), #network 1 = CON
                   avg_le_net2=mean(le[grepl("_2", roi)], na.rm=TRUE), #network 2 = DMN
                   avg_le_net3=mean(le[grepl("_3", roi)], na.rm=TRUE), #network 3 = FPN
                   avg_le_net4=mean(le[grepl("_4", roi)], na.rm=TRUE), #network 4 = SMD
                   avg_le_net5=mean(le[grepl("_5", roi)], na.rm=TRUE), #network 5 = VIS
                   overall_le = mean(le, na.rm = TRUE)) 


# E4. save
write.csv(le_netavg, sprintf("CMflux_localeffsNETAVG_%s.csv", date), row.names=FALSE)
#save(le_netavg, file = sprintf("MEMFYP_localeffsNETAVG_%s.RData", date))








#### F. calculate average node dissociation index (NDI) for each network ####

# F1. prep for reshape
dat_ndi_byroi <- dat_ndi
col_for_reshape_ndi <- paste(colnames(dat_ndi_byroi)[4:217], "_", roi_comms_5corenets, sep = "")
colnames(dat_ndi_byroi)[4:217] <- paste(colnames(dat_ndi_byroi)[4:217], "_", roi_comms_5corenets, sep = "") # adds network affiliation into the column name

# F2. reshape from wide to long keeping sub, ses, and task as columns
dat_ndi_byroi
dat_ndi_byroi_long <- reshape(dat_ndi_byroi, 
                             direction = "long",
                             varying = list(names(dat_ndi_byroi)[4:217]),
                             v.names = "ndi",
                             idvar = c("sub", "ses", "task"),
                             timevar = "roi",
                             times = col_for_reshape_ndi)

# F3. calculate averages for each network
ndi_netavg <- ddply(dat_ndi_byroi_long, .(sub, task), summarize, 
                   avg_ndi_net1=mean(ndi[grepl("_1", roi)], na.rm=TRUE), #network 1 = CON
                   avg_ndi_net2=mean(ndi[grepl("_2", roi)], na.rm=TRUE), #network 2 = DMN
                   avg_ndi_net3=mean(ndi[grepl("_3", roi)], na.rm=TRUE), #network 3 = FPN
                   avg_ndi_net4=mean(ndi[grepl("_4", roi)], na.rm=TRUE), #network 4 = SMD
                   avg_ndi_net5=mean(ndi[grepl("_5", roi)], na.rm=TRUE)) #network 5 = VIS

# F4. save
write.csv(ndi_netavg, sprintf("CMflux_ndisNETAVG_%s.csv", date), row.names=FALSE)
#save(ndi_netavg, file = sprintf("MEMFYP_ndisNETAVG_%s.RData", date))










#### NOTE: PROBABLY NEED TO QUANTIFY HOW MANY ROIS ARE USABLE IN EACH NETWORK FOR EACH TASK FOR EACH PERSON ####


