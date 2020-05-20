#### Data Management and Network Construction for KKI Task Analyses for Flux Conference
### Adapted from Mackenzie Mitchell's script for ADHD Med Challenge Analyses for Flux Conference and First-Year Project

rm(list=ls())

library(pracma)

files <- list.files("./bigbrain", full.names = T)
files <- files[grepl("bigbrain.csv", files, fixed = T)]
file_info <- do.call( "rbind",strsplit(basename(files), "_", fixed = T))[,1:4]
file_info = as.data.frame(file_info)
file_info$file = files
#QC = read.csv("QC_data_unblind.csv", stringsAsFactors = F)
colnames(file_info) = c("sub", "ses", "task", "run")
#file_info = merge(file_info, QC, by = c("sub", "ses", "task", "run"))

file_info$reject = FALSE
names(file_info)[5] = "roifiles"

file_info$reject[which((file_info[,1] == "sub-0941"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-0958"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1012"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1069"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1166"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1207"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1271"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1303"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1307"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1325"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1340"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1386"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1395"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1413"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1425"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1443"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1453"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1459"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1461"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1475"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1481"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1490"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1504"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1507"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1517"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1526"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1529"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1542"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1543"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1553"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1561"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1595"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1611"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1633"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1660"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1666"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1722"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1769"))] = TRUE
file_info$reject[which((file_info[,1] == "sub-1853"))] = TRUE


sub_list = unique(file_info[,1])
sub_list = levels(sub_list)[as.numeric(sub_list)]
sub_datalist <- list()


require(progress)
pb = progress_bar$new(total = length(sub_list), format = "[:bar] :current/:total :eta")


for(i in 1:length(sub_list)){
  #print(sub_list[i])
  sub_info <- file_info[which(file_info[,1] == sub_list[i]),]
  if(length(dim(unique(sub_info))) >1){
    ses_list <- unique(sub_info[, 2])
  }else{
    ses_list <- unique(sub_info[2])
  }
  ses_list = levels(ses_list)[as.numeric(ses_list)]
  sub_datalist[[sub_list[i]]] <- list()
  for(j in 1:length(ses_list)){
    sub_ses_info = file_info[which((file_info[,1] == sub_list[i]) & (file_info[,2] == ses_list[j])),]
    if(length(dim(unique(sub_ses_info))) >1){
      task_list <- unique(sub_ses_info[, 3])
    }else{
      task_list <- unique(sub_ses_info[3])
    }
    task_list = levels(task_list)[as.numeric(task_list)]
    sub_datalist[[sub_list[i]]][[ses_list[[j]]]] <- list()
    for(l in 1:length(task_list)){
      file_inds = which((file_info[,1] == sub_list[i]) & (file_info[,2] == ses_list[j]) & (file_info[,3] == task_list[l]) &  !(file_info$reject))
      goodtps = 0 # goodtps is specified per task
      goodtps_scan = 0
      #print(files[file_inds])
      
      
      ## For different scrubbing percentages across rest/GNG tasks, use the following code
      #if(task_list[l] == 'task-rest'){
      #  filter = 0.5  # and then goodtps_scan >= filter
      #}
      #if(task_list[l] == 'task-gng'){
      #  filter = 0.5  # and then goodtps_scan >= filter
      #}
      
      
      if(length(file_inds) >= 1){ ## asks how many runs in this subject/session/task
        data_list <- list()
        for(k in 1:length(file_inds)){ ## looping through each run of the task in the loop above
          data <- read.csv(file_info$roifiles[file_inds[k]], stringsAsFactors = F, header = F)
          goodtps_scan = sum(!is.na(data[,178]))  # changed from 1 to 178 to go to frontal ROI 178
          toRemove = which(is.na(data[,178]))  # changed from 1 to 178 to go to frontal ROI 178
          if(goodtps_scan >= 50){ #dont want to remove runs that have less than 78; want to take all tps from each run and add them together; then cutoff the total tps from the concatenated task runs t >=78
            if(length(toRemove) > 0){
              data_list[[k]] = scale(data[-toRemove,])
            }else{
              data_list[[k]] = scale(data)
            }
            goodtps = goodtps + goodtps_scan # adding the goodtps for each run to the goodtps for the task
          }
        }
        if(goodtps >= 78){  # i think use number 78? ##NOTE: probably move out one level; so trimming on the task level instead of the run level
          if(length(data_list) >1){
            sub_datalist[[sub_list[i]]][[ses_list[[j]]]][[task_list[l]]] <- do.call("rbind", data_list)
          }else{
            sub_datalist[[sub_list[i]]][[ses_list[[j]]]][[task_list[l]]] <- data_list[[1]]
          }
        }
      }
    }
    pb$tick()
  }
}
ROI_TS_list = unlist(unlist(sub_datalist, recursive = F), recursive = F)


sub_datalist_50tpsperrun = sub_datalist
ROI_TS_list_50tpsperrun = ROI_TS_list
save(sub_datalist_50tpsperrun, ROI_TS_list_50tpsperrun, file = "ConcatenatedROItimeseries_50tpsperrun.RData")


# NETWORK CONSTRUCTION ## loops through the subjects > sessions > tasks inside sub_datalist and runs correlations and outputs csv files for each subjects > sessions > tasks ## loop through those files for the network stats files (string split operation)

load(file = "ConcatenatedROItimeseries_50tpsperrun.RData")
sub_datalist = sub_datalist_50tpsperrun

sub_netlist <- list() # lists of lists: subject > session ##NEED TO ADD A TASK LAYER HERE
for(i in 1:length(sub_datalist)){
  sub_netlist[[names(sub_datalist)[i]]] <- list()
  for(j in 1:length(sub_datalist[[i]])){
    sub_netlist[[names(sub_datalist)[i]]][[names(sub_datalist[[i]])[j]]] <- list()
    for(l in 1:length(sub_datalist[[i]][[j]])){
      if(length(sub_datalist[[i]][[j]]) >= 1){
        sub_netlist[[names(sub_datalist)[i]]][[names(sub_datalist[[i]])[j]]][[names(sub_datalist[[i]][[j]])[l]]] <- list()
        tryCatch({
          cor_mat = cor(sub_datalist[[i]][[j]][[l]], use = "pairwise.complete.obs")
          sub_netlist[[i]][[j]][[l]] <- cor_mat
          file_name = paste("./correlation_matrices50tps/", names(sub_datalist)[i], "_", names(sub_datalist[[i]])[j], "_", names(sub_datalist[[i]][[j]])[l], ".csv", sep = "")
          write.csv(cor_mat, file_name, row.names = FALSE)
        },error = function(e){print("Issue")} )
      }
      
    }
  }
  #pb$tick()
}










