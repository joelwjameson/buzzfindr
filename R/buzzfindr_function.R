#BUZZ DETECTION FUNCTION

#Author: Joel Jameson
#Date: September 2022
#Makes use of the Bioacoustics R packages (Marchal et al. 2022).

#start of function
buzzfindr <- function(path, passes=FALSE, exp=1, buzzprob=0.8, clickfiltr=FALSE, channel="left", minfreq=20, maxfreq=45, threshold=4) {

#Load click filter
system.file("inst/clickfilter function.R", package="buzzfindr")
source("inst/clickfilter function.R")

#Load buzz detection model
system.file("inst/rf_model(shortfiles_SM2,SM4,MINI)10_7_22.rds", package="buzzfindr") #Model built from 3 detector types
model.buzz <- readRDS("inst/rf_model(shortfiles_SM2,SM4,MINI)10_7_22.rds") #Model built from 3 detector types

#load libraries
library(seewave)
library(bioacoustics)
library(tuneR)
library(runner)
library(svMisc)
library(rlist)
library(dplyr)
library(randomForest)
require(caTools)

#working directory
setwd <- path

#Create a txt file that will take the results
d <- data.frame(t(c("Time", "File", "Probability")))
write.table(d,paste0(path,"/Buzz_Results.txt"),row.names = FALSE, col.names = FALSE)

#load data
files <- list.files(path)
files <- files[grep("wav", files)]

#subset just bat passes
if (isTRUE(passes==FALSE)){
  files <- files
} else
  if (isTRUE(passes[1]!=FALSE) && length(passes)>0){
    files <- passes
  }

#create a dataframe for files
bat_data <- as.data.frame(files)
names(bat_data)[1] <- "File"
bat_data$Path <- paste0(path,"/",files)
#bat_data <- bat_data[c(1:100),]
Buzz_data <- NULL

#set the time expansion factor
exp <- exp

#Set buzz probability
buzzprob <- buzzprob

#Set other parameters
minfreq=minfreq*1000
maxfreq=maxfreq*1000
threshold=threshold

#Load libraries for the click filter
library(seewave)
library(bioacoustics)
library(tuneR)
library(runner)
library(svMisc)
library(rlist)
library(dplyr)

#Set the click filter threshold
if(isTRUE(clickfiltr)==TRUE){
  clickpercent <- 0.5} else if (is.numeric(clickfiltr)==TRUE){
    clickpercent <- clickfiltr} else {
      clickpercent <- NULL
    }


#Object that will take potentially corrupt files
cor <- NULL

#Specify the time window size (number of consecutive detections to analyze in the recording)
sk <- 4

for (i in 1:nrow(bat_data)){

  #______________________________________________
  #Progress bar
  svMisc::progress(i, max.value = nrow(bat_data))
  Sys.sleep(0.01)
  if (i==nrow(bat_data)) {cat("Done!  !  !  ! ! !!!!!!\n")}
  #_______________________________________________

  #read wave file (control for corrupt files)
  if(file.size(bat_data$Path[i])>0){w <- readWave(bat_data$Path[i])} else{
    corname <- bat_data$File[i]
    cor <- rbind(cor, corname)
    next
  }

  #List of frequencies delimiting bandwidths to scan
  LPF.list <- seq(minfreq+5000, maxfreq, by=5000)
  HPF.list <- seq(minfreq, maxfreq-5000, by=5000)

  #List that will take parameter extraction results
  w.list <- list()

  #Create the dataframe to take the parameters
  paramsall <- NULL

  #Loop through the frequency windows
  for (z in 1:length(LPF.list)){

    capture.output(ts <- threshold_detection(w,
                                             channel=channel,
                                             min_dur=0.5,
                                             min_TBE=2,
                                             threshold=threshold,
                                             LPF=LPF.list[z],
                                             HPF=HPF.list[z],
                                             time_exp = exp,
                                             FFT_size = 256,
                                             FFT_overlap = 0.85,
                                             SNR_thr=6,
                                             duration_thr=80,
                                             angle_thr=40,
                                             start_thr=40,
                                             end_thr=20,
                                             NWS=50),
                   file='NUL', type="message")
    w.list[z] <- assign(paste0("t",z), ts)


    if(length(ts$data$event_data$starting_time)>1){
      tstart <- ts$data$event_data$starting_time
      tstart <- substr(tstart, start = 7, stop = 12)
      tstart <- as.character(tstart)
      tstart <- as.numeric(tstart)
      #get the differences
      time <- as.data.frame(tstart)
      time$duration <- ts$data$event_data$duration/1000
      time$slope <- ts$data$event_data$slope
      tstart2 <- time$tstart[2:length(time$tstart)]
      time <- time[1:(nrow(time)-1),]
      time$start2 <- tstart2
      time$diff <- time$start2-time$tstart

      #Extract the parameters from Time (frequency window) using a moving window
      params <- as.data.frame(t(runner::runner(
        x = time,
        k = sk,
        f = function(time) {
          c(mean(time$diff),
            sd(time$slope),
            mean(time$duration),
            sd(time$duration))
        }
      )))

      params <- cbind(time$tstart, params)
      names(params) <- c("tstart","avgt","sdslope","avgdur","sddur")
      params <- params[-c(1:4),]#remove first lines. Often gives rsquared=1

      #Combine params from each frequency bandwidth
      paramsall <- rbind(paramsall, params)

    } else {next}
  }

  #Rename the combined parameter data as params
  params <- paramsall

  #remove sd so no issues with sd function
  rm(sd)

  if(isTRUE(nrow(params)>0)==TRUE){

  #Add file name to params
  params$file <- bat_data$File[i]

  #Make prediction based on model
  #Combine
  test.data <- params[c(complete.cases(params)),]

  # Make predictions on test data
  predictions <- model.buzz %>% predict(test.data)

  #Add the prediction to the dataframe as a variable
  test.data$result <- predictions

  #Change the prediction to a grouping variable. Here we use a Discriminant Probability of 95%
  test.data$group <- NA
  test.data$group[test.data$result>=buzzprob] <- "buzz"
  test.data$group[test.data$result<buzzprob] <- NA

  #subset just buzzes
  test.data <- test.data[c(complete.cases(test.data)),]

  #Go to next file if test.data is empty
  if (nrow(test.data)==0){
    next
  } else {
    test.data <- test.data
  }
  } else {
    next
  }

  #Have to remove object sd before going further because need to re-use aggregate and call on sd function
  #but the sd function from aggregate won't work if sd is also an object
  rm(sd)


  #-------------------------------------------------------------------------------
  #RUN THROUGH THE CLICK FILTER
  if(is.numeric(clickpercent)==TRUE & nrow(test.data>0)){

    #Format the time data in w.list before pasing to the clickfilter (ie.make start time numeric)
    for (s in 1:length(LPF.list)){
      w.list[[s]]$event_data$starting_time <- substr(w.list[[s]]$event_data$starting_time, start = 7, stop = 12)
      w.list[[s]]$event_data$starting_time <- as.character(w.list[[s]]$event_data$starting_time)
      w.list[[s]]$event_data$starting_time <- as.numeric(w.list[[s]]$event_data$starting_time)
    }

    #Pass to clickfilter
    w.list <- filtr_clicks(w=w, exp=exp, w.list=w.list, LPF.list=LPF.list, clickpercent=clickpercent, channel=channel)

    #Create the timecheck vector of time intervals to filter out obvious non-buzzes in the next step
    timecheck <- NULL

    #Create the list of dataframes from which trend parameters will be calculated
    bands <- list()

    #dataframe for parameters
    paramsall <- NULL

    for (s in 1:length(LPF.list)){

      #Extract buzz-relevant parameters using a sliding time window along the file
      ts <- w.list[[s]]$event_data

      if(length(ts$starting_time)>1){
        #tstart
        tstart <- ts$starting_time
        #get the differences
        time <- as.data.frame(tstart)
        time$duration <- ts$duration/1000
        time$slope <- ts$slope
        tstart2 <- time$tstart[2:length(time$tstart)]
        time <- time[1:(nrow(time)-1),]
        time$start2 <- tstart2
        time$diff <- time$start2-time$tstart

        #Extract the parameters from Time (frequency window) using a moving window
        params <- as.data.frame(t(runner::runner(
          x = time,
          k = sk,
          f = function(time) {
            c(mean(time$diff),
              sd(time$slope),
              mean(time$duration),
              sd(time$duration))
          }
        )))

        params <- cbind(time$tstart, params)
        names(params) <- c("tstart","avgt","sdslope","avgdur","sddur")
        params <- params[-c(1:4),]#remove first lines. Often gives rsquared=1

        #Combine params from each frequency bandwidth
        paramsall <- rbind(paramsall, params)

      } else {next}
    }

    #Rename the combined parameter data as params
    params <- paramsall

    if(isTRUE(nrow(params)>0)==TRUE){
      #Add file name to params
      params$file <- bat_data$File[i]

      #Make prediction based on model
      #Combine
      test.data <- params[c(complete.cases(params)),]

      # Make predictions on test data
      predictions <- model.buzz %>% predict(test.data)

      #Add the prediction to the dataframe as a variable
      test.data$result <- predictions########################RF

      #Change the prediction to a grouping variable. Here we use a Discriminant Probability of 95%
      test.data$group <- NA
      test.data$group[test.data$result>=buzzprob] <- "buzz"
      test.data$group[test.data$result<buzzprob] <- NA

      #subset just buzzes
      test.data <- test.data[c(complete.cases(test.data)),]

      #After filtering, if test.data is now empty, go to next file
      if (nrow(test.data)==0){
        next} else {
          test.data <- test.data
        }
    } else {
      next
    }
  }


  #Identify the separate buzzes (how many are there?)
  test.data <- test.data[c(order(test.data$tstart)),]
  test <- test.data$tstart
  if(length(test)>1){
    test.data$diff[2:length(test.data$tstart)] <- test[2:length(test)] - test[1:(length(test)-1)]
  } else {
    test.data <- test.data}
  test.data$diff[1] <- 1
  test.data$buzz <- NA
  test.data$buzz[c(test.data$diff>0.2)] <- "buzz"

  #Get the unique file names
  test.data <- test.data[c(complete.cases(test.data)),]
  test.data <- test.data[,c(1,6,7)]

  #Export Results (will constantly update so nothing is lost if the function stops suddenly)
  write.table (test.data, paste0(path,"/Buzz_Results.txt"), append=TRUE, row.names = FALSE, col.names = FALSE)

  #Export Potentially corrupted files
  write.table (cor, paste0(path,"/Files_Not_Analyzed (May be Corrupted).txt"), append=TRUE, row.names = FALSE, col.names = FALSE)

  #Add to the results dataframe
  Buzz_data <- rbind(Buzz_data, test.data)

}

Buzz_data

}
