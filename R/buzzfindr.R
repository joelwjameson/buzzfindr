#BUZZ DETECTION FUNCTION

#Author: Joel Jameson
#Date: September 2022
#Makes use of the Bioacoustics R packages (Marchal et al. 2022).

#start of function
buzzfindr <- function(path, 
                      passes=FALSE, 
                      pass_thresh=FALSE,
                      exp=1, 
                      buzzprob=0.8, 
                      channel="left", 
                      minfreq=20, 
                      maxfreq=50, 
                      out.file=NULL,
                      num_cores=FALSE,
                      vote=FALSE) {
  
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
  library(doSNOW)
  library(utils)
  library(parallel)
  library(car)
  library(readr)
  library(entropy)
  
  #Load buzz detection model
  z1 <- system.file("rf_TEST18(17_3_24).rds", package="buzzfindr")
  model.buzz <- readRDS(z1)
  
  #Scale for variables
  z2 <- system.file("scaler.rds", package="buzzfindr")
  scaler <- readRDS(z2)
  
  #Detection parameters
  z3 <- system.file("det_parameters.rds", package="buzzfindr")
  xallsub <- readRDS(z3)
  
  #Model name
  mod.name <- "rf_TEST18(17_3_24)"
  
  #working directory
  setwd <- path
  
  #load data
  files <- list.files(path, pattern=".wav")
  
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
  
  #set the time expansion factor
  exp <- exp
  
  #Set buzz probability
  buzzprob <- buzzprob
  
  #Set other parameters
  minfreq<-minfreq*1000
  maxfreq<-maxfreq*1000
  
  #Set the threshold for pass filter
  if(pass_thresh>0 & isTRUE(pass_thresh)==TRUE){
    pass_thresh <- 8
  } else if(pass_thresh>0 & isTRUE(pass_thresh)==FALSE){
    pass_thresh <- pass_thresh} else{
      pass_thresh <- NULL
    }
  
  #Object that will take potentially corrupt files
  cor <- NULL
  
  #Specify the time window size (number of consecutive detections to analyze in the recording)
  sk <- 4
  
  #channel
  channel <- channel
  
  #num_cores
  num_cores <- num_cores

  if(isTRUE(num_cores==FALSE)==TRUE){
  ###########################BASIC LOOP#########################################
  #Run through files
  Buzz_data <- NULL
  tryCatch({#Ignore errors
  for(i in 1:length(bat_data$Path)){

                          #Progress bar
                          svMisc::progress(i, max.value = nrow(bat_data))
                          Sys.sleep(0.01)
                          if (i==nrow(bat_data)) {cat("Done!  !  !  ! ! !!!!!!\n")}
    
                             #read wave file (control for corrupt files)
                           if(file.size(bat_data$Path[i])>0){w <- readWave(bat_data$Path[i])} else{
                             corname <- bat_data$File[i]
                             cor <- rbind(cor, corname)
                             next
                           }
                           
                           #Initial filter for bat calls
                           if(isTRUE(pass_thresh>0)==TRUE){
                             tcheck <- threshold_detection(w, threshold = pass_thresh, time_exp=exp, NWS=50, SNR_thr=6)
                             if(length(tcheck$data$event_data$starting_time)==0){
                               next}
                           }
                           
                           #List of frequencies delimiting bandwidths to scan
                           LPF.list <- seq(minfreq+5000, maxfreq, by=5000)
                           HPF.list <- seq(minfreq, maxfreq-5000, by=5000)
                           
                           #Create the dataframe to take the parameters
                           paramsall <- NULL
                           
                           #Loop through the frequency bands
                           for (z in 1:length(LPF.list)){
                             
                             capture.output(ts <- threshold_detection(w,
                                                                      min_dur=xallsub$min_dur[1],
                                                                      #max_dur = 4,
                                                                      min_TBE=xallsub$min_TBE[1],
                                                                      max_TBE=20,
                                                                      threshold=xallsub$threshold[1],
                                                                      LPF=LPF.list[z],
                                                                      HPF=HPF.list[z],
                                                                      time_exp = exp,
                                                                      FFT_size = 256,
                                                                      FFT_overlap = xallsub$FFT_overlap[1],
                                                                      SNR_thr=xallsub$SNR_thr[1],
                                                                      duration_thr=xallsub$duration_thr[1],
                                                                      angle_thr=xallsub$angle_thr[1],
                                                                      start_thr=xallsub$start_thr[1],
                                                                      end_thr=xallsub$end_thr[1],
                                                                      NWS=xallsub$NWS[1]),
                                            file='NUL', type="message")
                             
                             if(length(ts$data$event_data$starting_time)>1){
                               tstart <- ts$data$event_data$starting_time
                               tstart <- substr(tstart, start = 7, stop = 12)
                               tstart <- as.character(tstart)
                               tstart <- as.numeric(tstart)
                               #get the differences
                               time <- as.data.frame(tstart)
                               time$snr <- ts$data$event_data$snr
                               time$duration <- ts$data$event_data$duration/1000
                               tstart2 <- time$tstart[2:length(time$tstart)]
                               time <- time[1:(nrow(time)-1),]
                               time$start2 <- tstart2
                               time$diff <- time$start2-time$tstart
                               time$File <- bat_data$Path[i]
                               
                               #Extract the parameters from Time (frequency window) using a moving window
                               params <- as.data.frame(t(runner::runner(
                                 x = time,
                                 k = sk,
                                 f = function(time) {
                                   c(max(time$diff),
                                     sd(time$diff),
                                     mean(time$duration),
                                     mean(time$snr),
                                     min(time$diff)
                                   )
                                 }
                               )))
                               
                               params <- cbind(time$tstart, params, HPF.list[z])
                               names(params) <- c("tstart", "maxt","sd", "avgdur", "snravg","mint","HPF")
                               params <- params[-c(1:2),]
                               
                             } else {params <- NULL}
                             
                             #Combine params from each frequency bandwidth
                             paramsall <- rbind(paramsall, params)
                           }
                           
                           #Rename the combined parameter data as params
                           params <- paramsall
                           
                           #Before predicting, additional filtering for improbable buzz parameters
                           params <- params[c(params$maxt<0.014),]
                           
                           #Scale
                           if(isTRUE(nrow(params)>0)==TRUE){
                             
                             #Scale
                             scaler <- scaler
                             params[,c(2:6)] <- scale(params[,c(2:6)], center = attr(scaler, "scaled:center")[c(4,6,14,17,3)], scale = attr(scaler, "scaled:scale")[c(4,6,14,17,3)])  
                             
                             #Add file name to params
                             params$file <- bat_data$File[i]
                             
                             #Combine
                             test.data <- params[c(complete.cases(params)),]
                             
                             # Make predictions on test data
                             test.data$result <- model.buzz %>% predict(test.data)#RF model
                             #test.data$result <- neuralnet::compute(model.buzz,test.data)$net.result#nn Model
                             
                             #Change the prediction to a grouping variable.
                             test.data$group <- NA
                             test.data$group[test.data$result>=buzzprob] <- "Buzz"
                             
                             #Go to next file if test.data is empty
                             if (length(na.omit(test.data$group))==0){
                               next
                             } else {
                               test.data <- test.data
                             }
                           } else {
                             next
                           }
                           
                           #-----------------------------------------------------------------------------
                           # Apply Vote procedure if called for
                           if(vote==TRUE){
                             if(nrow(test.data)==1){next}
                             test.data <- test.data[c(order(test.data$HPF, test.data$tstart)),]
                             test.datasub <- NULL
                             for(l in 1:(nrow(test.data)-1)){
                               if(is.na(test.data$group[l])==FALSE & is.na(test.data$group[l+1])==FALSE){
                                 x <- test.data[l,]
                               } else {next}
                               test.datasub <- rbind(test.datasub, x)
                             }
                             test.data <- test.datasub
                             if(is.null(test.data)==TRUE){next}
                           }
                           
                           #subset just buzzes
                           test.data <- test.data[c(complete.cases(test.data)),]
                           
                           #-----------------------------------------------------------------------------
                           #Identify the separate buzzes (how many are there?)
                           #Get time difference between buzz detections
                           test.data <- test.data[c(order(test.data$tstart)),]
                           test <- test.data$tstart
                           if(length(test)>1){
                             test.data$diff[2:length(test.data$tstart)] <- test[2:length(test)] - test[1:(length(test)-1)]
                           } else {
                             test.data <- test.data}
                           
                           #Apply threshold for time between buzz detections
                           test.data$diff[1] <- 1
                           test.data$buzz <- NA
                           test.data$buzz[c(test.data$diff>0.2)] <- "Buzz"
                           
                           #Get the unique file names
                           test.data <- test.data[c(complete.cases(test.data$buzz)),]
                           test.data <- test.data[,c("tstart","file","result","buzz")]
                           test.data <- test.data
                         
                           #Combine
                           Buzz_data <- rbind(Buzz_data, test.data)
                        }
  }, error=function(e){})#End of tryCatch
  
  #Dataframe for results
  if(isTRUE(nrow(Buzz_data)==0) | is.null(Buzz_data)){
    Buzz_data <- data.frame(tstart=NA,file=NA,result=NA,buzz=NA)
  }
  
  #Create results directory
  if(isTRUE(length(Buzz_data)>0)==TRUE){
    dirx <- paste0("Buzz_Results_",Sys.Date(),"_",format(Sys.time(),"%H-%M-%S"),"")
    dir.create(paste0(path,"/", dirx))
  }
  
  #--------------------------------------
  #Outputs
  if(isTRUE("png" %in% out.file)==TRUE){
    tryCatch({#Ignore errors
      #Output Results
      for(i in 1:nrow(Buzz_data)){
        file <- as.character(paste0(path,"/",Buzz_data$file[i]))
        start <- Buzz_data$tstart[i]-0.15
        end <- Buzz_data$tstart[i]+0.25
        label <- gsub(".wav",paste0("-",Buzz_data$tstart[i]), Buzz_data$file[i])
        sound <- readWave(file, from = start, to = end, units = "seconds")
        png(paste0(path,"/",dirx,"/",label,".png"),width = 500, height = 400)
        spectro(sound, FFT_size = 1024,col=gray.colors(50, 1, 0))
        dev.off()
      }
    }, error=function(e){})#End of tryCatch
  }
  
  if(isTRUE("wav" %in% out.file)==TRUE){
    tryCatch({#Ignore errors
      #Output Results
      for(i in 1:nrow(Buzz_data)){
        file <- as.character(paste0(path,"/",Buzz_data$file[i]))
        start <- Buzz_data$tstart[i]-0.15
        end <- Buzz_data$tstart[i]+0.25
        label <- gsub(".wav",paste0("-",Buzz_data$tstart[i]), Buzz_data$file[i])
        sound <- readWave(file, from = start, to = end, units = "seconds")
        writeWave(sound, paste0(path,"/",dirx,"/",label,".wav"))
      }
    }, error=function(e){})#End of tryCatch
  }
  
  #Add process info
  Buzz_data$path <- path
  Buzz_data$passes <- ifelse(is.null(passes) ,"FALSE", as.character(passes))
  Buzz_data$pass_thresh <- ifelse(is.null(pass_thresh) ,"FALSE", as.character(pass_thresh))
  Buzz_data$exp <- exp
  Buzz_data$buzzprob <- buzzprob
  Buzz_data$channel <- channel
  Buzz_data$minfreq <- minfreq
  Buzz_data$maxfreq <- maxfreq
  Buzz_data$model <- mod.name
  
  if(isTRUE("csv" %in% out.file)==TRUE){
    write.csv(Buzz_data, paste0(path,"/",dirx,"/Buzz_Results.csv"))
  }
  
  out <- Buzz_data
  
} 
 ############################################################################### 

  
  
  
  
  
  
  
  else {
  #########################PARALLEL PROCESSING##################################
  #Run through files
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)
  
  #Progress Bar
  Start.time <- Sys.time()
  pb <- txtProgressBar(max = length(bat_data$Path), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  #tryCatch({#Ignore errors
    Buzz_data <- foreach(bat_data.path=bat_data$Path,
                         bat_data.file = bat_data$File,
                         .combine = 'rbind',
                         .options.snow = opts,
                         .packages=c("plyr", "dplyr", "tidyr", "car", "glmmTMB", "readr","bioacoustics","seewave","tuneR","runner","randomForest","caTools", "entropy")) %dopar% {
                          xallsub=xallsub
                          exp=exp
                          buzzprob=buzzprob
                          minfreq=minfreq
                          maxfreq=maxfreq
                          cor=cor
                          files=files
                          model.buzz=model.buzz
                          sk=sk
                          channel=channel
                           
                           #read wave file (control for corrupt files)
                           if(file.size(bat_data.path)>0){w <- readWave(bat_data.path)} else{
                             corname <- bat_data.file
                             cor <- rbind(cor, corname)
                             return(NULL)
                           }
                           
                           #Initial filter for bat calls
                           if(isTRUE(pass_thresh>0)==TRUE){
                             tcheck <- threshold_detection(w, threshold = pass_thresh, time_exp=exp, NWS=50, SNR_thr=6)
                             if(length(tcheck$data$event_data$starting_time)==0){
                               return(NULL)}
                           }
                           
                           #List of frequencies delimiting bandwidths to scan
                           LPF.list <- seq(minfreq+5000, maxfreq, by=5000)
                           HPF.list <- seq(minfreq, maxfreq-5000, by=5000)
                           
                           #Create the dataframe to take the parameters
                           paramsall <- NULL
                           
                           #Loop through the frequency bands
                           for (z in 1:length(LPF.list)){
                             
                             capture.output(ts <- threshold_detection(w,
                                                                      min_dur=xallsub$min_dur[1],
                                                                      #max_dur = 4,
                                                                      min_TBE=xallsub$min_TBE[1],
                                                                      max_TBE=20,
                                                                      threshold=xallsub$threshold[1],
                                                                      LPF=LPF.list[z],
                                                                      HPF=HPF.list[z],
                                                                      time_exp = exp,
                                                                      FFT_size = 256,
                                                                      FFT_overlap = xallsub$FFT_overlap[1],
                                                                      SNR_thr=xallsub$SNR_thr[1],
                                                                      duration_thr=xallsub$duration_thr[1],
                                                                      angle_thr=xallsub$angle_thr[1],
                                                                      start_thr=xallsub$start_thr[1],
                                                                      end_thr=xallsub$end_thr[1],
                                                                      NWS=xallsub$NWS[1]),
                                            file='NUL', type="message")
 
                             if(length(ts$data$event_data$starting_time)>1){
                               tstart <- ts$data$event_data$starting_time
                               tstart <- substr(tstart, start = 7, stop = 12)
                               tstart <- as.character(tstart)
                               tstart <- as.numeric(tstart)
                               #get the differences
                               time <- as.data.frame(tstart)
                               time$snr <- ts$data$event_data$snr
                               time$duration <- ts$data$event_data$duration/1000
                               tstart2 <- time$tstart[2:length(time$tstart)]
                               time <- time[1:(nrow(time)-1),]
                               time$start2 <- tstart2
                               time$diff <- time$start2-time$tstart
                               time$File <- bat_data.path
                               
                               #Extract the parameters from Time (frequency window) using a moving window
                               params <- as.data.frame(t(runner::runner(
                                 x = time,
                                 k = sk,
                                 f = function(time) {
                                   c(max(time$diff),
                                     sd(time$diff),
                                     mean(time$duration),
                                     mean(time$snr),
                                     min(time$diff)
                                   )
                                 }
                               )))
                               
                               params <- cbind(time$tstart, params, HPF.list[z])
                               names(params) <- c("tstart", "maxt","sd", "avgdur", "snravg","mint","HPF")
                               params <- params[-c(1:2),]
                               
                             } else {params <- NULL}
                             
                             #Combine params from each frequency bandwidth
                             paramsall <- rbind(paramsall, params)
                           }
                           
                           #Rename the combined parameter data as params
                           params <- paramsall
                           
                           #Before predicting, additional filtering for improbable buzz parameters
                           params <- params[c(params$maxt<0.014),]
                           
                           #Scale
                           if(isTRUE(nrow(params)>0)==TRUE){
                             
                             #Scale
                             scaler <- scaler
                             params[,c(2:6)] <- scale(params[,c(2:6)], center = attr(scaler, "scaled:center")[c(4,6,14,17,3)], scale = attr(scaler, "scaled:scale")[c(4,6,14,17,3)])  
                           
                             #Add file name to params
                             params$file <- bat_data.file
                             
                             #Combine
                             test.data <- params[c(complete.cases(params)),]
                             
                             # Make predictions on test data
                             test.data$result <- model.buzz %>% predict(test.data)#RF model
                             #test.data$result <- neuralnet::compute(model.buzz,test.data)$net.result#nn Model
                             
                             #Change the prediction to a grouping variable.
                             test.data$group <- NA
                             test.data$group[test.data$result>=buzzprob] <- "Buzz"

                             #Go to next file if test.data is empty
                             if (length(na.omit(test.data$group))==0){
                               return(NULL)
                             } else {
                               test.data <- test.data
                             }
                           } else {
                             return(NULL)
                           }
                           
#-----------------------------------------------------------------------------
                           # Apply Vote procedure if called for
                           if(vote==TRUE){
                             if(nrow(test.data)==1){return(NULL)}
                             test.data <- test.data[c(order(test.data$HPF, test.data$tstart)),]
                             test.datasub <- NULL
                             for(l in 1:(nrow(test.data)-1)){
                               if(is.na(test.data$group[l])==FALSE & is.na(test.data$group[l+1])==FALSE){
                                 x <- test.data[l,]
                               } else {next}
                               test.datasub <- rbind(test.datasub, x)
                             }
                             test.data <- test.datasub
                             if(is.null(test.data)==TRUE){return(NULL)}
                           }
                           
                           #subset just buzzes
                           test.data <- test.data[c(complete.cases(test.data)),]
                           
#-----------------------------------------------------------------------------
                           #Identify the separate buzzes (how many are there?)
                           #Get time difference between buzz detections
                           test.data <- test.data[c(order(test.data$tstart)),]
                           test <- test.data$tstart
                           if(length(test)>1){
                             test.data$diff[2:length(test.data$tstart)] <- test[2:length(test)] - test[1:(length(test)-1)]
                           } else {
                             test.data <- test.data}
                           
                           #Apply threshold for time between buzz detections
                           test.data$diff[1] <- 1
                           test.data$buzz <- NA
                           test.data$buzz[c(test.data$diff>0.2)] <- "Buzz"
                           
                           #Get the unique file names
                           test.data <- test.data[c(complete.cases(test.data$buzz)),]
                           test.data <- test.data[,c("tstart","file","result","buzz")]
                           test.data <- test.data
                         }
  #}, error=function(e){})#End of tryCatch
  
  #Stop cluster
  stopCluster(cl)

  #Dataframe for results
  if(isTRUE(nrow(Buzz_data)==0) | is.null(Buzz_data)){
    Buzz_data <- data.frame(tstart=NA,file=NA,result=NA,buzz=NA)
  }
  
  #Create results directory
  if(isTRUE(length(Buzz_data)>0)==TRUE){
    dirx <- paste0("Buzz_Results_",Sys.Date(),"_",format(Sys.time(),"%H-%M-%S"),"")
    dir.create(paste0(path,"/", dirx))
  }
  
  #--------------------------------------
  #Outputs
  if(isTRUE("png" %in% out.file)==TRUE){
    tryCatch({#Ignore errors
      #Output Results
      for(i in 1:nrow(Buzz_data)){
        file <- as.character(paste0(path,"/",Buzz_data$file[i]))
        start <- Buzz_data$tstart[i]-0.15
        end <- Buzz_data$tstart[i]+0.25
        label <- gsub(".wav",paste0("-",Buzz_data$tstart[i]), Buzz_data$file[i])
        sound <- readWave(file, from = start, to = end, units = "seconds")
        png(paste0(path,"/",dirx,"/",label,".png"),width = 500, height = 400)
        spectro(sound, FFT_size = 1024,col=gray.colors(50, 1, 0))
        dev.off()
      }
    }, error=function(e){})#End of tryCatch
  }
  
  if(isTRUE("wav" %in% out.file)==TRUE){
    tryCatch({#Ignore errors
      #Output Results
      for(i in 1:nrow(Buzz_data)){
        file <- as.character(paste0(path,"/",Buzz_data$file[i]))
        start <- Buzz_data$tstart[i]-0.15
        end <- Buzz_data$tstart[i]+0.25
        label <- gsub(".wav",paste0("-",Buzz_data$tstart[i]), Buzz_data$file[i])
        sound <- readWave(file, from = start, to = end, units = "seconds")
        writeWave(sound, paste0(path,"/",dirx,"/",label,".wav"))
      }
    }, error=function(e){})#End of tryCatch
  }
  
  #Add process info
  Buzz_data$path <- path
  Buzz_data$passes <- ifelse(is.null(passes) ,"FALSE", as.character(passes))
  Buzz_data$pass_thresh <- ifelse(is.null(pass_thresh) ,"FALSE", as.character(pass_thresh))
  Buzz_data$exp <- exp
  Buzz_data$buzzprob <- buzzprob
  Buzz_data$channel <- channel
  Buzz_data$minfreq <- minfreq
  Buzz_data$maxfreq <- maxfreq
  Buzz_data$model <- mod.name
  
  if(isTRUE("csv" %in% out.file)==TRUE){
    write.csv(Buzz_data, paste0(path,"/",dirx,"/Buzz_Results.csv"))
  }
  
  out <- Buzz_data
  
  }
}
################################################################################
