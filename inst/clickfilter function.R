########################################################################################
#Input is a wave file named 'w', time expansion factor, list of detection data for each freq band (w.list)
filtr_clicks <- function(w,exp,w.list,LPF.list,clickpercent, channel="left", fn=fn) {

  #Load model
  system.file("inst/qda_modelclick(shortfiles).rds", package="buzzfindr")
  model <- readRDS("inst/qda_modelclick(shortfiles).rds")

  #load function to use when subsetting the non-click data form the detected signals in buzz_findr
  fn <- function(x,y){
    if(x >= (y-0.001) && x <=(y+0.003)){ #Can play around with these values
      NA} else {
        x}
  }
  fn <- Vectorize(fn)

  #List of frequency filters
  LPF <- seq(13000, 50000, by=2000)
  HPF <- seq(12000, 49000, by=2000)

  #detect signals
  test <- NULL

for (z in 1:length(LPF)){

  capture.output(t <- threshold_detection(w,
                                          channel = channel,
                                          min_dur=0.5,
                                          min_TBE=2,
                                          threshold=4,
                                          LPF=LPF[z],
                                          HPF=HPF[z],
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
  t$data$event_data$freqband <- z
  if (is.null(ncol(t$data$event_data))==FALSE){
    test <- rbind(test,t$data$event_data)
  } else {
    next
  }

}

#Change starting time to numeric
tstart <- test$starting_time
tstart <- substr(tstart, start = 7, stop = 12)
test$starting_time <- tstart
test$starting_time <- as.character(test$starting_time)
test$starting_time <- as.numeric(test$starting_time)


#simplify the dataframe
test <- test[,c(2,29,7,8,10,20,3,14,16,17)]

#order by time
test <- test[c(order(test$starting_time)),]

#change all times within 1 msec to the same time
consecs <- runner::runner(
  x = test$starting_time,
  k = 2,
  f = function(x) {
    h <- x[2]-x[1]
    if(isTRUE(h < 0.002) | isTRUE(h == 0)){
      x[1]
    } else {
      x[2]
    }
  }
)

test$consecs <- consecs
test$consecs[1] <- 0

#Use aggregate to extract information
test$num <- 1
test$consecs <- as.factor(test$consecs)

#Aggregate to get variables
clickdata <- test %>%
  group_by(consecs) %>%
  summarise(num = sum(num),
            min_freq = min(freq_start),
            max_freq = max(freq_start),
            var_amp = var(bin_max_amp),
            avg_pcfreqmax = mean(pc_freq_max),
            avg_pcfreqmin = mean(pc_freq_min))

#bandwidth
clickdata$bandwidth <- clickdata$max_freq-clickdata$min_freq

#SUBSET BEFORE MODELLING
#Remove single detections
clickdata.sub1 <- clickdata[c(clickdata$num>11),]

#MORE SUBSET
#Remove single detections
clickdata.sub3 <- clickdata[c(clickdata$num>5) & c(clickdata$num<=9),]
clickdata.sub3 <- clickdata.sub3[c(clickdata.sub3$bandwidth>30),]
clickdata.sub3 <- clickdata.sub3[c(clickdata.sub3$min_freq>=0),]
clickdata.sub3 <- clickdata.sub3[c(clickdata.sub3$min_freq<16000),]

#USE LDA MODEL
#Subset for the model
clickdata.sub2 <- clickdata[c(clickdata$num<=9),]
if (nrow(clickdata.sub2)>0){
#Make predictions on test data
predictions <- model %>% predict(clickdata.sub2)
#just take predicted clicks
#Add the prediction to the dataframe as a variable
clickdata.sub2$result <- predictions$posterior[,2]
#Change the prediction to a grouping variable
clickdata.sub2$group <- NA
clickdata.sub2$group[clickdata.sub2$result>=(1-clickpercent)] <- "click"
#subset just clicks
clickdata.sub2 <- clickdata.sub2[c(is.na(clickdata.sub2$group)==FALSE),]
#remove new comlumns and combine with above clickdata
clickdata.sub2 <- clickdata.sub2[,-c(9,10)]
clickdata <- rbind(clickdata.sub1, clickdata.sub2, clickdata.sub3)
} else {
  clickdata <- clickdata.sub1
}

#get vector of detections
cd <- clickdata$consecs
cd <- as.vector(cd)
cd <- as.character(cd)
cd <- as.numeric(cd)

#Removing the click time data
if (length(cd)>0){
for (s in 1:length(LPF.list)){
  truedets <- NULL

    if (length(w.list[[s]]$event_data$starting_time)==0){
    next
  }

  x <- w.list[[s]]$event_data$starting_time

  test <- as.data.frame(outer(x, cd, FUN = fn))
  truedets <- as.data.frame(test[c(complete.cases(test)),])
  truedets <- truedets[,1]


  w.list[[s]]$event_data <- w.list[[s]]$event_data[c(w.list[[s]]$event_data$starting_time %in% truedets),]
    }
} else {
  }
 return(w.list)
}



