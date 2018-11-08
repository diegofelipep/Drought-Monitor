library(lfstat)
library(RcmdrPlugin.lfstat)
library(ggplot2)
library(forecast)
library(dplyr)
library(lubridate)
library(tibble)
library(plyr)

#################### Part 1

files <- list.files("PATH")
files = files[1:(length(files)-1)]
#files

minima <- array(dim = length(files))
quantile70 <- array(dim = length(files))
quantile90 <- array(dim = length(files))
duration <- array(dim = length(files))
volume <- array(dim = length(files))

for(i in 1:length(files)){
  data <- read.table(file = paste("PATH", files[i], sep = ""), header = F)
 
  timeseries = ts(data[,5]) # convert flow to time series
  startdate = as.Date(paste(data[1,2], "/", data[1,3], "/", data[1,4], sep = "")) # set start date
  
  obj = createlfobj(timeseries, hyearstart = 11, startdate = startdate) # convert data to low flow object
  flowunit(obj) = "m³/s" # set unit of time series
  
  fdc(obj,yearly = T, legend = T, separate=T)
 
  min = MAM(obj, n = 7, yearly = T) # compute annual minima, adjust n
  minima[i] = min$hyear[which.min(min$MAn)] 

  quant1 = Q70(obj, yearly = T) # compute quantiles, adjust
  quantile70[i] = quant1$hyear[which.min(quant1$flow)] 
  
  quant2 = Q90(obj, yearly = T) # compute quantiles, adjust
  quantile90[i] = quant2$hyear[which.min(quant2$flow)] 
  
  finddroughts = find_droughts(obj,threshold= "Q70")#, threshold = ..., varying = ...) compute droughts
  droughts = pool_it(finddroughts,tmin = 7) # pool droughts
  sum = summary(droughts)
  dur = sum$time[which.max(sum$duration)]
  vol = sum$time[which.max(sum$volume)]
  duration[i] = substr(dur,1,4)
  volume[i] = substr(vol,1,4)
  
  report <- lfstat:::streamdefRcmdr(lfobj = obj, pooling = 'MA', threslevel = , thresbreaks = 'fixed', breakdays = c(), MAdays = 7, tmin = 5,IClevel = 0.1, mindur = '0', minvol = '0', table = 'all', plot = 1)
  write.csv(report, paste("PATH",files[i],".csv"))
 }

volume = as.numeric(volume)
duration = as.numeric(duration)

histmin <- hist(minima, breaks = max(minima, na.rm = T) - min(minima, na.rm = T) + 2)
histq70 <- hist(quantile70, breaks = max(quantile70, na.rm = T) - min(quantile70, na.rm = T) + 2)
histq90 <- hist(quantile90, breaks = max(quantile90, na.rm = T) - min(quantile90, na.rm = T) + 2)
histvol <- hist(volume, breaks = max(volume, na.rm = T) - min(volume, na.rm = T) + 2)
histdur <- hist(duration, breaks = max(duration, na.rm = T) - min(duration, na.rm = T) + 2)


#################### Part 2


for(i in 1:length(files)){
  data <- read.table(file = paste("PATH", files[i], sep = ""), header = F)
  dataselect <- data[data[,2] %in% 1975:1977,]
  
  if(nrow(dataselect) == 0) next
  
  Station <- array(dim = nrow(dataselect))
  Year <- array(dim=nrow(dataselect))
  Month <- array(dim=nrow(dataselect))
  Day <- array(dim=nrow(dataselect))
  Q <- array(dim = nrow(dataselect))
  deficit <- array(dim = nrow(dataselect))
  duration <- array(dim = nrow(dataselect))
  sumdeficit <- array(dim = nrow(dataselect))
  sumduration <- array(dim = nrow(dataselect))
  
  threshold <- quantile(data[,5], probs = 0.1, na.rm = T)
  
  for(ii in 1:nrow(dataselect)){
    
    Station[ii] = dataselect[ii,1]
    Year[ii] = dataselect[ii,2]
    Month[ii] = dataselect[ii,3]
    Day[ii] = dataselect[ii,4]
    Q[ii] = dataselect[ii,5]
    deficit[ii] = dataselect[ii,5] - threshold
  }
  
  deficit[deficit > 0] = 0
  duration = deficit
  duration[duration < 0] = 1
  
  
  sumdeficit <- ave(deficit,cumsum(deficit==0),FUN=cumsum)
  sumduration <- ave(duration,cumsum(duration==0),FUN=cumsum)
  
  
   cat("Station Year Month Day Q Deficit Duration SumDeficit SumDuration \n", file = paste("PATH", files[i], sep = ""), append = F)
  for(ii in 1:length(Q)){
    cat(Station[ii],Year[ii], Month[ii], Day[ii],Q[ii], deficit[ii], duration[ii], sumdeficit[ii], sumduration[ii], "\n", file = paste("PATH", files[i], sep = ""), append = T)
  }
}
  
#################### Part 3

# Variables that you can change to fit best to your purpose
pathTostations = "PATH"                  # the path where the stations files are located
start.date = c(1975,11,01)     # year, month,day
end.date   = c(1977,10,31)     # year, month,day, hour, minute
output.path = "PATH"		     # the path where you want the final table to be saved - together with the file name
time.resolution = 24*60		     # [min] in what temporal resolution is your data - i suppose in 5 min?
missing.value  = -999                # the value for recognizing missing timesteps

#-------------------- START OF THE COMPUTATIONS ----------------------------#

# reading the names of the stations
stations.fullnames = list.files(pathTostations,pattern="*.txt", full.names=T)
stations.names = list.files(pathTostations,pattern="*.txt", full.names=F)
no.of.stats = length(stations.names)

# converting the dates in the respective format into time, vector and character 
start.asDate = as.POSIXct(paste0(as.character(start.date), collapse="/"), format= "%Y/%m/%d", 
                          origin="1900-01-01 00:00", tz="UTC")
end.asDate = as.POSIXct(paste0(as.character(end.date), collapse="/"), format= "%Y/%m/%d", 
                        origin="1900-01-01 00:00", tz="UTC")
dates.Sequence = seq(start.asDate, end.asDate, time.resolution*60)
dates.Sequence.vec = do.call(rbind, lapply(dates.Sequence, function(i) c(year(i), month(i), day(i))))
dates.Sequence.char = sapply(dates.Sequence, function(i) paste0(year(i),"/", month(i), "/", day(i)))
no.timesteps = length(dates.Sequence)

# reading the data for each station according to the start and end time
# iterating over each station from 1 to no.of.stats
combined.Data = lapply(1:no.of.stats, function(stat){
  # reading the data at the stat-th station
  stat.data = read.table(stations.fullnames[stat], header=T)
  # finding the date indexes that correspond to the given start and end time
  start.index = which(stat.data[,2]==start.date[1] & stat.data[,3]==start.date[2] &
                        stat.data[,4]==start.date[3] )
  end.index = which(stat.data[,2]==end.date[1] & stat.data[,3]==end.date[2] &
                      stat.data[,4]==end.date[3])
  # first condition: if the desired start and end time are not present in the data then it returns missing values (here notes as -999)
  if(length(start.index)==0 | length(end.index)==0) rain.stat = as.data.frame(matrix(missing.value, nrow = no.timesteps, ncol=6)) ####
  # otherwise if the start and date time are present in the station data:
  if(length(start.index)>0 & length(end.index)>0){
    cut.data = stat.data[start.index:end.index,]
    # second condition are there any missing timesteps within the selected data? If yes the missing timesteps information is noted as -999 
    if(dim(cut.data)[1]== no.timesteps) rain.stat = cut.data[,c(1,5:9)]
    else{
      cut.in.dates = apply(cut.data[,2:4], function(i) collapse="/")
      rain.stat = rep(missing.value, no.timesteps)
      rain.stat[which(dates.Sequence.char %in% cut.in.dates)] <- cut.data[which(cut.in.dates%in% dates.Sequence.char),c(1,5:9)]
    }
  }
  names(rain.stat) = c("Station" ,   "Q"     ,     "Deficit" ,   "Duration"   ,"SumDeficit", "SumDuration")
  rain.stat[,1] = stat.data[1,1]
  print(paste(stations.fullnames[stat], "has been read!"))
  # return only the precipitation vector for the given duration 
  return(rain.stat)
})

# combining the dates information together with the Discharge from each station in a dataframe		
Final.Data = cbind(dates.Sequence.vec, do.call(rbind, combined.Data))
Final.Data = as.data.frame(Final.Data)
names(Final.Data) = c("Year", "Month", "Day", "Station", "Q", "Deficit", "Duration", "SumDeficit", "SumDuration")

for(i in 1:length(dates.Sequence.char)){
  subset.daily = which(Final.Data[,1]==dates.Sequence.vec[i,1] & Final.Data[,2]==dates.Sequence.vec[i,2]
                       & Final.Data[,3]==dates.Sequence.vec[i,3])
  daily.data = Final.Data[subset.daily,]
  output.filename = paste0(output.path,daily.data[1,1], "_",daily.data[1,2],"_",daily.data[1,3],".txt")
  write.table(daily.data, file = output.filename, quote=F, row.names = F, sep="\t")
  print(i)
}








