library('chron')
library('dplyr')


OceanoMain <- function(csvurl) {
#Main function
#takes in CSV and returns NowCast AQI for Oceano based on Air Visual PM10 values
#
#Args:
#URL for CSV file containing:
   #Date/Time
   #Oceano BAM moniter PM10 Values
   #Air Visual Ambient Temperature Values
   #Air Visual PM10 Values
   #Air Visual PM2.5 Values
   #Air Visual Relative Humidity Values
#  
#Returns:
#12 Hour Nowcast AQI for PM10  
#NA if criteria isnt met for computing 12-Hour-Nowcast AQI
  data        <- read.csv(csvurl,
                          as.is = TRUE, # to keep dates as characters, not factors
                          header = FALSE, 
                          skip = 3)
  
  names(data) <- c("Date","BAMPM10", "AVAmbientTemp", "AVPM1", "AVPM10","AVPM2.5","AVRelHum")
  data$Date   <-  as.POSIXct(data$Date, tz = "Etc/GMT+8", format = "%d-%b-%Y %H:%M")
  data        <- data[nrow(data):1, ]
  conc        <- OceanoPM10AQI(data)
  if(is.na(conc)) {
  print("PM10 Nowcast AQI not available")
  }
  return(conc)
}



#DataFrame for comparing NowCast concentrations
#to AQI values
#
# Used in PM10ConcToAQInew
AQILookUpTable <- data.frame(ozone = c(0, 55, 71, 86, 106, 201, NA, NA), 
                             pm2.5 = c(0, 12.1, 35.5, 55.5, 150.5, 250.5, 350.5, 500.4),
                             pm10  = c(0, 55, 155, 255, 355, 425, 505, 604),
                             so2   = c(0, 36, 76, 186, 305, 605, 805, 1004),
                             AQILo = c(0, 51, 101, 151, 201, 301, 401, NA),
                             AQIHi = c(50, 100, 150, 200, 300, 400, 500, NA))




Converter <- function(x, table = AQILookUpTable) {
  #Takes in a NowCast concentration and the AQIlookuptable
  #Returns the NowCast AQI value
  #
  #Args:
  # Nowcast AQI concentration for PM10 and AQIlookuptable
  #
  #Returns:
  # 
  stopifnot(is.numeric(x), length(x) == 1)
  if(is.na(x)) return (NA)
  AQILo  <- table$AQILo[max(which(x >= table$pm10))]
  AQIHi  <- table$AQIHi[max(which(x >= table$pm10))]
  ConcLo <- table$pm10[max(which(x >= table$pm10))]
  ConcHi <- table$pm10[max(which(x >= table$pm10)) + 1] - 1
  AQI    <- round((((AQIHi-AQILo) / (ConcHi-ConcLo)) * (x-ConcLo))+AQILo)
  return(AQI)
}


OceanoPM10AQI<-function (df) {
  #Takes a dataframe of AirVisual PM10 values from the Oceano sensor and returns the 12-hour NowCast AQI score
  #Returns NA if criteri is not met for calculating NowCast AQI
  #Args:
  # A single datatable of PM10 by location of sensor
  #
  #Returns:
  # PM10 NowCast AQI and location of sensor
  # "Nowcast not available" if 2/3 of most recent hours are NA values
  num <- 0
  dem <- 0
  LastTwelveHours <- df[1:12, ]
  PM10ConcVar <- LastTwelveHours[["AVPM10"]]
  if(sum(is.na(PM10ConcVar)) >= 3) {
    return(NA)
  } else {
    for(i in 1:12) { 
      if(is.na(df[["AVPM10"]][i] == TRUE)) {
        PM10ConcVar[i] <- 0
      }
    }
  }
  WeightFactor <- 1 - (min(PM10ConcVar) / max(PM10ConcVar))
  if( WeightFactor < 1/2) {
    WeightFactor <- 1/2
  }
  NowCastConc <- for(i in 1:12) {
    num <- num + PM10ConcVar[i] * WeightFactor ^ (i-1)
    dem <- dem + WeightFactor ^ (i-1)
  }
  AQI <- Converter(num / dem)
  return(AQI)
}



