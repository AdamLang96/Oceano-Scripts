########## Libraries #######################

library(chron)
library(dplyr)
library(gbm)          
library(caret)     
load("mymodel.rda")



############# Script #####################

url <- "https://storage.googleapis.com/slocleanair-org/pages/air-quality/OceanoCDF.CSV"
df <- suppressWarnings(OceanoData(url)) %>%
  DataCorrection



############## Functions ###################


OceanoData <- function(url) {
   # Takes in googleapi url for the 3 Air Visual Sensors in Oceano
   # and returns list of 3 dfs pertaining to each sensor with cols:
   #    Site
   #    Date
   #    AVPM10        (PM10 values)
   #    AVAmbientTemp (ambient temperature)
   #    AVRelHum      (relative humidity)
   #    Latitude
   #    Longitude
   #
   #
   # Args:
   #   "https://storage.googleapis.com/slocleanair-org/pages/air-quality/OceanoCDF.CSV"
   #
   # Returns:
   #   List of 3 dataframes
   #
  data <- read.csv(url,
                   strip.white = TRUE)
  
  #Cleaning and isolating CSD data
  OceanoCSD <- data.frame("Site"=rep("OceanoCSD", length(nrow(data))),
                          "Date"          = data$X, 
                          "AVPM10"        = data$Oceano.CSD.7, 
                          "AVAmbientTemp" = data$Oceano.CSD.3, 
                          "AVRelHum"      = data$Oceano.CSD.9,
                          "Latitude"      = data$Oceano.CSD.4,
                          "Longitude"     = data$Oceano.CSD.5)
  
  OceanoCSD$Date    <- as.POSIXct(OceanoCSD$Date, 
                                  tz = "Etc/GMT+8", 
                                  format = "%d-%b-%Y %H:%M")
  
  OceanoCSD$AVPM10        <- as.numeric(as.character(OceanoCSD$AVPM10))
  OceanoCSD$AVAmbientTemp <- as.numeric(as.character(OceanoCSD$AVAmbientTemp))
  OceanoCSD$AVRelHum      <- as.numeric(as.character(OceanoCSD$AVRelHum))
  OceanoCSD               <- na.omit(OceanoCSD)
  rownames(OceanoCSD) <- seq(length=nrow(OceanoCSD))
  OceanoCSD <- OceanoCSD[order(nrow(OceanoCSD):1), ]
  
  #Cleaning and isolating 22nd St. data
  Oceano22nd<-data.frame("Site"          =rep("Oceano22nd", length(nrow(data))),
                         "Date"          = data$X,
                         "AVPM10"        = data$Oceano.22nd.St.4,
                         "AVAmbientTemp" = data$Oceano.22nd.St,
                         "AVRelHum"      = data$Oceano.22nd.St.6,
                         "Latitude"      = data$Oceano.22nd.St.1,
                         "Longitude"     = data$Oceano.22nd.St.2)
  Oceano22nd$Date    <- as.POSIXct(Oceano22nd$Date, 
                                   tz = "Etc/GMT+8", 
                                   format = "%d-%b-%Y %H:%M")
  
  Oceano22nd$AVPM10        <- as.numeric(as.character(Oceano22nd$AVPM10))
  Oceano22nd$AVAmbientTemp <- as.numeric(as.character(Oceano22nd$AVAmbientTemp))
  Oceano22nd$AVRelHum      <- as.numeric(as.character(Oceano22nd$AVRelHum))
  Oceano22nd               <- na.omit(Oceano22nd)
  rownames(Oceano22nd)     <- seq(length=nrow(Oceano22nd))
  Oceano22nd <- Oceano22nd[order(nrow(Oceano22nd):1), ]
  
  OceanoEast<-data.frame("Site" = rep("OceanoEast", length(nrow(data))),
                         "Date" = data$X,
                         "AVPM10" = data$Oceano.East.4,
                         "AVAmbientTemp" = data$Oceano.East,
                         "AVRelHum" = data$Oceano.East.6,
                         "Latitude" = data$Oceano.East.1,
                         "Longitude" = data$Oceano.East.2)
  
  OceanoEast$Date          <- as.POSIXct(OceanoEast$Date, 
                                         tz = "Etc/GMT+8", 
                                         format = "%d-%b-%Y %H:%M")
  OceanoEast$AVPM10        <- as.numeric(as.character(OceanoEast$AVPM10))
  OceanoEast$AVAmbientTemp <- as.numeric(as.character(OceanoEast$AVAmbientTemp))
  OceanoEast$AVRelHum      <- as.numeric(as.character(OceanoEast$AVRelHum))
  OceanoEast               <- na.omit(OceanoEast)
  rownames(OceanoEast)     <- seq(length=nrow(OceanoEast))
  OceanoEast <- OceanoEast[order(nrow(OceanoEast):1), ]
  
  
  AllSensors <-list( OceanoCSD, 
                     Oceano22nd, 
                     OceanoEast)
  return(AllSensors)
}



################################


DataCorrection <- function(datalist) {
  # Takes in output from OceanoData
  # corrects PM10 values based on regression model
  # and calculates AQI based off corrected data
  #
  # Args:
  #  Output from OceanoData
  #
  # Returns:
  #  Dataframe with cols:
  #    Site (OceanoCSD, Oceano 22nd, OceanoEast)
  #    Date
  #    AQI
  #    Latitude
  #    Longitude
  data<- data.frame("Site" = NA,
                    "Date" = NA,
                    "AQI"  =  NA,
                    "Latitude"  = NA,
                    "Longitude" = NA,
                    stringsAsFactors = FALSE)
  for(i in 1:length(datalist)) {
    dataframe     <- datalist[[i]]
    correcteddata <- predict(oceanofullmodel, newdata = dataframe)
    AQI <- AQICalc(correcteddata)
    vec <- data.frame( "Site" = as.character(dataframe$Site[1]),
                       "Date" = dataframe$Date[1],
                       "AQI"  = AQI,
                       "Latitude"  = as.numeric(as.character(dataframe$Latitude[1])),
                       "Longitude" = as.numeric(as.character(dataframe$Longitude[1])))
    data<-rbind(data,vec)
    attributes(data$Date) <- attributes(dataframe$Date[1])
    data <- na.omit(data)
    rownames(data) <- seq(length = nrow(data))
 }
return(data)
}

###################################

AQICalc <- function(vect) {
  # Takes a vector of PM10 values and
  # calculates the 12 hour NowCast AQI
  # according to EPA specifications
  #
  # Args:
  #  vector of PM10 values
  #
  # Returns:
  #   NowCast AQI
  #
  num <- 0
  dem <- 0
  PM10 <- vect[1:12]
  if(sum(is.na(PM10)) >= 3) {
    return(NA)
  } else {
    for(i in 1:12) { 
      if(is.na(PM10[i])) {
        PM10[i] <- 0
      }
    }
  }
  WeightFactor <- 1 - (min(PM10) / max(PM10))
  if( WeightFactor < 1/2) {
    WeightFactor <- 1/2
  }
  NowCastConc <- for(i in 1:12) {
    num <- num + (PM10[i] * WeightFactor ^ (i-1))
    dem <- dem + (WeightFactor ^ (i-1))
  }
  AQI <- Converter(num / dem)
  return(AQI)
}

################################################

Converter <- function(x, table = AQILookUpTable) {
  #Takes in a NowCast concentration and the AQIlookuptable
  #Returns the NowCast AQI value
  #Used in AQICalc
  #
  #Args:
  # Nowcast AQI concentration for PM10 and AQIlookuptable
  #
  #Returns:
  # NowCastAQI
  stopifnot(is.numeric(x), length(x) == 1)
  if(is.na(x)) return (NA)
  AQILo  <- table$AQILo[max(which(x >= table$pm10))]
  AQIHi  <- table$AQIHi[max(which(x >= table$pm10))]
  ConcLo <- table$pm10[max(which(x >= table$pm10))]
  ConcHi <- table$pm10[max(which(x >= table$pm10)) + 1] - 1
  AQI    <- round((((AQIHi-AQILo) / (ConcHi-ConcLo)) * (x-ConcLo))+AQILo)
  return(AQI)
}


##############################################################

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