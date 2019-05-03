library(rsample)      
library(gbm)          
library(caret)        

oceanodata1 <- read.csv("https://storage.googleapis.com/slocleanair-org/pages/air-quality/OceanoCDF.CSV", 
                       stringsAsFactors = FALSE)
CSDdata1 <- data.frame(oceanodata1$X, 
                      oceanodata1$Oceano.CSD.2, 
                      oceanodata1$Oceano.CSD.7, 
                      oceanodata1$Oceano.CSD.3, 
                      oceanodata1$Oceano.CSD.9)
names(CSDdata1) <- c("Date","BAMPM10", "AVPM10", "AVAmbientTemp", "AVRelHum")
CSDdata1$Date <- as.POSIXct(CSDdata1$Date, tz = "Etc/GMT+8", format = "%d-%b-%Y %H:%M")
CSDdata1$BAMPM10 <- as.numeric(as.character(CSDdata1$BAMPM10))
CSDdata1$AVPM10 <- as.numeric(as.character(CSDdata1$AVPM10))
CSDdata1$AVAmbientTemp <- as.numeric(as.character(CSDdata1$AVAmbientTemp))
CSDdata1$AVRelHum <- as.numeric(as.character(CSDdata1$AVRelHum))
CSDdata1 <- na.omit(CSDdata1)

oceanodata2 <- read.csv("/Users/adamlang/Desktop/OceanoAVData.csv")
CSDdata2 <- data.frame("Date" = oceanodata2$X, 
                       "AVAmbientTemp" = oceanodata2$Oceano.CSD.2,
                       "AVRelHum" = oceanodata2$Oceano.CSD.5,
                       "AVPM10" = oceanodata2$Oceano.CSD.3,
                       "BAMPM10" = oceanodata2$Oceano.CSD)
CSDdata2$Date <- as.POSIXct(CSDdata2$Date, tz = "Etc/GMT+8", format = "%d-%b-%Y %H:%M")
CSDdata2$BAMPM10 <- as.numeric(as.character(CSDdata2$BAMPM10))
CSDdata2$AVPM10 <- as.numeric(as.character(CSDdata2$AVPM10))
CSDdata2$AVAmbientTemp <- as.numeric(as.character(CSDdata2$AVAmbientTemp))
CSDdata2$AVRelHum <- as.numeric(as.character(CSDdata2$AVRelHum))
CSDdata2 <- na.omit(CSDdata2)

CSDdata <- rbind(CSDdata2, CSDdata1)

mtrain <- CSDdata[1:160, ]
mtest <- CSDdata[161:222, ]

oceanoAVnew.model1 <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                         distribution = 'gaussian',
                         data = mtrain,
                         interaction.depth = 3,
                         n.trees = 500,
                         shrinkage = 0.1,
                         cv.folds = 5,
                         n.cores = NULL,
                         verbose = FALSE)

sqrt(min(oceanoAVnew.model1$cv.error))
newtestpred <- predict(oceanoAVnew.model1, n.trees = oceanoAVnew.model1$n.trees, mtest)
caret::RMSE(newtestpred, mtest$BAMPM10)


oceanohypergrid <- expand.grid(
  shrinkage = c(.005,.01, .05, .1),
  interaction.depth = c(1, 3, 5, 7),
  n.minobsinnode = c(5, 10, 15, 20),
  bag.fraction = c(.5, .65, .8, 1), 
  optimal_trees = 0,
  min_RMSE = 0)  

for(i in 1:nrow(oceanohypergrid)) {
  
  # reproducibility
  set.seed(123)
  
  
  # train model
  oceanonew.tune <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                        distribution = 'gaussian',
                        data = mtrain,
                        n.trees = 500,
                        interaction.depth = oceanohypergrid$interaction.depth[i],
                        shrinkage = oceanohypergrid$shrinkage[i],
                        n.minobsinnode = oceanohypergrid$n.minobsinnode[i],
                        bag.fraction = oceanohypergrid$bag.fraction[i],
                        train.fraction = 1,
                        n.cores = NULL,
                        verbose = FALSE)
  
  # add min training error and trees to grid
  oceanohypergrid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  oceanohypergrid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

##ordering by least RMSE
oceanohypergrid <- oceanohypergrid[order(oceanohypergrid$min_RMSE), ]

##top ten results
oceanohypergrid[1:10,]

oceanoAVnew.modelbest <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                          distribution = 'gaussian',
                          data = mtrain,
                          interaction.depth = 1,
                          n.trees = 500,
                          shrinkage = 0.1,
                          n.minobsinnode = 5,
                          bag.fraction = 0.5,
                          cv.folds = 5,
                          n.cores = NULL,
                          verbose = FALSE)

sqrt(min(oceanoAVnew.modelbest$cv.error))
newtestpredbest <- predict(oceanoAVnew.modelbest, n.trees = oceanoAVnew.modelbest$n.trees, mtest)
caret::RMSE(newtestpredbest, mtest$BAMPM10)


plot(mtest$BAMPM10,newtestpredbest, 
     xlim = c(0,105),
     ylim = c(0,105),
     xlab = "BAMPM10",
     ylab = "Predicted Values")

abline(0,1)


FullData<- read.csv("/Users/adamlang/Desktop/BasicDataExportReport0426.csv")
names(FullData) <- c("Date",
                 "BAMPM10", 
                 "AVAmbientTemp", 
                 "AVPM1", 
                 "AVPM10",
                 "AVPM2.5",
                 "AVRelHum")
FullData<-FullData[3:nrow(FullData), ]
FullData$AVPM1<-NULL
FullData$AVPM2.5<-NULL
FullData$Date<-as.POSIXct(FullData$Date, tz = "Etc/GMT+8", format = "%d-%b-%Y %H:%M")
FullData$BAMPM10<-as.numeric(as.character(FullData$BAMPM10))
FullData$AVAmbientTemp<-as.numeric(as.character(FullData$AVAmbientTemp))
FullData$AVPM10<-as.numeric(as.character(FullData$AVPM10))
FullData$AVRelHum<-as.numeric(as.character(FullData$AVRelHum))
FullData<-na.omit(FullData)
FullTrain<-FullData[1:380, ]
FullTest<-FullData[381:505, ]



oceanofullmodel <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                             distribution = 'gaussian',
                             data = FullTest,
                             interaction.depth = 1,
                             n.trees = 500,
                             shrinkage = 0.1,
                             n.minobsinnode = 5,
                             bag.fraction = 0.5,
                             cv.folds = 5,
                             n.cores = NULL,
                             verbose = FALSE)

sqrt(min(oceanofullmodel$cv.error))
fulltestpred <- predict(oceanofullmodel, n.trees = oceanofullmodel$n.trees, FullTest)
caret::RMSE(fulltestpred, FullTest$BAMPM10)

save(oceanofullmodel, file = "mymodel.rda")

