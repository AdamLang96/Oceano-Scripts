library(rsample)      
library(gbm)          
library(caret)        

data <- read.csv("/Users/adamlang/Desktop/airvisualcollocation2.csv", 
                 as.is = TRUE, # to keep dates as characters, not factors
                 header = FALSE, 
                 skip = 3)

names(data) <- c("Date",
                 "BAMPM10", 
                 "AVAmbientTemp", 
                 "AVPM1", 
                 "AVPM10",
                 "AVPM2.5",
                 "AVRelHum")

data$Date <- as.POSIXct(data$Date, tz = "Etc/GMT+8", format = "%d-%b-%Y %H:%M")

oceanoAV.model <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                distribution = 'gaussian',
                data = rtree_train,
                n.trees = 5000,
                interaction.depth = 3,
                shrinkage = 0.1,
                cv.folds = 5,
                n.cores = NULL,
                verbose = FALSE)

preds.oceanoAV <- predict(oceanoAV.model)
plot(rtree_train$BAMPM10, preds.oceanoAV)
abline(0, 1)
sqrt(min(oceanoAV.model$cv.error))

### 4^4=256 trees computationally expensive
hyper_grid <- expand.grid(
  shrinkage = c(.005,.01, .05, .1),
  interaction.depth = c(1, 3, 5, 7),
  n.minobsinnode = c(15, 20, 25, 30),
  bag.fraction = c(.5, .65, .8, 1), 
  optimal_trees = 0,
  min_RMSE = 0)  



#####tuning over grid of hyper parameters
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(123)
  
  
  # train model
  gbm.tune <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                  distribution = 'gaussian',
                  data = rtree_train,
                  n.trees = 5000,
                  interaction.depth = hyper_grid$interaction.depth[i],
                  shrinkage = hyper_grid$shrinkage[i],
                  n.minobsinnode = hyper_grid$n.minobsinnode[i],
                  bag.fraction = hyper_grid$bag.fraction[i],
                  train.fraction = .75,
                  n.cores = NULL,
                  verbose = FALSE)
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

##ordering by least RMSE
hyper_grid <- hyper_grid[order(hyper_grid$min_RMSE), ]

##top ten results
hyper_grid[1:10,]

#fitting best model
bestfit.oceano <- gbm(formula = BAMPM10 ~ AVPM10 + AVRelHum + AVAmbientTemp,
                      distribution = 'gaussian',
                      data = rtree_train,
                      n.trees = 5000,
                      interaction.depth = 5,
                      shrinkage = 0.1,
                      n.minobsinnode = 25,
                      bag.fraction = .5,
                      train.fraction = 1,
                      cv.folds = 5,
                      n.cores = NULL,
                      verbose = FALSE)


##RMSE on test data very similar to training data
#8.725 on test vs 8.748 on train
bestfit.pred <-predict(bestfit.oceano, n.trees = bestfit.oceano$n.trees, rtree_test)
caret::RMSE(bestfit.pred, rtree_test$BAMPM10)
plot(rtree_test$BAMPM10,bestfit.pred)
abline(0, 1)



