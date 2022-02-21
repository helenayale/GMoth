library(ggplot2)
library(sp)
library(rgdal)
library(raster)
library(gtools)
library(stringr)
library(randomForest)
library(dplyr)
library(caret)
library(sf)


setwd('D:/Haiyin/Final')

## read files

## responsible variable

basic <- read.csv(file = "basic2800.csv", header = T, row.names = 1)
head(basic)
#x_em <- def[,names(def) %in% c( "plot","mean_em", "def2017", "def2018")]

## defoliation rate
def <- read.csv('def_final2800.csv')
y_def <- def[,names(def) %in% c( "plot","def_ac")]
head(y_def)


## predictor variable
# 
x_nmds1 <- read.csv('nmds1_final2800.csv')
x_nmds1 <- x_nmds1[, -1]
head(x_nmds1)
#x_nmds1$nmds1_abs <- abs(x_nmds1$nmds1_pred)
#x_nmds1 <- nmds1[,names(nmds1) %in% c( "plot","nmds1_pred","SR_pred")]
#write.csv(x_nmds1, 'nmds1_final2800.csv')


str <- read.csv('F5_final2800.csv')
x_str <- str[,names(str) %in% c( "plot","LiD_PR_CANUND_pred")]
x_str$plot <- basic$plot

# temperature
temp_mean <- read.csv('temp/temp_mean2800.csv', header = T, row.names = 1)
temp_min <- read.csv('temp/temp_min2800.csv', header = T, row.names = 1)
temp_max <- read.csv('temp/temp_max2800.csv', header = T, row.names = 1)
prec <- read.csv('temp/temp_prec2800.csv', header = T, row.names = 1)

temp_diff <- temp_max - temp_min
temp_diff <- temp_diff[,-1]

x_temp <- data.frame(plot = basic$plot)
#x_temp$t_mean <- rowMeans(temp_mean [,!names(temp_mean ) %in% c( "id")])
x_temp$t_mean <- rowMeans(temp_mean [,3:5])
x_temp$t_min <- apply(temp_min[, 3:5], 1, FUN = min)
x_temp$t_diff <- rowMeans(temp_diff[2:4])
x_temp$prec_mean <- rowMeans(prec[3:5])


#cater <- read.csv('Figure3_Caterpillars_Defoliation.csv')
#cater$plot <- cater$Plot
#x_cater <- cater[,names(cater) %in% c( "plot","Caterpillars")]


### merge preditor variables
pred <- basic
pred <- left_join(pred, x_nmds1, by = 'plot')  
pred <- left_join(pred, x_str, by = 'plot')
pred <- left_join(pred, x_temp, by = 'plot')
#pred <- left_join(pred, x_cater, by = 'plot')

pred <- na.omit(pred)


var_all <- left_join(pred, y_def, by = 'plot')

var_pred <- pred[,-1]

######################################################
################## random forest #####################
######################################################


colnames(var_pred)


v_pred <- var_pred[,!names(var_pred) %in% c("t_mean","t_diff")]
#v_pred <- var_pred[,names(var_pred) %in% c("nmds1_pred","LiD_PR_CANUND_pred", "SR_pred","def2017","def2018" )]
#v_pred <- var_pred[,!names(var_pred) %in% c("t_mean","t_diff")]
colnames(v_pred)


#############  compare models  ##############

library(foreach)
library(doParallel)


colnames(var_pred)

v_list <- list()

v_list[[1]] <- var_pred[,!names(var_pred) %in% c("t_min","t_mean")]
v_list[[2]] <- var_pred[,!names(var_pred) %in% c("t_diff","t_mean")]
v_list[[3]] <- var_pred[,!names(var_pred) %in% c("t_min","t_diff")]
v_list[[4]] <- var_pred[,names(var_pred) %in% c("LiD_PR_CANUND_pred","X", "Y", "def2017","def2018","nmds1_pred")]
v_list[[5]] <- var_pred[,names(var_pred) %in% c("LiD_PR_CANUND_pred","X", "Y", "def2017","def2018","SR_pred")]
v_list[[6]] <- var_pred[,names(var_pred) %in% c("LiD_PR_CANUND_pred","X", "Y", "def2017","def2018")]


#name_list <- c('t_diff', 't_min', 't_mean','nmds','SR','XY')

# A grid can be generated to specify candidate hyper-parameter values for inclusion into the models training
rf.grid <- expand.grid(mtry=1:24) # number of variables available for splitting at each tree node, can be adjusted to improve model

# Set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                   number = 10, # number of folds
                   repeats = 5, # number of repeats
                   allowParallel = TRUE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE, #print training log
                   p=0.7) # the training percentage

modellist <- list()

#Define how many cores you want to use
UseCores <- detectCores()-1

#Register CoreCluster
cl <- makeCluster(UseCores)
registerDoParallel(cl)


modellist <- foreach(i = 1:6) %dopar% {
  library(caret)
  
  #train random forest model 
  set.seed(1) #make it reproducible
  
  fit <- caret::train(x = v_list[[i]], y = var_all[,'def_ac'],
                      method = "rf", metric="Rsquared", trainControl = tc, 
                      tuneGrid = rf.grid)
  #key <- name_list[i]
  return(fit)
}


#end cluster
stopCluster(cl)

# turn parallel processing off and run sequentially again:
registerDoSEQ()


# compare results
results <- resamples(modellist)

summary(results)

modellist[[1]]$finalModel
modellist[[2]]$finalModel
modellist[[3]]$finalModel
modellist[[4]]$finalModel
modellist[[5]]$finalModel
modellist[[6]]$finalModel

max(modellist[[3]]$results$Rsquared)
dotplot(results, metric = 'Rsquared')
dotplot(results, metric = 'RMSE')
dotplot(results, metric = 'MAE')
min(modellist[[3]]$results$RMSE)

modellist[[3]]$bestTune
importance(modellist[[3]]$finalModel)
varImpPlot(modellist[[3]]$finalModel, main = "n=2800", pch =17, cex = 1.5, n.var = 10, color = "mediumpurple")


max(modellist[[3]]$results$Rsquared)
min(modellist[[3]]$results$RMSE)

modellist[[3]]

for (i in 1:6){
  print(max(modellist[[i]]$results$Rsquared))
  print(min(modellist[[i]]$results$RMSE))
}

################ single model ##################


# A grid can be generated to specify candidate hyper-parameter values for inclusion into the models training
rf.grid <- expand.grid(mtry=1:24) # number of variables available for splitting at each tree node, can be adjusted to improve model

# Set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                   number = 10, # number of folds
                   repeats = 5, # number of repeats
                   allowParallel = TRUE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE, #print training log
                   p=0.7) # the training percentage

##############

#train random forest model 
set.seed(1) #make it reproducible
rf_model<- caret::train(x = v_pred, y = var_all[,'def_ac'],
                        method = "rf", metric="Rsquared", trainControl = tc, 
                        tuneGrid = rf.grid)

rf_model$results
max(rf_model$results$Rsquared)
rf_model$bestTune


rf_model$finalModel

importance(rf_model$finalModel)
varImpPlot(rf_model$finalModel)



#######################
##### Spearman's Correlation


# Spearman correlation 1+2+3
# PR ~ def
pv <- cor.test(var_all[,'LiD_PR_CANUND_pred'], var_all[,'def_ac'], na.action = "na.exclude", method="spearman")
if (pv$p.value <0.001) {pvalue <- "***"} else if (pv$p.value <0.05) {pvalue <- "**"} else if (pv$p.value <0.01) {pvalue <- "*"}else(pvalue <- "")
pv
pvalue


# Linear Regression
# Apply the lm() function
y <- var_all[,'def_ac']
x <- var_all[,'LiD_PR_CANUND_pred']

relation <- lm(y ~ x)
#print(relation)
print(summary(relation))

# multi-lm
x1 <- var_all[,'LiD_PR_CANUND_pred']
x2 <- var_all[,'def2018']
x3 <- var_all[,'def2017']
x4 <- var_all[,'Y']
x5 <- var_all[,'X']
x6 <- var_all[,'nmds1_pred']
x7 <- var_all[,'SR_pred']
x8 <- var_all[,'t_mean']
x9 <- var_all[,'prec_mean']
x10 <- var_all[,'mean_em']

relation <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10)
#print(relation)
print(summary(relation))

