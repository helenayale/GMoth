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

basic <- read.csv('basic.csv')
head(basic)
basic <- basic[, -1]


# load caterpillar data
cater_file <- read.csv(file = 'Figure3_Caterpillars_Defoliation.csv')
cater <- data.frame(plot = cater_file$Plot, caterpillars = cater_file$Caterpillars)
head(cater)

## defoliation rate
def <- read.csv('def_final.csv')
y_def <- def[,names(def) %in% c( "plot","def_new","def_old")]

## predictor variable
nmds1 <- read.csv('nmds1_final.csv')
nmds1
x_nmds1 <- nmds1[,names(nmds1) %in% c( "plot","nmds1_pred", "SR_pred")]

str <- read.csv('F5_final.csv')
x_str <- str[,names(str) %in% c( "plot","LiD_PR_CANUND_pred")]



# temperature
temp_mean <- read.csv('temp/temp_mean.csv', header = T, row.names = 1)
temp_min <- read.csv('temp/temp_min.csv', header = T, row.names = 1)
temp_max <- read.csv('temp/temp_max.csv', header = T, row.names = 1)
prec <- read.csv('temp/temp_prec.csv', header = T, row.names = 1)

temp_diff <- temp_max[, -1] - temp_min[, -1]


x_temp <- data.frame(plot = basic$plot)
#x_temp$t_mean <- rowMeans(temp_mean [,!names(temp_mean ) %in% c( "id")])
x_temp$t_mean <- rowMeans(temp_mean [,3:5])
x_temp$t_min <- apply(temp_min[, 3:5], 1, FUN = min)
x_temp$t_diff <- rowMeans(temp_diff[2:4])
x_temp$prec_mean <- rowMeans(prec[3:5])





### merge preditor variables
pred <- basic
pred <- left_join(pred, cater, by = 'plot')
pred <- left_join(pred, x_nmds1, by = 'plot')
pred <- left_join(pred, x_str, by = 'plot')
pred <- left_join(pred, x_temp, by = 'plot')


pred <- na.omit(pred)


var_all <- left_join(pred, y_def, by = 'plot')

var_pred <- pred[,!names(pred) %in% c( "plot")]


######################################################
################## random forest #####################
######################################################

colnames(var_all)
colnames(var_pred)

#v_pred <- var_pred

v_pred <- var_pred[,!names(var_pred) %in% c("t_mean","t_diff")]
# v_pred <- var_pred[,!names(var_pred) %in% c("t_min","t_diff")]
# v_pred <- var_pred[,!names(var_pred) %in% c("t_mean","t_min")]
#v_pred <- var_pred[,!names(var_pred) %in% c("t_mean","t_diff","t_min", "X", "Y")]

# original 4 var
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars","nmds1_pred","t_mean", "LiD_PR_CANUND_pred")]

# replace nmds1 with SR
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars","SR_pred","t_mean", "LiD_PR_CANUND_pred")]

# add SR
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars","nmds1_pred" ,"SR_pred","t_mean", "LiD_PR_CANUND_pred")]

# add def2017, def2018
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars", "Def2017", "Def2018", "LiD_PR_CANUND_pred", "nmds1_pred", "t_mean")]

# add def2017, def2018, SR
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars", "Def2017", "Def2018", "LiD_PR_CANUND_pred", "nmds1_pred", "SR_pred", "t_mean")]

# add def1718, SR, SY
#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars", "Def2017", "Def2018", "LiD_PR_CANUND_pred", "nmds1_pred", "SR_pred", "t_mean", "X", "Y")]


#v_pred <- var_pred[,names(var_pred) %in% c("caterpillars", "LiD_PR_CANUND_pred")]

colnames(v_pred)


# A grid can be generated to specify candidate hyper-parameter values for inclusion into the models training
rf.grid <- expand.grid(mtry=1:24) # number of variables available for splitting at each tree node, can be adjusted to improve model

# Set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                   number = 10, # number of folds
                   repeats = 5, # number of repeats
                   allowParallel = TRUE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE, #print training log
                   p=0.7) # the training percentage

#########

#train random forest model 
set.seed(1) #make it reproducible
rf_model<- caret::train(x = v_pred, y = var_all[,'def_old'],
                        method = "rf", metric="Rsquared", trainControl = tc, 
                        tuneGrid = rf.grid)

rf_model$results
max(rf_model$results$Rsquared)
rf_model$bestTune


rf_model$finalModel

importance(rf_model$finalModel)
varImpPlot(rf_model$finalModel)


print(rf_model)
plot(rf_model)

max_rsq <- round(max(rf_model$finalModel$rsq),4)
max_rsq

RMSE <- min(sqrt(rf_model$finalModel$mse))
RMSE

min(rf_model$finalModel$mse)


#######################################################
################## compare models #####################
#######################################################
tune <- c('t_mean', 't_min', 't_diff')

v_list <- list()
v_list[['t_min']] <- var_pred[,!names(var_pred) %in% c("t_mean","t_diff")]
v_list[['t_mean']] <- var_pred[,!names(var_pred) %in% c("t_min","t_diff")]
v_list[['t_diff']] <- var_pred[,!names(var_pred) %in% c("t_mean","t_min")]

modellist <- list()


# A grid can be generated to specify candidate hyper-parameter values for inclusion into the models training
rf.grid <- expand.grid(mtry=1:24) # number of variables available for splitting at each tree node, can be adjusted to improve model

# Set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                   number = 10, # number of folds
                   repeats = 5, # number of repeats
                   allowParallel = TRUE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE, #print training log
                   p=0.7) # the training percentage


#train with different ntree parameters
for (n in tune){
  set.seed(1) #make it reproducible
  rf_model<- caret::train(x = v_list[[n]], y = var_all[,'def_old'],
                          method = "rf", metric="Rsquared", trainControl = tc, 
                          tuneGrid = rf.grid)
  
  key <- toString(n)
  modellist[[key]] <- rf_model
}

#Compare results
results <- resamples(modellist)
summary(results)
dotplot(results)
dotplot(results, metric = c("Rsquared", "RMSE"))
dotplot(results, metric = "Rsquared")
dotplot(results, metric = "RMSE")
dotplot(results, metric = "MAE")

max(modellist[['t_min']]$results$Rsquared)
min(modellist[['t_min']]$results$RMSE)


rf_model <- modellist[['t_min']]
rf_model$finalModel

importance(rf_model$finalModel)
varImpPlot(rf_model$finalModel, main = "LiD_PR_CANUND ~ S1 model", pch =17, cex = 1.5, n.var = 10, color = "mediumpurple")

#######################
##### Spearman's Correlation

colnames(v_list[['t_min']])

# Spearman correlation 1+2+3
# PR ~ def
pv <- cor.test(var_all[,'LiD_PR_CANUND_pred'], var_all[,'def_old'], na.action = "na.exclude", method="spearman")
if (pv$p.value <0.001) {pvalue <- "***"} else if (pv$p.value <0.05) {pvalue <- "**"} else if (pv$p.value <0.01) {pvalue <- "*"}else(pvalue <- "")
pv
pvalue

pv <- cor.test(var_all[,'caterpillars'], var_all[,'def_old'], na.action = "na.exclude", method="spearman")
if (pv$p.value <0.001) {pvalue <- "***"} else if (pv$p.value <0.05) {pvalue <- "**"} else if (pv$p.value <0.01) {pvalue <- "*"}else(pvalue <- "")
pv
pvalue


########################
#######    partial plot

rf_model <- modellist[['t_min']]
rf_model$finalModel$xNames


imp <- importance(rf_model$finalModel)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
set.seed(1)

op <- par(mfrow=c(2, 5))
for (i in seq_along(impvar)) {
  partialPlot(rf_model$finalModel, var_all, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]),
              ylim=c(30, 70))
}
par(op)

