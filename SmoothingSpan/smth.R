library(shiny)
library(stringr)
library(tidyr)
library(dplyr)

smth_df <- function(file_name, band_name, span1, span2, span3){
  band_file <- read.csv(file = file_name, header = T, row.names = 1)
  
  band_days <- as.numeric(str_extract_all(colnames(band_file),"[0-9]+[0-9]"))
  band_days <- band_days[-1]
  
  int_days <- seq(from = 60, to = 270)
  date_name <- paste('Day', int_days, sep = '')
  
  eval(parse(text=paste('colnames(band_file)','<- paste("', band_name, '_", colnames(band_file), sep = "")', sep = '')))
  eval(parse(text=paste('colnames(band_file)[colnames(band_file) == "', band_name, '_id"] <- "id"', sep = '')))
  
  data.matrix<-matrix(nrow = 1, ncol = length(int_days))
  int <- data.frame(data.matrix)
  
  colnames(int) <- date_name
  int
  
  ## interpolate
  
  band_file_t <- as.data.frame(t(band_file))
  colnames(band_file_t) <-band_file_t[1,]
  band_file_t <- band_file_t[-1,]
  band_file_t$days <- band_days
  
  int_t <- as.data.frame(t(int))
  int_t$days <- int_days
  join_t <- left_join(int_t, band_file_t , by="days")
  join_t <- join_t[,-1]
  join_t
  
  join_t <- join_t[,colSums(is.na(join_t)) < nrow(join_t) - 1]
  
  cor_t <- join_t
  
  
  for (i in 2:ncol(join_t)){
    inter <- approx (join_t[,i], y = NULL, method = "linear", n = nrow(join_t), ties = mean)
    cor_t[,i] <-  inter$y
  }
  
  rownames(cor_t) <- date_name
  orig <- cor_t
  
  # smoothing 1
  
  for (i in 2:ncol(cor_t)){
    loessMod <- loess(cor_t[,i] ~ cor_t[,1], span = span1) # smoothing span
    smoothed <- predict(loessMod) 
    if (i == 2){
      smth_t <- as.data.frame(cor_t[,1])
      smth_t <- cbind(smth_t, smoothed)
    }else{
      smth_t <- cbind(smth_t, smoothed)
    }
  }
  
  colnames(smth_t) <- colnames(cor_t)
  s1 <- smth_t
  
  # smoothing 2
  
  for (i in 2:ncol(cor_t)){
    loessMod <- loess(cor_t[,i] ~ cor_t[,1], span = span2) # smoothing span
    smoothed <- predict(loessMod) 
    if (i == 2){
      smth_t <- as.data.frame(cor_t[,1])
      smth_t <- cbind(smth_t, smoothed)
    }else{
      smth_t <- cbind(smth_t, smoothed)
    }
  }
  
  colnames(smth_t) <- colnames(cor_t)
  s2 <- smth_t
  
  # smoothing 3
  
  for (i in 2:ncol(cor_t)){
    loessMod <- loess(cor_t[,i] ~ cor_t[,1], span = span3) # smoothing span
    smoothed <- predict(loessMod) 
    if (i == 2){
      smth_t <- as.data.frame(cor_t[,1])
      smth_t <- cbind(smth_t, smoothed)
    }else{
      smth_t <- cbind(smth_t, smoothed)
    }
  }
  
  colnames(smth_t) <- colnames(cor_t)
  s3 <- smth_t
  
  df <- data.frame(days = orig[,1], origin = orig[,2], s1 = s1[,2], s2 = s2[,2], s3 = s3[,2])
  col_s1 <- paste('s', span1, sep = '')
  col_s2 <- paste('s', span2, sep = '')
  col_s3 <- paste('s', span3, sep = '')
  
  colnames(df) <- c('DOY','origin',col_s1,col_s2,col_s3)
  eval(parse(text = paste('df2 <- df %>%
    gather(key = "span", value = ', band_name, ',-DOY)',sep = '')))
  return(df2)
}