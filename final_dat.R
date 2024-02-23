library(ggplot2)
library(sp)
library(rgdal)
library(raster)
library(gtools)
library(stringr)


library(dplyr)
library(sf)




setwd('D:/DataRepository')



# load def former year
def2017_file <- read.csv(file = 'Def_09_Plot44_2017_4TimeSteps_4NCDIs_211109.csv')
def2018_file <- read.csv(file = 'Def_09_Plot44_2018_4TimeSteps_4NCDIs_211109.csv')
head(def2017_file)
def2017 <- def2017_file[,names(def2017_file) %in% c( "Plot","Defoliation")]
def2018 <- def2018_file[,names(def2018_file) %in% c( "Plot","Defoliation")]
colnames(def2017) <- c('plot', 'Def2017')
colnames(def2018) <- c('plot', 'Def2018')

# load XY
study_plots<- readOGR("Pt44_XY.shp")
head(study_plots@data)
XY <- data.frame(plot = study_plots$Name, X = study_plots$X, Y = study_plots$Y)


basic <- left_join(XY, def2017, by = 'plot')
basic <- left_join(basic, def2018, by = 'plot')
head(basic)


write.csv(basic, 'basic.csv')

##################### load S1 data ######################

fun_int <- function(band_name, span) {
  eval(parse(text=paste('band_file','<- read.csv(file = "final38/final38_', band_name,'.csv")', sep = '')))
  band_file<- band_file[,-1]
  
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
  
  cor_t 
  cor <- as.data.frame(t(cor_t))
  colnames(cor) <- date_name
  cor
  
  for (i in 2:ncol(cor_t)){
    loessMod <- loess(cor_t[,i] ~ cor_t[,1], span= span) # smoothing span
    smoothed <- predict(loessMod) 
    if (i == 2){
      smth_t <- as.data.frame(cor_t[,1])
      smth_t <- cbind(smth_t, smoothed)
    }else{
      smth_t <- cbind(smth_t, smoothed)
    }
    
  }
  
  colnames(smth_t) <- colnames(cor_t)
  smth_t
  #smth <- as.data.frame(t(smth_t))
  
  return(smth_t)
  
}


s1_VV_t <- fun_int('VV', 0.5)
s1_VH_t  <- fun_int('VH', 0.5)

s1_VV <- as.data.frame(t(s1_VV_t))
s1_VH <- as.data.frame(t(s1_VH_t))

plot <- rownames(s1_VH)
plot <- plot[-1]

write.csv(s1_VV, 'n38_bd/s1_VV_d_final.csv')
write.csv(s1_VH, 'n38_bd/s1_VH_d_final.csv')


s1_CDI <- s1_VV - s1_VH
s1_CDI[1,] <- s1_VV[1,]
# first 60 columns, Day60-119
leaf_off <- s1_CDI[, 1:60]
leaf_off

min_lo <- apply(leaf_off,1,min,na.rm = TRUE)
min_lo <- min_lo[-1]
min_lo

write.csv(min_lo, 'min_lo_final.csv')
write.csv(s1_CDI, 'n38_bd/s1_CDI_d_final.csv')


#### calculate daily indices

## Sentinel-1


# NCDI
s1_NCDI <- s1_CDI[-1,]

for (i in 1 : (ncol(s1_NCDI))){
  
  s1_NCDI[,i] <- s1_NCDI[,i] / min_lo
  
}

s1_NCDI
a <- s1_CDI[1,]
aa <- rbind(a, s1_NCDI)
s1_NCDI <- aa


## calculate Ratio
s1_ratio <- s1_VV / s1_VH
s1_ratio



## calculate RVI

for (i in 1:nrow(s1_VV)){
  b <- s1_VV[i,] + s1_VH[i,]
  a <- 4 * s1_VH[i,] / b
  if (i == 1){
    s1_RVI <- a
  }else{
    s1_RVI <- bind_rows(s1_RVI, a)
  }
}


s1_RVI



## calculate RFDI

s1_RFDI <- (s1_VV - s1_VH)/(s1_VV + s1_VH)

s1_RFDI


# save files
#write.csv(s1_CDI, "s1_CDI_d_final.csv")
write.csv(s1_NCDI, "n38_bd/s1_NCDI_d_final.csv")
write.csv(s1_ratio, "n38_bd/s1_Ratio_d_final.csv")
write.csv(s1_RVI, "n38_bd/s1_RVI_d_final.csv")
write.csv(s1_RFDI, "n38_bd/s1_RFDI_d_final.csv")




# monthly data

fun_bm <- function(band_file){
  days <- as.numeric(band_file[1,])
  ymd_days <- as.Date(days, origin = "2018-12-31")
  ymd_days
  
  month <- format(ymd_days,format="%m")
  month
  
  ## remove the days
  band_num <- band_file[-1,]

  band_t <- as.data.frame(t(band_num))
  #colnames(band_t) <- band_t[1,]
  band_t$month <- month
  
  n <- ncol(band_t)
  m <- n - 1
  
  
  
  for (i in 1:m){
    by_month <- aggregate(band_t[,i], list(band_t[,n]), mean, na.rm = TRUE)
    colnames(by_month) <- c('month', plot[i])
    if (i == 1){
      band_bm <- by_month
    }else{
      band_bm <- left_join(band_bm, by_month, by = 'month')
    }
    
  }
  
  return(band_bm)
}


s1_VH_bm <- fun_bm (s1_VH)
s1_VV_bm <- fun_bm (s1_VV)

## save VV-VH monthly composite data
write.csv(s1_VV_bm, "n38_bm/s1_VV_final.csv")
write.csv(s1_VH_bm, "n38_bm/s1_VH_final.csv")






##################### load S2 data ######################

fun_int2 <- function(band_name, span) {
  eval(parse(text=paste('band_file','<- read.csv(file = "final38/final38_', band_name,'.csv")', sep = '')))
  band_file<- band_file[,-1]
  
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
  
  cor_t <- join_t
  
  
  for (i in 2:ncol(join_t)){
    inter <- approx (join_t[,i], y = NULL, method = "linear", n = nrow(join_t), ties = mean)
    cor_t[,i] <-  inter$y
  }
  
  cor_t 
  cor <- as.data.frame(t(cor_t))
  colnames(cor) <- date_name
  cor
  
  for (i in 2:ncol(cor_t)){
    loessMod <- loess(cor_t[,i] ~ cor_t[,1], span = span) # smoothing span
    smoothed <- predict(loessMod) 
    if (i == 2){
      smth_t <- as.data.frame(cor_t[,1])
      smth_t <- cbind(smth_t, smoothed)
    }else{
      smth_t <- cbind(smth_t, smoothed)
    }
    
  }
  
  colnames(smth_t) <- colnames(cor_t)
  smth_t
  smth <- as.data.frame(t(smth_t))
  

  return(smth)
  
}


s2_RED <- fun_int2('RED',0.5)
s2_NIR <- fun_int2('NIR', 0.5)
s2_BLUE <- fun_int2('BLUE',0.5)
s2_RE2 <- fun_int2('RE2',0.5)
s2_SWIR1 <- fun_int2('SWIR1',0.5)
s2_GREEN <- fun_int2('GREEN',0.5)
s2_RE1 <- fun_int2('RE1',0.5)
s2_SWIR2 <- fun_int2('SWIR2',0.5)
s2_RE3 <- fun_int2('RE3',0.5)
s2_RE4 <- fun_int2('RE4',0.5)


write.csv(s2_RED, "n38_bd/s2_RED_d_final.csv")
write.csv(s2_NIR, "n38_bd/s2_NIR_d_final.csv")
write.csv(s2_BLUE, "n38_bd/s2_BLUE_d_final.csv")
write.csv(s2_RE2, "n38_bd/s2_RE2_d_final.csv")
write.csv(s2_SWIR1, "n38_bd/s2_SWIR1_d_final.csv")
write.csv(s2_GREEN, "n38_bd/s2_GREEN_d_final.csv")
write.csv(s2_RE1, "n38_bd/s2_RE1_d_final.csv")
write.csv(s2_SWIR2, "n38_bd/s2_SWIR2_d_final.csv")
write.csv(s2_RE3, "n38_bd/s2_RE3_d_final.csv")
write.csv(s2_RE4, "n38_bd/s2_RE4_d_final.csv")


#### calculate daily indices

## sentinel-2


s2_NDVI <- (s2_NIR -s2_RED)/(s2_NIR + s2_RED)
s2_NDWI <- (s2_NIR - s2_SWIR1)/(s2_NIR + s2_SWIR1)
s2_MSI <- s2_NIR/s2_RED
s2_PSRI <- (s2_RED - s2_BLUE)/s2_RE2
s2_RERVI <- s2_NIR/s2_RE2
s2_RENDVI <- (s2_NIR -s2_RE2)/(s2_NIR + s2_RE2)
s2_REEVI2 <- 2.5 * (s2_NIR -s2_RE2)/(s2_NIR + 2.4 * s2_RE2 + 1)


write.csv(s2_NDVI, "n38_bd/s2_NDVI_d_final.csv")
write.csv(s2_NDWI, "n38_bd/s2_NDWI_d_final.csv")
write.csv(s2_MSI, "n38_bd/s2_MSI_d_final.csv")
write.csv(s2_PSRI, "n38_bd/s2_PSRI_d_final.csv")
write.csv(s2_RERVI, "n38_bd/s2_RERVI_d_final.csv")
write.csv(s2_RENDVI, "n38_bd/s2_RENDVI_d_final.csv")
write.csv(s2_REEVI2, "n38_bd/s2_REEVI2_d_final.csv")


### monthly data

s2_RED_bm <- fun_bm (s2_RED)
s2_NIR_bm <- fun_bm (s2_NIR)
s2_BLUE_bm <- fun_bm (s2_BLUE)
s2_RE2_bm <- fun_bm (s2_RE2)
s2_SWIR1_bm <- fun_bm (s2_SWIR1)
s2_GREEN_bm <- fun_bm (s2_GREEN)
s2_RE1_bm <- fun_bm (s2_RE1)
s2_SWIR2_bm <- fun_bm (s2_SWIR2)
s2_RE3_bm <- fun_bm (s2_RE3)
s2_RE4_bm <- fun_bm (s2_RE4)

write.csv(s2_RED_bm, "n38_bm/s2_RED_final.csv")
write.csv(s2_NIR_bm, "n38_bm/s2_NIR_final.csv")
write.csv(s2_BLUE_bm, "n38_bm/s2_BLUE_final.csv")
write.csv(s2_RE2_bm, "n38_bm/s2_RE2_final.csv")
write.csv(s2_SWIR1_bm, "n38_bm/s2_SWIR1_final.csv")
write.csv(s2_GREEN_bm, "n38_bm/s2_GREEN_final.csv")
write.csv(s2_RE1_bm, "n38_bm/s2_RE1_final.csv")
write.csv(s2_SWIR2_bm, "n38_bm/s2_SWIR2_final.csv")
write.csv(s2_RE3_bm, "n38_bm/s2_RE3_final.csv")
write.csv(s2_RE4_bm, "n38_bm/s2_RE4_final.csv")


s2_NDVI <- (s2_NIR -s2_RED)/(s2_NIR + s2_RED)
s2_NDWI <- (s2_NIR - s2_SWIR1)/(s2_NIR + s2_SWIR1)
s2_MSI <- s2_NIR/s2_RED
s2_PSRI <- (s2_RED - s2_BLUE)/s2_RE2
s2_RERVI <- s2_NIR/s2_RE2
s2_RENDVI <- (s2_NIR -s2_RE2)/(s2_NIR + s2_RE2)
s2_REEVI2 <- 2.5 * (s2_NIR -s2_RE2)/(s2_NIR + 2.4 * s2_RE2 + 1)

write.csv(s2_NDVI, "n38_bd/s2_NDVI_d_final.csv")
write.csv(s2_NDWI, "n38_bd/s2_NDWI_d_final.csv")
write.csv(s2_MSI, "n38_bd/s2_MSI_d_final.csv")
write.csv(s2_PSRI, "n38_bd/s2_PSRI_d_final.csv")
write.csv(s2_RERVI, "n38_bd/s2_RERVI_d_final.csv")
write.csv(s2_RENDVI, "n38_bd/s2_RENDVI_d_final.csv")
write.csv(s2_REEVI2, "n38_bd/s2_REEVI2_d_final.csv")

#####################################
######### calculate indices #########
#####################################


#### load files

min_lo <- read.csv('min_lo_final.csv', header = T, row.names = 1)

month_name <- c('Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')

fun_prep <- function(satellite, band){
  eval(parse(text=paste('file','<- read.csv(file = "n38_bm/', satellite, '_', band, '_final.csv")', sep = '')))
  file <- file[,-1]
  t <- as.data.frame(t(file))
  colnames(t) <- paste(month_name, '_', band, sep='')
  t <- t[-1,]
  return(t)
}

s1_VV <- fun_prep('s1','VV')
s1_VH <- fun_prep('s1','VH')

s2_RED <- fun_prep('s2','RED')
s2_NIR <- fun_prep('s2','NIR')
s2_BLUE <- fun_prep('s2','BLUE')
s2_RE2 <- fun_prep('s2','RE2')
s2_SWIR1 <- fun_prep('s2','SWIR1')


#### calculate indices

## Sentinel-1
# CDI
s1_CDI <- s1_VV - s1_VH
colnames(s1_CDI) <- paste(month_name, '_CDI', sep='')

# NCDI
s1_NCDI <- s1_CDI

for (i in 1 : ncol(s1_NCDI)){
  
  s1_NCDI[,i] <- s1_NCDI[,i] / min_lo
  
}

s1_NCDI

colnames(s1_NCDI) <- paste(month_name, '_NCDI', sep = '')

## calculate Ratio
s1_ratio <- s1_VV / s1_VH
s1_ratio
colnames(s1_ratio) <- paste(month_name, '_Ratio', sep = '')


## calculate RVI

for (i in 1:nrow(s1_VV)){
  b <- s1_VV[i,] + s1_VH[i,]
  a <- 4 * s1_VH[i,] / b
  if (i == 1){
    s1_RVI <- a
  }else{
    s1_RVI <- bind_rows(s1_RVI, a)
  }
}


s1_RVI
colnames(s1_RVI) <- paste(month_name, '_RVI', sep = '')


## calculate RFDI

s1_RFDI <- (s1_VV - s1_VH)/(s1_VV + s1_VH)

s1_RFDI
colnames(s1_RFDI) <- paste(month_name, '_RFDI', sep = '')

# save files
write.csv(s1_CDI, "n38_bm/s1_CDI_final.csv")
write.csv(s1_NCDI, "n38_bm/s1_NCDI_final.csv")
write.csv(s1_ratio, "n38_bm/s1_Ratio_final.csv")
write.csv(s1_RVI, "n38_bm/s1_RVI_final.csv")
write.csv(s1_RFDI, "n38_bm/s1_RFDI_final.csv")

## sentinel-2


s2_NDVI <- (s2_NIR -s2_RED)/(s2_NIR + s2_RED)
s2_NDWI <- (s2_NIR - s2_SWIR1)/(s2_NIR + s2_SWIR1)
s2_MSI <- s2_NIR/s2_RED
s2_PSRI <- (s2_RED - s2_BLUE)/s2_RE2
s2_RERVI <- s2_NIR/s2_RE2
s2_RENDVI <- (s2_NIR -s2_RE2)/(s2_NIR + s2_RE2)
s2_REEVI2 <- 2.5 * (s2_NIR -s2_RE2)/(s2_NIR + 2.4 * s2_RE2 + 1)


colnames(s2_NDVI) <- paste(month_name, '_NDVI', sep = '')
colnames(s2_NDWI) <- paste(month_name, '_NDWI', sep = '')
colnames(s2_MSI) <- paste(month_name, '_MSI', sep = '')
colnames(s2_PSRI) <- paste(month_name, '_PSRI', sep = '')
colnames(s2_RERVI) <- paste(month_name, '_RERVI', sep = '')
colnames(s2_RENDVI) <- paste(month_name, '_RENDVI', sep = '')
colnames(s2_REEVI2) <- paste(month_name, '_REEVI2', sep = '')

write.csv(s2_NDVI, "n38_bm/s2_NDVI_final.csv")
write.csv(s2_NDWI, "n38_bm/s2_NDWI_final.csv")
write.csv(s2_MSI, "n38_bm/s2_MSI_final.csv")
write.csv(s2_PSRI, "n38_bm/s2_PSRI_final.csv")
write.csv(s2_RERVI, "n38_bm/s2_RERVI_final.csv")
write.csv(s2_RENDVI, "n38_bm/s2_RENDVI_final.csv")
write.csv(s2_REEVI2, "n38_bm/s2_REEVI2_final.csv")


################# Defoliation Index #################

## 0.26(delta CDI)+0.52(dleta RERVI)
## delta CDI

CDI_d <- read.csv('n38_bd/S1_CDI_d_final.csv')
RERVI_d <- read.csv('n38_bd/S2_RERVI_d_final.csv')


days <- CDI_d[1,]
plot <- CDI_d[,1]
days <- days[-1]
plot <- plot[-1]

CDI_d <- CDI_d[-1,]
CDI_d <- CDI_d[,-1]

RERVI_d <- RERVI_d[-1,]
RERVI_d <- RERVI_d[,-1]

defoli <- data.frame(plot = plot)
defoli$A <- apply(CDI_d[,49:60],1,mean,na.rm = TRUE)
defoli$B <- apply(CDI_d[,73:108],1,max,na.rm = TRUE)
defoli$C <- apply(CDI_d[,109:180],1,min,na.rm = TRUE)
defoli$D <- apply(CDI_d[,193:211],1,max,na.rm = TRUE)

defoli$def_old <- defoli$A - defoli$C

defoli$def_ac <- defoli$A - defoli$C
defoli$def_bc <- defoli$B - defoli$C
defoli$def_dc <- defoli$D - defoli$C
defoli$def_db <- defoli$D - defoli$B

defoli$RA <- apply(RERVI_d[,49:60],1,mean,na.rm = TRUE)
defoli$RB <- apply(RERVI_d[,73:108],1,max,na.rm = TRUE)
defoli$RC <- apply(RERVI_d[,109:180],1,min,na.rm = TRUE)
defoli$RD <- apply(RERVI_d[,193:211],1,max,na.rm = TRUE)

defoli$def_new <- 0.26 * (defoli$A - defoli$C) + 0.52 * (defoli$RA - defoli$RC)
write.csv(defoli, "def_final.csv")
