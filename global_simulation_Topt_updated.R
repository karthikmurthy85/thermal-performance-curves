#setwd('/media/karthik/karthik_cluster/')
source('thermo_ph/codes/simulation_functions_simplified.R')
library(raster); library(rgdal);
ndvi_veg <- brick('thermo_ph/rasters/global_veg9km.grd')
ndvi_veg_df <- data.frame(rasterToPoints(ndvi_veg))
ndvi_veg_df$layer <- NULL

ndvi_mon <- brick('thermo_ph/rasters/global_NDVI_GIMMS_9km.grd')
temp_ras <- brick('thermo_ph/rasters/global_monthly_Tavg.grd')
wb_ras <- brick('thermo_ph/rasters/water_balance_monthly.grd')
ids <- seq(6, 414, 12) 
ll <- 10; yrs <- 1:34
wb_val1 <- data.frame(Var1 = rep(1,6), Var2 = rep("layer.1",6), value = rep(NA,6))

##To run all the missing files
fils <- list.files('thermo_ph/global_veg_results/')
head(fils)

library(stringr)
fils1 <- as.numeric(str_sub(fils, 5, 12))
nn <- range(fils1)
afils <- nn[1] : nn[2] 
fls1 <- afils[-fils1]



library(foreach)
library(doParallel)
registerDoParallel(8*8)
#nrow(ndvi_veg_df)
#foreach(j = 1200001:nrow(ndvi_veg_df))%dopar%{
foreach(j = fls1)%dopar%{
#for(j in fls1){
  library(raster); library(reshape2); library(stringr)
  ndvi_val <- melt(extract(ndvi_mon, ndvi_veg_df[j,]))
  temp_val <- melt(extract(temp_ras, ndvi_veg_df[j,]))
  wb_val <-  rbind(wb_val1, melt(extract(wb_ras, ndvi_veg_df[j,])))
  
  topt_df <- data.frame( TEMP = temp_val$value, NDVI = ndvi_val$value, WB = wb_val$value)
  topt_rs <- matrix(NA, ncol = 10, nrow = (length(yrs)-ll))
  for(i in 1:(length(yrs)-ll)){
    topt_df1 <- rbind(topt_df[(ids[i]-5):ids[i],],
                      topt_df[(ids[i]+1):ids[ll+i+1],])
    topt_df1$NDVI[topt_df1$NDVI < 0] <- NA
    topt_df1$NDVI <-  topt_df1$NDVI/10000
    topt_df1$TEMP <- topt_df1$TEMP - 273.15
    wb_mean <- mean(topt_df1$WB, na.rm=T)
    topt_df1$WB <- NULL
    topt_df2 <- topt_df1[complete.cases(topt_df1),]
    if(nrow(topt_df2) < 55){
      next
    }
    plot(topt_df2)
    
    npp_qn <- matrix(NA, ncol=4, nrow=1)
    kss <- seq(1,4,1)
    for(k in kss){
      npp_qn <- rbind(npp_qn, quan_fun(k, topt_df2))
    }
    
    npp_qn <- data.frame(npp_qn[complete.cases(npp_qn),])
    npp_qn <- npp_qn[order(npp_qn$X1),]
    npp_q <- subset(npp_qn, X1 > min(topt_df2$TEMP))
    
    np1 <- npp_q[complete.cases(npp_q),]
    np2 <- melt(with(np1, tapply(X3, X1, median)))
    colnames(np2) <- c('tmp', 'npp')
    points(np2, ty='l', col = 'blue')
    
    rs1 <- cc_l1 <- cc_l2 <- md1 <- md2 <- lrsq <- NA
    try(rs1 <- schoolfield_fun(np2))
    try(cc_l1 <- logistic_inv(np2$npp))
    try(cc_l2 <- logistic(np2$npp))
    try(md1 <- summary(lm(cc_l1[-1] ~ rs1[[2]][,1])))
    try(md2 <- summary(lm(cc_l2[-1] ~ rs1[[2]][,1])))
    try(md1$r.squared[md1$r.squared %in% NaN] <- 0)
    try(md2$r.squared[md2$r.squared %in% NaN] <- 0)
    
    try(lrsq <- ifelse(md1$r.squared > md2$r.squared, md1$r.squared, md2$r.squared))
    
    try(topt_rs[i,] <- c(ndvi_veg_df[j,1], ndvi_veg_df[j,2], 
                     rs1[[1]], lrsq, mean(topt_df2$TEMP), median(topt_df2$TEMP), wb_mean))
    print(i)
  }
  topt_rs <- round(topt_rs, 3)
  j1 <- format(j, scientific = F)
  fid <- str_pad(j1, pad='0', width=8)
  fnam <- paste0('thermo_ph/global_veg_results/pix_',
                 fid, '.csv')
  write.csv(topt_rs,fnam)
  print(paste("pixel id---->", j1, "finished"))
}
