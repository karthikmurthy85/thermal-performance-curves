#setwd('/media/karthik/karthik_cluster/')
source('thermo_ph/codes/simulation_functions.R')
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


library(foreach)
library(doParallel)
registerDoParallel(64)
#nrow(ndvi_veg_df)
foreach(j = 109400:200000)%dopar%{
  library(raster); library(reshape2); library(stringr)
  ndvi_val <- melt(extract(ndvi_mon, ndvi_veg_df[j,]))
  temp_val <- melt(extract(temp_ras, ndvi_veg_df[j,]))
  wb_val <-  rbind(wb_val1, melt(extract(wb_ras, ndvi_veg_df[j,])))
  
  topt_df <- data.frame( TEMP = temp_val$value, NDVI = ndvi_val$value, WB = wb_val$value)
  topt_rs <- matrix(NA, ncol = 9, nrow = (length(yrs)-ll))
  for(i in 1:(length(yrs)-ll)){
    topt_df1 <- rbind(topt_df[(ids[i]-5):ids[i],],
                      topt_df[(ids[i]+1):ids[ll+i+1],])
    topt_df1$NDVI[topt_df1$NDVI < 0] <- NA
    topt_df1$NDVI <-  topt_df1$NDVI/10000
    topt_df1$TEMP <- topt_df1$TEMP - 273.15
    wb_mean <- mean(topt_df1$WB, na.rm=T)
    topt_df1$WB <- NULL
    topt_df2 <- topt_df1[complete.cases(topt_df1),]
    if(nrow(topt_df2) < 75){
      next
    }
    #plot(topt_df2)
    
    
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
    #points(np2, ty='l', col = 'blue')
    
    minval <- min(np2$npp)
    np2$npp <- np2$npp - minval
    
    if(min(np2$tmp) > 0){
      #ky <- 1
      rs1 <- schoolfield_fun(np2, minval)
      cc_tropical <- logistic_inv(np2$npp)
      md <- summary(lm(cc_tropical[-1] ~ rs1[[2]][,1]))
      topt_rs[i,] <- c(ndvi_veg_df[j,1], ndvi_veg_df[j,2], 
                       rs1[[1]], md$r.squared, mean(topt_df2$TEMP), wb_mean)
    }
    
    else{
      nnp <- np2
      nnp$tmp <- np2$tmp + abs(min(np2$tmp)) + 5
      rs2 <- list(); rs2[[1]] <- rep(NA, 5)
      try(rs2 <- schoolfield_fun(nnp, minval))
      if(length(rs2[[1]])> 4){
        res_1 <- c(NA,NA,NA,NA,NA,NA,NA,NA, NA)
      }
      else{
        rs2[[1]][1] <- rs2[[1]][1] - (abs(min(np2$tmp)) + 5)
        cc_tropical <- logistic_inv(nnp$npp)
        md1 <- summary(lm(cc_tropical[-1] ~ rs2[[2]][,1]))
        res_1 <- c(ndvi_veg_df[j,1], ndvi_veg_df[j,2], 
                   rs2[[1]], md1$r.squared, mean(topt_df2$TEMP), wb_mean)
      }
      
      np3 <- np2
      np3$tmp <- (np2$tmp*-1)
      tmp_add <- abs(min(np3$tmp)) + 5
      np3$tmp <- np3$tmp + tmp_add
      np4 <- np3[order(np3$tmp),]
      rs1 <- list(); rs1[[1]] <- rep(NA, 5)
      try(rs1 <- schoolfield_fun(np3, minval))
      if(length(rs1[[1]])> 4){
        res_2 <- c(NA,NA,NA,NA,NA,NA,NA,NA, NA)
      }
      else{
        rs1[[1]][1] <- ((rs1[[1]][1]) - tmp_add)*-1
        logrs <- logistic_inv(np4$npp)[-1]
        pred <- rs1[[2]]
        pred <- pred[order(pred$temp),]
        md <- summary(lm(logrs ~ pred$per))
        res_2 <- c(ndvi_veg_df[j,1], ndvi_veg_df[j,2], 
                   rs1[[1]], md$r.squared, mean(topt_df2$TEMP), wb_mean)
      }
      res_2[res_2 %in% NA| res_2 == -Inf | res_2 == Inf] <- 10000
      res_1[res_1 %in% NA| res_1 == -Inf | res_1 == Inf] <- 10000
      ifelse(sum(res_1[4:5]) < sum(res_2[4:5]), 
             topt_rs[i,] <- res_1, topt_rs[i,] <- res_2)
    }
    #print(i)
  }
  j1 <- format(j, scientific = F)
  fid <- str_pad(j1, pad='0', width=8)
  fnam <- paste0('thermo_ph/global_veg_results/pix_',
                 fid, '.csv')
  write.csv(topt_rs,fnam)
  print(paste("pixel id---->", j1, "finished"))
}
