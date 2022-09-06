therm_acclim_fun <- function(x, y){
  x = x - 273.15; y= y - 273.15
  xmin <- min(c(x,y)); xmax <- max(c(x,y)) 
  x1 <- seq(xmin-0.5, xmax+0.5, length.out=1000)
  
  dxy1 = data.frame(x=x, y=y)
  dxy = dxy1[order(dxy1$x),]
  txy = data.frame(x=x1, y=x1)
  
  dist <- NA
  for(i in 1:nrow(dxy)){
    dist[i] <- min(apply(txy, 1, function(x) 
      (sqrt((x[1] - dxy[i,1])^2 + (x[2] - dxy[i,2])^2)) ))
  }
  par(mfrow=c(1,2))
  plot(dxy, xlim=c(xmin, xmax), ylim=c(xmin, xmax), ty='o')
  abline(a=0, b=1)
  plot(dxy$x, dist, ty='o')
  
  acc <- c(mean(dist), cor(dxy$x, dist, method='kendall'))
  return(acc)
}


resfn <- list.files('thermo_ph/global_veg_results/', full.names = T)
head(resfn)

library(raster)
rntavg <- brick('thermo_ph/rasters/global_monTavg_range.grd')
sdtavg <- brick('thermo_ph/rasters/global_monTavg_sd.grd')


library(foreach)
library(doParallel)
registerDoParallel(4)

foreach(i = 1:length(resfn))%dopar%{

  library(zyp); library(raster); library(reshape2);
  
  a1 <- read.csv(resfn[i])
  a1$X <- NULL
  
  a1 <- a1[complete.cases(a1),]
  
  if(nrow(a1) > 0){
    
    colnames(a1) <- c('x','y','Topt','P90', 'P75', 'CC_fit', 'Rlgt', 'Tmn', 'Tmed', 'WB')
    rs <- subset(a1, Rlgt < 0.97 & CC_fit > 0.5)
    rs1 <- subset(rs, P90 <= 15)
    head(rs1)
    
    rs1$Tmn <- rs1$Tmn + 273.15
    rs1$Tmed <- rs1$Tmed + 273.15
    
    if(nrow(rs1) > 0){
    
      
    mean_Topt_dat <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$Topt, na.rm=T)), 
                            ncol =3, nrow=1)
    mean_P90_dat <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$P90, na.rm=T)),
                           ncol =3, nrow=1)
    mean_P75_dat <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$P75, na.rm=T)),
                           ncol =3, nrow=1)
    
    rownames(mean_Topt_dat) <- rownames(mean_P90_dat) <- rownames(mean_P75_dat) <- i
    
    write.table(mean_Topt_dat, 'thermo_ph/results/mean_Topt_dat.csv',
                append = T, col.names=F, sep=',')
    write.table(mean_P90_dat, 'thermo_ph/results/mean_P90_dat.csv',
                append = T, col.names=F, sep=',')
    write.table(mean_P75_dat, 'thermo_ph/results/mean_P75_dat.csv',
                append = T, col.names=F, sep=',')
    
    if(nrow(rs1)>10){ ##6 for the second rule
      topt_rst <-melt(zyp.trend.vector(rs1$Topt)[c(2, 5, 6)])
      
      p90_rst <- melt(zyp.trend.vector(rs1$P90)[c(2, 5, 6)])
      
      p75_rst <- melt(zyp.trend.vector(rs1$P75)[c(2, 5, 6)])
      
      tmp_ts <- matrix(c(rs1$x[1], rs1$y[1], 
                      topt_rst$value, p90_rst$value, p75_rst$value),
                      ncol=11, nrow=1)
      
      nTopt1 <- summary(lm(rs1$Topt ~ rs1$Tmn))
      sltopt <- ifelse(nTopt1$coefficients[8] <= 0.1, nTopt1$coefficients[2], 0)
      acc <- therm_acclim_fun(rs1$Tmn, rs1$Topt)
      
      rngTavg <- melt(extract(rntavg, rs1[1,1:2]) )
      rngTavg <- rngTavg[as.numeric(rownames(rs1)),]
      
      P90rng1 <- summary(lm(rngTavg$value ~ rs1$P90))
      P90rng <- ifelse(P90rng1$coefficients[8] <= 0.1, P90rng1$coefficients[2], 0)
      
      P75rng1 <- summary(lm(rngTavg$value ~ rs1$P75))
      P75rng <- ifelse(P75rng1$coefficients[8] <= 0.1, P75rng1$coefficients[2], 0)
      
      
      sdTavg1 <- melt(extract(sdtavg, rs1[1,1:2]) )
      sdTavg <- sdTavg1[as.numeric(rownames(rs1)),]
      
      P90sd1 <- summary(lm(sdTavg$value ~ rs1$P90))
      P90sd <- ifelse(P90sd1$coefficients[8] <= 0.1, P90sd1$coefficients[2], 0)
      
      P75sd1 <- summary(lm(sdTavg$value ~ rs1$P75))
      P75sd <- ifelse(P75sd1$coefficients[8] <= 0.1, P75sd1$coefficients[2], 0)
      
      tmp_ac <- matrix(c(rs1$x[1], rs1$y[1], 
                      sltopt, acc, P90rng, P90sd, P75rng, P75sd),
                      ncol=9, nrow=1)
      rownames(tmp_ts) <- rownames(tmp_ac) <- i
      
      write.table(tmp_ts, 'thermo_ph/results/tmp_ts.csv',
                  append = T, col.names=F, sep=',')
      write.table(tmp_ac, 'thermo_ph/results/tmp_ac.csv',
                  append = T, col.names=F, sep=',')
      
      
      print(i)
      }
    
    }
  }
 
}

