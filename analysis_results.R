sen_slp_segs <- function(y1, bcp){
  bcp1 <- round(bcp) ##Get the integer value of the breakpoint
  yt1 <- y1[1:bcp1] ##Get the data for the first segment
  yt2 <- y1[(bcp1+1):length(y1)] ##Get the data for the second segment
  
  mod1 <-  zyp.trend.vector(yt1) ## senslope model for the first segement
  mod2 <- zyp.trend.vector(yt2) ## senslope model for the second segement
  
  return(c(mod1[2], mod1[6],mod2[2], mod2[6])) ##1st segment slope, 1st segment pvalue, 2nd segment slope, 2nd segment pvalue,
  
}

seg_regr_calc <- function(y1){
  if(length(y1[is.na(y1)]) > 0| var(y1[!is.na(y1)]) == 0){
    return(c(NA,NA,NA,NA,NA,NA))
  }
  else{
    t <- 1:length(y1)
    mod <- lm(y1~t) #Linear model
    s1 <- summary(mod)
    
    library(segmented)
    for(i in 1:20){
      smod1 <- NA
      try(smod1 <- segmented.lm(mod, seg.Z = ~t,
                                control = seg.control(display = FALSE, n.boot=5), K=1))
      if(class(smod1)[1] %in% "segmented"){break}
      #print(i)
    }
    
    if(class(smod1)[1] %in% 'lm'){
      return(c(NA,NA,NA,NA,NA,NA))
    }
    
    else{
      s2 <- summary(smod1)
      rel_rsq <- (s2$r.squared - s1$r.squared)
      sen_res <- sen_slp_segs(y1, s2$psi[2])
      return(c(s2$psi[2], rel_rsq, sen_res)) ##Breakpoint, relative Rsq
    }
    
  }
}

therm_acclim_fun <- function(x, y){
  x = x - 273.15; y= y - 273.15
  xmin <- min(c(x,y)); xmax <- max(c(x,y)) 
  x1 <- seq(xmin-0.5, xmax+0.5, length.out=10000)
  
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

# therm_acclim_fun <- function(x, y){
#   
#   library(sf); library(geosphere); library(rgdal)
#   x = x - 273.15; y= y - 273.15
#   xmin <- min(c(x,y)); xmax <- max(c(x,y)) 
#   x1 <- seq(xmin-0.5, xmax+0.5, length.out=10000)
#   dxy1 = data.frame(x=x, y=y)
#   dxy = dxy1[order(dxy1$x),]
#   txy = data.frame(x=x1, y=x1)
#   
#   points.sf <- SpatialPointsDataFrame(dxy, dxy, match.ID = TRUE) 
# 
#   tl <- Line(txy)
#   tl1 <- Lines(list(tl), ID="L1")
#   tl2 <- SpatialLines(list(tl1))
#   df <- data.frame(len = 1); rownames(df) <- 'L1'
#   line.sf <- SpatialLinesDataFrame(tl2, df, match.ID = TRUE)
# 
#   dist <- dist2Line(p = points.sf@coords, line = line.sf)
#   
#   par(mfrow=c(1,2))
#   plot(dxy, xlim=c(xmin, xmax), ylim=c(xmin, xmax), ty='o')
#   abline(a=0, b=1)
#   plot(dist[,1], ty='o')
#   return(zyp.trend.vector(dist[,1])[c(2,5,6,11)])
# }

setwd('/media/karthik/karthik_cluster/')
resfn <- list.files('thermo_ph/global_veg_results/', full.names = T)
head(resfn)

library(raster)
rntavg <- brick('thermo_ph/rasters/global_monTavg_range.grd')
sdtavg <- brick('thermo_ph/rasters/global_monTavg_sd.grd')

mean_Topt_dat <- matrix(NA, nrow = length(resfn), ncol =3)
mean_P90_dat <- matrix(NA, nrow = length(resfn), ncol =3)
mean_P75_dat <- matrix(NA, nrow = length(resfn), ncol =3)
# tmp_ts <- matrix(NA, nrow = length(resfn), ncol = 11)
# tmp_ac <- matrix(NA, nrow = length(resfn), ncol = 9)


library(zyp); library(raster); library(reshape2)


library(foreach)
library(doParallel)
registerDoParallel(4)

for(i in 1:length(resfn)){
  a1 <- read.csv(resfn[i])
  a1$X <- NULL
  
  a1 <- a1[complete.cases(a1),]
  
  if(nrow(a1) > 0){
    
    colnames(a1) <- c('x','y','Topt','P90', 'P75', 'CC_fit', 'Rlgt', 'Tmn', 'Tmed', 'WB')
    rs <- subset(a1, Rlgt < 0.96 & CC_fit > 0.55)
    rs1 <- subset(rs, P90 <=25)
    head(rs1)
    
    rs1$Tmn <- rs1$Tmn + 273.15
    rs1$Tmed <- rs1$Tmed + 273.15
    
    if(nrow(rs1) > 0){
    
      
    mean_Topt_dat[i,] <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$Topt, na.rm=T)), 
                            ncol =3, nrow=1)
    mean_P90_dat[i,] <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$P90, na.rm=T)),
                           ncol =3, nrow=1)
    mean_P75_dat[i,] <- matrix(c(rs1$x[1], rs1$y[1], mean(rs1$P75, na.rm=T)),
                           ncol =3, nrow=1)
    
    # rownames(mean_Topt_dat) <- rownames(mean_P90_dat) <- rownames(mean_P75_dat) <- i
    # 
    # write.table(mean_Topt_dat, 'thermo_ph/results/mean_Topt_dat.csv',
    #             append = T, col.names=F, sep=',')
    # write.table(mean_P90_dat, 'thermo_ph/results/mean_P90_dat.csv',
    #             append = T, col.names=F, sep=',')
    # write.table(mean_P75_dat, 'thermo_ph/results/mean_P75_dat.csv',
    #             append = T, col.names=F, sep=',')
    # 
    # if(nrow(rs1)>10){
    #   topt_rst <-melt(zyp.trend.vector(rs1$Topt)[c(2, 5, 6)])
    #   
    #   p90_rst <- melt(zyp.trend.vector(rs1$P90)[c(2, 5, 6)])
    #   
    #   p75_rst <- melt(zyp.trend.vector(rs1$P75)[c(2, 5, 6)])
    #   
    #   tmp_ts <- matrix(c(rs1$x[1], rs1$y[1], 
    #                   topt_rst$value, p90_rst$value, p75_rst$value),
    #                   ncol=11, nrow=1)
    #   
    #   nTopt1 <- summary(lm(rs1$Topt ~ rs1$Tmn))
    #   sltopt <- ifelse(nTopt1$coefficients[8] <= 0.1, nTopt1$coefficients[2], 0)
    #   acc <- therm_acclim_fun(rs1$Tmed, rs1$Topt)
    #   
    #   rngTavg <- melt(extract(rntavg, rs1[1,1:2]) )
    #   rngTavg <- rngTavg[as.numeric(rownames(rs1)),]
    #   
    #   P90rng1 <- summary(lm(rngTavg$value ~ rs1$P90))
    #   P90rng <- ifelse(P90rng1$coefficients[8] <= 0.1, P90rng1$coefficients[2], 0)
    #   
    #   P75rng1 <- summary(lm(rngTavg$value ~ rs1$P75))
    #   P75rng <- ifelse(P75rng1$coefficients[8] <= 0.1, P75rng1$coefficients[2], 0)
    #   
    #   
    #   sdTavg1 <- melt(extract(sdtavg, rs1[1,1:2]) )
    #   sdTavg <- sdTavg1[as.numeric(rownames(rs1)),]
    #   
    #   P90sd1 <- summary(lm(sdTavg$value ~ rs1$P90))
    #   P90sd <- ifelse(P90sd1$coefficients[8] <= 0.1, P90sd1$coefficients[2], 0)
    #   
    #   P75sd1 <- summary(lm(sdTavg$value ~ rs1$P75))
    #   P75sd <- ifelse(P75sd1$coefficients[8] <= 0.1, P75sd1$coefficients[2], 0)
    #   
    #   tmp_ac <- matrix(c(rs1$x[1], rs1$y[1], 
    #                   sltopt, acc, P90rng, P90sd, P75rng, P75sd),
    #                   ncol=9, nrow=1)
    #   rownames(tmp_ts) <- rownames(tmp_ac) <- i
    #   
    #   write.table(tmp_ts, 'thermo_ph/results/tmp_ts.csv',
    #               append = T, col.names=F, sep=',')
    #   write.table(tmp_ac, 'thermo_ph/results/tmp_ac.csv',
    #               append = T, col.names=F, sep=',')
    #   
    #   
    #   print(i)
    #   }
    print(i)
    }
  }
 
}

# write.csv(mean_Topt_dat[1:20000, ], 'thermo_ph/results/mean_Topt_dat.csv')
# write.csv(mean_P90_dat[1:20000, ], 'thermo_ph/results/mean_P90_dat.csv')
# write.csv(mean_P75_dat[1:20000, ], 'thermo_ph/results/mean_P75_dat.csv')
# write.csv(tmp_ts[1:20000, ], 'thermo_ph/results/tmp_ts.csv')
# write.csv(tmp_ac[1:20000, ], 'thermo_ph/results/tmp_ac.csv')




write.csv(mean_Topt_dat, 'thermo_ph/results/mean_topt.csv')
write.csv(mean_P90_dat, 'thermo_ph/results/mean_P90.csv')
write.csv(mean_P75_dat, 'thermo_ph/results/mean_P75.csv')




mn_topt_dt <- read.csv('thermo_ph/results/mean_topt.csv')
head(mn_topt_dt)
plot(mn_topt_dt$V2 ~ mn_topt_dt$V3, ty='l')

mn_p90_dt <- read.csv('thermo_ph/results/mean_P90.csv')
mn_p90_dt <- subset(mn_p90_dt, V3 < 25)
head(mn_p90_dt)
plot(mn_p90_dt$V2 ~ mn_p90_dt$V3, ty='l')

##plot latitudinal variation in Topt
lat_rng <- round(range(mn_topt_dt$V2, na.rm=T),2)
intv <- seq(lat_rng[1], lat_rng[2], 0.25)

library(reshape2)
latrdf <- matrix(NA, ncol = 3, nrow = (length(intv)-1))
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(mn_topt_dt, V2 >= intv[i] & V2 < intv[i+1])
  midlatr <- mean(c(intv[i], intv[i+1]))
  topts <- melt(quantile(latrdf1$V3, probs=c(0.5, 0.9)))
  latrdf[i,] <- c(midlatr, topts$value)
  print(i)
}

library(ggplot2)
latrdf <- data.frame(latrdf)
colnames(latrdf) <- c('latitude', 'T050', 'T090')

ttl <- expression('Latitudinal variation of'~T[opt])
g11 <- ggplot(latrdf, aes(y = T050 - 273.15, x = latitude))+
  geom_point(size=0.5) + geom_smooth(method='gam', span=0.1)+
  theme_bw() +
  xlab('latitude') + ylab('Topt') +
  ggtitle(ttl) + 
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )
g11

##Plot latitudinal variation in P75
mn_p90_dt <- data.frame(mean_P75_dat)
colnames(mn_p90_dt) <- c("V1", "V2", "V3")
library(reshape2)
latrdf <- matrix(NA, ncol = 3, nrow = (length(intv)-1))
mn_p90_dt$V3[mn_p90_dt$V3 > 10] <- NA
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(mn_p90_dt, V2 >= intv[i] & V2 < intv[i+1])
  midlatr <- mean(c(intv[i], intv[i+1]))
  latrdf1[] <- Map(function(x) replace(x, is.infinite(x), NA), latrdf1)
  P75s <- melt(quantile(latrdf1$V3, probs=c(0.5, 0.9), na.rm=T))
  latrdf[i,] <- c(midlatr, P75s$value)
  print(i)
}

library(ggplot2)
latrdf <- data.frame(latrdf)
colnames(latrdf) <- c('latitude', 'T050', 'T090')

ttl <- expression('Latitudinal variation of'~T[0.9])
g22 <- ggplot(latrdf, aes(y = T050, x = latitude))+
  geom_point(size=0.5) + geom_smooth(method='loess', span=0.08)+
  theme_bw() +
  xlab('latitude') + ylab('High performance range (90%)') +
  ggtitle(ttl) + 
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )
g22
gridExtra::grid.arrange(g11, g22, ncol=2)


##Mapping Topt##
ttl <- expression("Optimum temperature"~"("~T[opt]~")")
mn_topt_dt <- data.frame(mean_Topt_dat)
mn_topt_dt$X3[mn_topt_dt$X3 < 260] <- NA
hist(mn_topt_dt$X3-273.15, breaks=40)
cols = c("grey70", "green1", "green2", "green3",
         "royalblue1","royalblue2", "royalblue4",
         "yellow1", "yellow2", 
         "tomato1", "tomato2", "tomato3", "red4")
contnts <- readOGR('/media/karthik/ADATA HD720/global_veg/continent_shapefile/continent shapefile/','continent')
contnts <- subset(contnts, ! CONTINENT %in% 'Oceania' )
contnt_shp <- fortify(contnts)
g1 <- ggplot() + 
  geom_polygon(data = contnt_shp,
               aes(x = long, y = lat, group = group),
               fill = 'white', colour= 'black', size=0.1)+
  geom_tile(data = mn_topt_dt, aes(x=X1, y=X2, fill = (X3-273.15))) + 
  coord_equal() +
  scale_fill_gradientn("Temperature (°C)", colours = cols) +
  theme_bw()+xlab("")+ylab("")+ylim(-55, 85)+
  ggtitle(ttl)+
  theme(
    axis.text = element_text(size =15, color = 'black'),
    axis.ticks = element_blank(),
    panel.grid = element_blank(), 
    panel.background = element_rect(colour= 'black', fill = 'lightblue1'),
    legend.text = element_text(size = 16, color = 'black'),
    legend.title =  element_text(size =16, color = 'black'),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(size=20, colour='black', hjust=0.5)
  )##multi_returnrate 800 500

##plot latitudinal variation in Topt
lat_rng <- round(range(mn_topt_dt$X2, na.rm=T),2)
intv <- seq(lat_rng[1], lat_rng[2], 0.25)

library(reshape2)
latrdf <- matrix(NA, ncol = 3, nrow = (length(intv)-1))
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(mn_topt_dt, X2 >= intv[i] & X2 < intv[i+1])
  midlatr <- mean(c(intv[i], intv[i+1]))
  topts <- melt(quantile(latrdf1$X3, probs=c(0.5, 0.9), na.rm=T))
  latrdf[i,] <- c(midlatr, topts$value)
  print(i)
}

library(ggplot2)
latrdf <- data.frame(latrdf)
colnames(latrdf) <- c('latitude', 'T050', 'T090')

ttl <- expression('Latitudinal variation of'~T[opt])
g11 <- ggplot(latrdf, aes(y = T050 - 273.15, x = latitude))+
  geom_point(size=0.5) + geom_smooth(method='gam', span=0.1)+
  geom_vline(xintercept = 0, linetype=2) +
  theme_bw() +
  xlab('latitude') + ylab('Optimum temperature(°C)') +
  ggtitle(ttl) + 
  coord_flip()+
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )
g11
library(gridExtra)
lmat <- matrix(c(rep(1, 9),0,2,2), nrow=1, byrow = T)
grid.arrange(g1, g11, layout_matrix = lmat)

##Mapping Pejur 90%##
mean_P75_dat <- read.csv('thermo_ph/results/mean_P75.csv', h=T)
head(mean_P75_dat)
mn_P75_dt <- data.frame(mean_P75_dat)
quantile(mn_P75_dt$V3, probs=c(0.9, 0.95, 0.99), na.rm=T)
mn_P75_dt$V3[mn_P75_dt$V3> 20] <- NA

cols = c( "green1", "green2", "green3",
         "royalblue1","royalblue2", "royalblue4",
         "yellow1", "yellow2", 
         "tomato1", "tomato2", "tomato3", "red4")
contnts <- readOGR('/media/karthik/ADATA HD720/global_veg/continent_shapefile/continent shapefile/','continent')
contnts <- subset(contnts, ! CONTINENT %in% 'Oceania' )
contnt_shp <- fortify(contnts)
ttl <- expression(~Temperature~range~"75% ("~T[0.75]~")")
g2 <- ggplot() + 
  geom_polygon(data = contnt_shp,
               aes(x = long, y = lat, group = group),
               fill = 'white', colour= 'black', size=0.1)+
  geom_tile(data = mn_P75_dt, aes(x=V1, y=V2, fill = (V3))) + 
  coord_equal() +
  scale_fill_gradientn("Temperature (°C)", colours = cols) +
  theme_bw()+xlab("")+ylab("")+ylim(-55, 85)+
  ggtitle(ttl)+
  theme(
    axis.text = element_text(size =15, color = 'black'),
    axis.ticks = element_blank(),
    panel.grid = element_blank(), 
    panel.background = element_rect(colour= 'black', fill = 'lightblue1'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title =  element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(size=20, colour='black', hjust=0.5)
  )##multi_returnrate 800 500

g2


##Plot latitudinal variation in P75
mn_p90_dt <- data.frame(mean_P75_dat)
colnames(mn_p90_dt) <- c("V1", "V2", "V3")
library(reshape2)
latrdf <- matrix(NA, ncol = 3, nrow = (length(intv)-1))
mn_p90_dt$V3[mn_p90_dt$V3 > 20] <- NA
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(mn_p90_dt, V2 >= intv[i] & V2 < intv[i+1])
  midlatr <- mean(c(intv[i], intv[i+1]))
  latrdf1[] <- Map(function(x) replace(x, is.infinite(x), NA), latrdf1)
  P75s <- melt(quantile(latrdf1$V3, probs=c(0.5, 0.9), na.rm=T))
  latrdf[i,] <- c(midlatr, P75s$value)
  print(i)
}

library(ggplot2)
latrdf <- data.frame(latrdf)
colnames(latrdf) <- c('latitude', 'T050', 'T090')

ttl <- expression('Latitudinal variation of'~T[0.75])
g22 <- ggplot(latrdf, aes(y = T050, x = latitude))+
  geom_point(size=0.5) + geom_smooth(method='loess', span=0.08)+
  theme_bw() +
  xlab('latitude') + ylab('performance range (75%)') +
  coord_flip()+
  ggtitle(ttl) + 
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )
g22

##thermal acclimatization##
whb  <- raster('/media/karthik/ADATA HD720/global_veg/rasters/whittaker_worldclim_ecoregion.tif')
plot(whb)
acc_dat <- read.csv('/media/karthik/karthik_cluster/thermo_ph/results/tmp_ac1.csv', h=F)
head(acc_dat)
qq <- quantile(acc_dat$V5, probs=c(0.9, 0.95))
acc_dat$V5[acc_dat$V5 > qq[1]] <- NA
cols = c( "green1", "green2", "green3",
          "royalblue1","royalblue2", "royalblue4",
          "yellow1", "yellow2", rep("yellow2", 5),
          "tomato1", "tomato2", "tomato3", "red4", rep("red4", 10))
g3 <- ggplot() + 
  geom_polygon(data = contnt_shp,
               aes(x = long, y = lat, group = group),
               fill = 'white', colour= 'black', size=0.1)+
  geom_tile(data = acc_dat, aes(x=V2, y=V3, fill = (V5))) + 
  coord_equal() +
  scale_fill_gradientn("Temperature (°C)", colours = cols) +
  theme_bw()+xlab("")+ylab("")+ylim(-55, 85)+
  ggtitle(ttl)+
  theme(
    axis.text = element_text(size =15, color = 'black'),
    axis.ticks = element_blank(),
    panel.grid = element_blank(), 
    panel.background = element_rect(colour= 'black', fill = 'lightblue1'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title =  element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(size=20, colour='black', hjust=0.5)
  )##multi_returnrate 800 500

g3

acc_dat$whb <- extract(whb, acc_dat[,2:3] )
acc_dat1 <- acc_dat[,c(2,3,5,6,11)]
acc_dat1 <- acc_dat1[complete.cases(acc_dat1),]
ggplot(subset(acc_dat1, V5 < 10), aes(x = factor(whb), y = V5, fill = as.factor(whb)))+
  geom_hline(yintercept=0, linetype=2)+
  geom_boxplot(outlier.shape = NA, width=0.5) + theme_bw()+
  ylab("Temperature (°C)")+ xlab("")+
  ggtitle("Difference b/w ambient & optimum temperature")+
  scale_x_discrete(breaks = 1:10, 
                   label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_fill_manual(values= c("green4", "olivedrab3", "tan4", "red4", "cyan",
                              "royalblue4", "grey70", "yellow", "sandybrown", "steelblue3"))+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )###800, 400, water_stress_biomes


ggplot(subset(acc_dat, V5 < 10 & !whb %in% NA), aes(x = factor(whb), y = V6, 
                                                    fill = factor(whb)))+
  geom_hline(yintercept=0, linetype=2)+
  geom_boxplot(outlier.shape = NA, width=0.5) + theme_bw()+
  ylab(expression("Kendall" ~ tau))+ xlab("")+
  ggtitle("Acclimatization trend")+
  scale_x_discrete(breaks = 1:10, 
                   label = c("Tropical \n rainforest","Tropical \n deciduous", 
                             "Tropical \n grassland", "Hot \n desert",
                             "Temperate \n rainforest","Temperate \n deciduous",
                             "Cold \n desert","Temperate \n grassland", "Taiga", "Tundra"))+
  scale_fill_manual(values= c("green4", "olivedrab3", "tan4", "red4", "cyan",
                              "royalblue4", "grey70", "yellow", "sandybrown", "steelblue3"))+
  theme(
    legend.position = "none",
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle= 90, size=14, colour='black',vjust=0.70),
    axis.text.y = element_text(size=10, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )###800, 400, water_stress_biomes





library(reshape2)
p90rs <- read.csv('thermo_ph/results/P90_trend.csv')
head(p90rs)
p90rs$sig <- p90rs$V4
p90rs$sig[p90rs$sig < 0.1] <- 10
p90rs$sig[!p90rs$sig %in% 10] <- 0
p90rs$sig <- (p90rs$sig*p90rs$V3)/10

lat_rng <- round(range(mn_topt_dt$V2, na.rm=T),2)
intv <- ceiling(seq(lat_rng[1], lat_rng[2], 2))-1
latrdf <-  data.frame(matrix(NA, ncol = 4, nrow = 1))
colnames(latrdf) <- c( 'V1', 'V2', 'sig', 'midlat')
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(p90rs, V2 >= intv[i] & V2 <  intv[i+1])
  latrdf1$midlat <- mean(intv[i], intv[i+1])
  latrdf <- rbind(latrdf, latrdf1[,c(2,3,6,7)])
  print(i)
}

ylb <- expression(Delta~T[0.9]~'/'~yr^-1)
latrdf <- latrdf[complete.cases(latrdf),]
g1 <- ggplot(latrdf, aes(y = sig, x = factor(midlat)))+
  geom_hline(yintercept = 0, linetype=2) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  coord_cartesian(ylim=c(-0.5, 1))+
  xlab('Latitude') + ylab(ylb)+
  ggtitle("Temporal trends in high performance \n temperature range") +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )


toptrs <- read.csv('thermo_ph/results/topt_trend_res.csv')
head(toptrs)
toptrs$sig <- toptrs$V4
toptrs$sig[toptrs$sig < 0.1] <- 10
toptrs$sig[!toptrs$sig %in% 10] <- 0
toptrs$sig <- (toptrs$sig*toptrs$V3)/10

lat_rng <- round(range(mn_topt_dt$V2, na.rm=T),2)
intv <- ceiling(seq(lat_rng[1], lat_rng[2], 2))-1
latrdf <-  data.frame(matrix(NA, ncol = 4, nrow = 1))
colnames(latrdf) <- c( 'V1', 'V2', 'sig', 'midlat')
for(i in 1:(length(intv)-1)){
  latrdf1 <- subset(toptrs, V2 >= intv[i] & V2 <  intv[i+1])
  latrdf1$midlat <- mean(intv[i], intv[i+1])
  latrdf <- rbind(latrdf, latrdf1[,c(2,3,6,7)])
  print(i)
}


latrdf <- latrdf[complete.cases(latrdf),]
ylb <- expression(Delta~T[opt]~'/'~yr^-1)
g2 <- ggplot(latrdf, aes(y = sig, x = factor(midlat)))+
  geom_hline(yintercept = 0, linetype=2) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  coord_cartesian(ylim=c(-0.2, 0.31))+ 
  xlab('Latitude') + ylab(ylb)+
  ggtitle("Temporal trends in Topt") + 
  theme_bw()+
  theme(
    axis.text = element_text(size = 16, colour= 'black'),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16, colour= 'black'),
    panel.grid = element_blank()
  )
gridExtra::grid.arrange(g2, g1, ncol=1)
