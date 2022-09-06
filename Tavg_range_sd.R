library(raster)
setwd('/media/karthik/karthik_cluster/')
ndvi_veg <- brick('thermo_ph/rasters/global_veg9km.grd')
ndvi_veg_df <- data.frame(rasterToPoints(ndvi_veg))
ndvi_veg_df$layer <- NULL
temp_ras <- brick('thermo_ph/rasters/global_monthly_Tavg.grd')
ll <- 10; yrs <- 1:34

ids <- seq(6, 414, 12) 
res_temp_ras <- temp_ras[[1]]
var_temp_ras <- temp_ras[[1]]
for(i in 1:(length(yrs)-ll)){
  lids <- c((ids[i]-5):ids[i], (ids[i]+1):ids[ll+i+1])
  vr_tmp1 <- temp_ras[[lids]]
  vr_tmp1[vr_tmp1 < 270] <- NA
  vr_tmp <- calc(vr_tmp1, function(x) sd(x, na.rm=T))
  
  var_temp_ras <- stack(var_temp_ras, vr_tmp)
  
  ##For calculating range##
  # temp_ras_iter <- calc(temp_ras[[lids]], function(x) range(x))
  # rg1 <- temp_ras_iter[[1]]
  # rg1[rg1 < 270] <- 273.15
  # rg <- temp_ras_iter[[2]] - rg1
  # rg[rg < 0] <- NA
  # res_temp_ras <- stack(res_temp_ras, rg)
  print(i)
  
}
#res_temp_ras <- res_temp_ras[[-1]]
var_temp_ras <- var_temp_ras[[-1]]
writeRaster(res_temp_ras, 'thermo_ph/rasters/global_monTavg_range.grd')


##SD of the monthly temperatures##

sd_temp_ras <- temp_ras[[1]]

sq <- seq(0, 24, 4)
sq

for(j in 1:(length(sq)-1)){
  
  id1 <- sq[j] + 1
  id2 <- sq[j+1]
  
  library(foreach)
  library(doParallel)
  registerDoParallel(4)
  
  sd_temp_ras1 <- foreach(i = id1:id2, .inorder=T, .combine='stack')%dopar%{
    library(raster)
    lids <- c((ids[i]-5):ids[i], (ids[i]+1):ids[ll+i+1])
    vr_tmp1 <- temp_ras[[lids]]
    vr_tmp1[vr_tmp1 < 270] <- NA
    vr_tmp <- calc(vr_tmp1, function(x) sd(x, na.rm=T))
    vr_tmp
  }
  
  sd_temp_ras <- stack(sd_temp_ras, sd_temp_ras1)
  print(j)
  
}