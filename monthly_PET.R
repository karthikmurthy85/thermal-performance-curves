petfl <- list.files('thermo_ph/PET_daily/', full.names=T)

library(foreach)
library(doParallel)
registerDoParallel(39)

foreach(i = 1:length(petfl))%dopar%{

	library(raster)
	
	pt1 <- brick(petfl[i])
	dats <- getZ(pt1)
	yr <- format(dats, "%Y")[1]
	mpt1 <- zApply(pt1, by = format(dats,"%m"), fun=sum)
	fnm <- paste0('thermo_ph/PET_monthly_', yr, '.grd')
	writeRaster(mpt1, fnm)

}
	
