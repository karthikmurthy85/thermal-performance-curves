library(utils)

r = seq(0.5, 5.5, 0.5)
q = 1/seq(10000, 50000, 5000)
aC = seq(0.001, 0.5, 0.05)
eC = seq(0.01, 0.3, 0.05)
hC = seq(0.05, 0.4, 0.05)
mC = seq(0.1, 0.5, 0.1)

aP = seq(0.05, 0.5, 0.05)
eP = seq(0.05, 0.5, 0.05)
hP = seq(0.05, 0.5, 0.05)
mP = seq(0.1, 0.7, 0.1)

d1 <- expand.grid(r = r, q = q,
                  aC = aC, eC = eC, hC = hC, mC = mC,
                  aP = aP, eP = eP, hP = hP, mP = mP)

#length(r)*length(q)*length(aC)*length(eC)*length(hC)*length(mC)*length(aP)*length(eP)*length(hP)#*length(mP)

diff <- d1$eC - d1$eP
st1 <- (d1[-which(diff > 0),])
rm(d1)
#write.csv(st1, 'ss_parm.csv')

split_csv <- function(dat, n, fpath){
  library(stringr)
  nr = c(seq(0, nrow(dat),n),nrow(dat))
  for(z in 1:(length(nr)-1)){
    dta <- dat[((nr[z]+1):nr[z+1]),]
    id <- str_pad(z, pad='0', width = nchar(nrow(dat)))
    fnm <- paste0(fpath, id, '.csv')
    #print(fnm)
    write.csv(dta, fnm, row.names = F )
  }
  
}

split_csv(st1, 1000000, 'RCP_param_space/ss_parm')
