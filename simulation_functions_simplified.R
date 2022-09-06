#Essential finctions
cc_score <- function(dist, fitted){
  
  num<-sum((dist-fitted)^2)
  denom1<-sum((dist-mean(dist))^2)
  denom2<-sum((fitted-mean(fitted))^2)
  denom3<-length(dist)*((mean(dist)-mean(fitted))^2)
  CC_score<- 1-((num)/(denom1+denom2+denom3))
  return(CC_score)
}

min_max_fun <- function(x){
  return((x -  min(x, na.rm=T))/(max(x,na.rm=T) - min(x, na.rm=T)))
}


schoolfield_fun<- function(np2){
  
  Wx <- function(x){
    boltz <- 8.617e-5
    numerator <- x[1]*exp((-x[2]/boltz)*((1/temp)-(1/x[5])))
    denominator <- 1 + (x[2]/(x[3] -x[2]))*exp((x[3]/boltz)*((1/x[4]) - (1/temp)))
    erg1 <- numerator/denominator
    erg <- sum((erg1 - np2$npp)^2)
    return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
  }
  
  
  Wt <- function(B0, E, Ed, Tpk, Tref){
    boltz <- 8.617e-5
    numerator <- B0*exp((-E/boltz)*((1/temp)-(1/Tref)))
    denominator <- 1 + (E/(Ed -E))*exp((Ed/boltz)*((1/Tpk) - (1/temp)))
    res <- numerator/denominator
    return(res)
  }
  
  Wt1 <- function(B0, E, Ed, Tpk, Tref, temp1){
    boltz <- 8.617e-5
    numerator <- B0*exp((-E/boltz)*((1/temp1)-(1/Tref)))
    denominator <- 1 + (E/(Ed -E))*exp((Ed/boltz)*((1/Tpk) - (1/temp1)))
    res <- numerator/denominator
    return(res)
  }
  
  
  a1 <- min_max_fun(np2$npp)
  np2$npp <- a1 + 0.01
  
  np2$tmp <- np2$tmp + 273.15
  #plot(np2)
  temp <- np2$tmp 
  topt <- np2[order(np2$npp),]
  
  boltz <- 8.617e-5
  
  Tpk <- median(topt$tmp[(nrow(topt)-2): nrow(topt)]) 
  B0dat <- subset(np2, tmp < Tpk)
  B0dat1 <- B0dat[order(B0dat$npp),]
  B0 <- mean((B0dat1$npp)[1:3], na.rm=T)
  Tref <- mean((B0dat1$tmp)[1:3], na.rm=T) 
  
  tr <- abs(mean(np2$tmp) - 273.12)
  
  tt <- seq(min(np2$tmp),max(np2$tmp), 1)
  E = 0.1*tr; Ed = 0.1*tr
  
  optimal <- optim(c(B0, E, Ed, Tpk, Tref),fn=Wx,gr=NULL,method="L-BFGS-B",
                   lower=c(B0-0.01, 0.001, 0.001, (Tpk-1),(Tref -5)),
                   upper=c(B0+0.01, 40, 40, (Tpk+1),(Tref +5)),
                   control=list(maxit = 5000, pgtol = 1e-10,
                                ndeps = c(0.01, 1e-10, 1e-10, 0.1, 0.1),
                                lmm = 200))
  print(optimal$par)
  B01 <- optimal$par[1]
  E1 <- optimal$par[2]
  Ed1 <- optimal$par[3]
  Tpk1 <- optimal$par[4]
  Tref1 <- optimal$par[5]
  tt <- seq(min(np2$tmp),max(np2$tmp), 1)
  cc_score(np2$npp, Wt(B01, E1, Ed1, Tpk1, Tref1))
  
  delE <- E1/2.5; delEd <- Ed1/2.5
  optimal <- optim(c(B01, E1, Ed1, Tpk1, Tref1),fn=Wx,gr=NULL,method="L-BFGS-B",
                   lower=c(B0-0.01, E1 - delE, Ed1 - delEd, (Tpk-1),(Tref -5)),
                   upper=c(B0+0.01, E1 + delE, Ed1 + delEd, (Tpk+1),(Tref +5)),
                   control=list(maxit = 5000, pgtol = 1e-10,
                                ndeps = c(0.01, 1e-10, 1e-10, 0.1, 0.1),
                                lmm = 200))
  print(optimal$par)
  B01 <- optimal$par[1]
  E1 <- optimal$par[2]
  Ed1 <- optimal$par[3]
  Tpk1 <- optimal$par[4]
  Tref1 <- optimal$par[5]
  tt <- seq(min(np2$tmp),max(np2$tmp), 1)
  
  fit1 <- Wt1(B01, E1, Ed1, Tpk1, Tref1, tt) #255:310)
  #fit1
  # plot(np2)
  # points(  tt, fit1, col='blue', ty='l')
  cc_optim2 <- cc_score(np2$npp, Wt(B01, E1, Ed1, Tpk1, Tref1))
  
  temp = np2$tmp
  tst1 <- data.frame(per = Wt1(B01, E1, Ed1, Tpk1, Tref1, temp),
                     temp = temp)
  
  temp1 <- seq(253,323,by=0.01)
  tst <- data.frame(per = Wt1(B01, E1, Ed1, Tpk1,Tref1, temp1),
                    temp = temp1)
  
  np090 <- subset(tst, per > 0.9*max(tst$per, na.rm=T))
  diff(range(np090$temp))
  
  np075 <- subset(tst, per > 0.75*max(tst$per, na.rm=T))
  diff(range(np075$temp))
  
  tpt <- subset(tst, per %in% max(tst$per, na.rm=T))
  
  res1 <- c(tpt$temp,diff(range(np090$temp)),
            diff(range(np075$temp)), cc_optim2)
  res2 <- tst1
  res <- list(res1, res2)
  return(res)
  
}



logistic_inv <- function(ndvi){
  time <- 1:length(ndvi)
  
  Wx <- function(x){
    erg <- sum( ((x[1] + (x[2]/(1 + exp((x[3]-time)/x[4])))) - (ndvi))^2)
    return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
  }
  
  Wt <- function(Init, Asym, xmid, scal) {
    erg <-  Init + (Asym/(1 + exp((xmid-time)/scal)))
    return(erg)
  }
  
  Init <- ndvi[1]
  Asym <- max(ndvi)
  xmid <- 6
  scal <- -0.9
  
  xmids <- seq(2, length(ndvi)-1, by= 3)
  par_dat <- matrix(data=NA, nrow=length(xmids), ncol=6)
  for(z in 1:length(xmids)){
    optimal1 <- try(optim(c(Init, Asym, xmid, scal),fn=Wx,gr=NULL,
                          method="L-BFGS-B",
                          lower=c(0.001,0, xmids[z], -7), 
                          upper=c(1, 1,(length(ndvi)), -0.1),
                          control=list(maxit = 5000, pgtol = 1e-10, 
                                       ndeps = c(1e-10, 1e-10, 1e-10, 1e-10), 
                                       lmm = 200)))
    Init1 <- optimal1$par[1]
    Asym1 <- optimal1$par[2]
    xmid1 <- optimal1$par[3]
    scal1 <- optimal1$par[4]
    cc1 <- cc_score(ndvi, Wt(Init1, Asym1, xmid1, scal1))
    par_dat[z,] <- c(z, cc1, c(Init1, Asym1, xmid1, scal1))
    print(z)
  }
  
  mxdt <- subset(data.frame(par_dat), X2 %in% max(X2, na.rm=T))
  # plot(ndvi, ty ='o')
  fit <- Wt(mxdt[1,3], mxdt[1,4], mxdt[1,5], mxdt[1,6])
  # points( fit, ty = 'o', col = 'red')
  res <- c(max(par_dat[,2]),fit )
  return(res)
}

logistic <- function(ndvi){
  
  time <- 1:length(ndvi)
  
  Wx <- function(x){
    erg <- sum( ((x[1] + (x[2]/(1 + exp((x[3]-time)/x[4])))) - (ndvi))^2)
    return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
  }
  
  Wt <- function(Init, Asym, xmid, scal) {
    erg <-  Init + (Asym/(1 + exp((xmid-time)/scal)))
    return(erg)
  }
  
  Init <- ndvi[1]
  Asym <- max(ndvi)
  xmid <- 6
  scal <- 0.9
  
  xmids <- seq(2, length(ndvi)-1, by= 3)
  par_dat <- matrix(data=NA, nrow=length(xmids), ncol=6)
  for(z in 1:length(xmids)){
    optimal1 <- try(optim(c(Init, Asym, xmid, scal),fn=Wx,gr=NULL,
                          method="L-BFGS-B",
                          lower=c(0.001,0, xmids[z], 0.1), 
                          upper=c(1, 1,(length(ndvi)), 7),
                          control=list(maxit = 5000, pgtol = 1e-10, 
                                       ndeps = c(1e-10, 1e-10, 1e-10, 1e-10), 
                                       lmm = 200)))
    Init1 <- optimal1$par[1]
    Asym1 <- optimal1$par[2]
    xmid1 <- optimal1$par[3]
    scal1 <- optimal1$par[4]
    cc1 <- cc_score(ndvi, Wt(Init1, Asym1, xmid1, scal1))
    par_dat[z,] <- c(z, cc1, c(Init1, Asym1, xmid1, scal1))
    print(z)
  }
  prdt <- subset(data.frame(par_dat), X2 %in% max(X2))
  # plot(ndvi, ty='o')
  fit <- Wt(prdt$X3[1], prdt$X4[1], prdt$X5[1], prdt$X6[1])
  # points(fit, ty='o', col='red')
  res <- c(max(par_dat[,2]), fit)
  return(res)
}

quan_fun <- function(tmp_intr, df1){
  tmp_int <- seq(-40, 70, tmp_intr)
  npp_q <- matrix(NA, ncol = 4, nrow = (length(tmp_int)-1))
  for(i in 1:(length(tmp_int)-1)){
    tmpint <- subset(df1, TEMP > tmp_int[i] & TEMP < tmp_int[i+1])
    if(nrow(tmpint)> 1){
      npp_q[i,] <- c(mean(c(tmp_int[i], tmp_int[i+1])), 
                     quantile(tmpint$NDVI, probs = c(0.85, 0.9, 0.95)))
    }
    else{
      npp_q[i,] <- c(NA,NA,NA,NA)
    }
    print(i)
  }
  return(npp_q)
}



