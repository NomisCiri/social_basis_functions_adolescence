## define fitting functions ---- 
ll_ddm2 <- function(x, dataBeh2, vd, hd, bd) {
  sourceCpp(here::here("tddm","tDDM_Rcpp_serial.cpp"))
  
  d_v = x[1] # weigh for primBF
  d_h = x[2] # weight for secBF
  d_b = x[3] # weight for bonus
  thres = x[4] # threshold
  nDT = x[5] # non-decision time
  tHin = x[6] # relative starting time of secBF vs. primBF
  bias = x[7] # starting point bias
  
  probs = NULL
  for (i in 1:length(vd)) {
    rts = ddm2_serial(d_v,d_h,d_b,thres,nDT,tHin,bias,vd[i],hd[i],bd[i],3000) # simulates rts
    rts = rts[rts!=0] #remove all rts that are 0
    xdens = density(rts, from=-5, to=5, n=1024, bw=0.11) #matches the prob. space of RTs from above
    idx = which(dataBeh2$primBF==vd[i] & dataBeh2$secBF==hd[i] & dataBeh2$bonus==bd[i])
    probs = c(probs, dt*xdens$y[dataBeh2$RTddm_pos[idx]])
  }
  
  probs[probs==0] = 1e-100
  return (-sum(log(probs)))
  
}




fitSub_partner <- function(s, dataBeh) {
  fits=matrix(0,1,10) 
  idx = which(dataBeh$subject_id==s)
  ##
  ##Adding this basis functions stuff
  ##
  ##
  #create primary basis functions
  dataBeh$primBF <- 0.5*dataBeh$S + 0.5*dataBeh$P - 0.5*dataBeh$O1 - 0.5*dataBeh$O2

  #create secondary basis functions partner
  dataBeh$secBF <- 0.5*dataBeh$P - 0.5*dataBeh$S - 0.5*dataBeh$Or + 0.5*dataBeh$Oi
  
  #log-transform reaction times
  dataBeh$logRT = log(dataBeh$rt)
  
  #assign negative RTs to choice == 0
  dataBeh$RTddm = dataBeh$rt
  dataBeh$RTddm[dataBeh$response_n == 0] = dataBeh$RTddm[dataBeh$response_n == 0] * -1
  
  #get number of trials
  ntrials = length(dataBeh$rt)
  
  # bin the prob. density space of RTs:
  # based on xpos, we can see the size of the time bins
  # if you want to adjust the bins, change the sequence length (-5,5) or 
  # length.out, which tells you how many divisions are put into the sequence
  
  xpos = seq(-5,5,length.out=1024)
  dt = xpos[2] - xpos[1]
  dataBeh$RTddm_pos = 0
  for (i in 1:ntrials) {
    dataBeh$RTddm_pos[i] = which.min(abs(xpos - dataBeh$RTddm[i]))
  }
  
  dataBeh2 = dataBeh[idx,]
  
  # surrogate data for ddm
  data1 = ddply(dataBeh2, .(primBF, secBF, bonus), summarize, acc = mean(response_n)) 
  
  vd = data1$primBF # performance score for primBF
  hd = data1$secBF # performance score for secBF
  bd = data1$bonus # performance score for bonus
  
  # boundaries for parameters
  lower <- c(-2,-2,-2,0.6,0.01,-1,-1)
  upper <- c(2,2,2,3,1,1,1)

  #
  #actually fit model
  #
  
  fit_s = DEoptim(ll_ddm2, lower, upper, DEoptim.control(itermax = 150), dataBeh2=dataBeh2, vd=vd, hd=hd, bd=bd) 
  
  
  fits[1,1:7] = fit_s$optim$bestmem #fitted parameters
  fits[1,8] = fit_s$optim$bestval #LL
  fits[1,9] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) #BIC
  fits[1,10] = 2*fit_s$optim$bestval + 2*length(lower) #AIC
  
  # tidy up
  fits = as.data.frame(matrix(unlist(fits),ncol=10,byrow=TRUE))
  names(fits)<-c("d_primBF", "d_secBF", "d_bonus", "thres", "nDT", "timesecBFIn", "bias", "LL", "BIC", "AIC")
  fits$ppt<-s
  return(fits)
}










fitSub_self <- function(s, dataBeh) {
  fits=matrix(0,1,10) 
  idx = which(dataBeh$subject_id==s)
  ##
  ##Adding this basis functions stuff
  ##
  ##
  #create primary basis functions
  dataBeh$primBF <- 0.5*dataBeh$S + 0.5*dataBeh$P - 0.5*dataBeh$O1 - 0.5*dataBeh$O2
  
  #create secondary basis functions partner
  dataBeh$secBF <- 0.5*dataBeh$S - 0.5*dataBeh$P - 0.5*dataBeh$Or + 0.5*dataBeh$Oi
  
  #log-transform reaction times
  dataBeh$logRT = log(dataBeh$rt)
  
  #compute RTs in seconds and assign negative RTs to choice == 0
  dataBeh$RTddm = dataBeh$rt/1000
  dataBeh$RTddm[dataBeh$response_n == 0] = dataBeh$RTddm[dataBeh$response_n == 0] * -1
  
  #get number of trials
  ntrials = length(dataBeh$rt)
  
  # bin the prob. density space of RTs:
  # based on xpos, we can see the size of the time bins
  # if you want to adjust the bins, change the sequence length (-5,5) or 
  # length.out, which tells you how many divisions are put into the sequence
  
  xpos = seq(-5,5,length.out=1024)
  dt = xpos[2] - xpos[1]
  dataBeh$RTddm_pos = 0
  for (i in 1:ntrials) {
    dataBeh$RTddm_pos[i] = which.min(abs(xpos - dataBeh$RTddm[i]))
  }
  # subset for one participant (probably a wierd matlab way to do it.) 
  dataBeh2 = dataBeh[idx,]

  # surrogate data for ddm
  data1 = ddply(dataBeh2, .(primBF, secBF, bonus), summarize, acc = mean(response_n)) 
  
  vd = data1$primBF # performance score for primBF
  hd = data1$secBF # performance score for secBF
  bd = data1$bonus # performance score for bonus
  
  # boundaries for parameters
  lower <- c(-2,-2,-2,0.6,0.01,-1,-1)
  upper <- c(2,2,2,3,1,1,1)
  
  #
  #actually fit model
  #
  
  fit_s = DEoptim(ll_ddm2, lower, upper, DEoptim.control(itermax = 150), dataBeh2=dataBeh2, vd=vd, hd=hd, bd=bd) 
  
  
  fits[1,1:7] = fit_s$optim$bestmem #fitted parameters
  fits[1,8] = fit_s$optim$bestval #LL
  fits[1,9] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) #BIC
  fits[1,10] = 2*fit_s$optim$bestval + 2*length(lower) #AIC
  
  # tidy up
  fits = as.data.frame(matrix(unlist(fits),ncol=10,byrow=TRUE))
  names(fits)<-c("d_primBF", "d_secBF", "d_bonus", "thres", "nDT", "timesecBFIn", "bias", "LL", "BIC", "AIC")
  fits$ppt<-s
  return(fits)
}

