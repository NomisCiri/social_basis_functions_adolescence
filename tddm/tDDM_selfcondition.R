###This script fits a tDDM to the self condition data

# clear workspace
rm(list=ls())

# set current wd and seed
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(007)

# Install required packages if necessary:
want = c("DEoptim", "Rcpp", "plyr", "parallel", "BH", "tidyverse")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }
# Now load them all
lapply(want, require, character.only = TRUE)
here::i_am("./tddm/tDDM_selfcondition.R")

RcppParallel::setThreadOptions(numThreads = 1) #this is critical for running on Mac OS.
### --- how to run the tDDM:
# load the c++ file containing the functions to simulate the time-varying DDM
sourceCpp("tDDM_Rcpp.cpp")

# read in behavioral data
dataBeh=read_csv(here::here("data","results","sod_data.csv"))%>%filter(!is.na(version))%>%
  mutate(group=ifelse(age>18,"adults","adolescents"))%>%rowwise()%>%
  mutate(self_partner=case_when(
    (decision_type==1 | decision_type == 2)~1, # 1 and 2 are self trials
    (decision_type==3 | decision_type == 4)~2, # 2 and 3 are partner trials
    decision_type==5 ~ 3 #. 5 are group trials
  ))%>%
  mutate(Or=case_when(
    (decision_type==1 | decision_type == 3)~O1, #i.e., 'S vs. O1' or  'P vs. O1', makes O1, relevant opponent
    (decision_type==2 | decision_type == 4)~O2, # i.e., 'P vs.O2' or 'P vs. 02' makes O2 the relevant opponent
    decision_type==5 ~ O1# here it does not matter
  ),
  Oi=case_when(
    (decision_type==1 | decision_type == 3)~O2, #i.e., 'S vs. O1' or  'P vs. O1', makes O2, irrelevant opponent
    (decision_type==2 | decision_type == 4)~O1, # i.e., 'P vs.O2' or 'P vs. 02' makes O1 the irrelevant opponent
    decision_type==5 ~ O2# does not matter
  )
  )%>%#select(-S_perf,-P_perf,-O1_perf,-O2_perf)%>%
  mutate(
    relevant=case_when(self_partner==1 ~ (S-Or),#+bonus,
                       self_partner==2 ~ (P-Or)#+bonus
    ),
    irrelevant=case_when(self_partner==1 ~ P-Oi,
                         self_partner==2 ~ S-Oi
    )
  )

#filter for self condition
dataBeh <- filter(dataBeh, self_partner == 1)

#create primary basis functions
dataBeh$primBF <- 0.5*dataBeh$S + 0.5*dataBeh$P - 0.5*dataBeh$O1 - 0.5*dataBeh$O2#dataBeh$relevant#0.5*dataBeh$S + 0.5*dataBeh$P - 0.5*dataBeh$O1 - 0.5*dataBeh$O2

#create secondary basis functions
dataBeh$secBF <- 0.5*dataBeh$S - 0.5*dataBeh$P - 0.5*dataBeh$Or + 0.5*dataBeh$Oi#dataBeh$irrelevant#0.5*dataBeh$S - 0.5*dataBeh$P - 0.5*dataBeh$Or + 0.5*dataBeh$Oi

#log-transform reaction times
dataBeh$logRT = log(dataBeh$rt)

#assign negative RTs to choice == 0
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

## define fitting functions ---- 
ll_ddm2 <- function(x, dataBeh2, vd, hd, bd) {
  
  d_v = x[1] # weight for primBF
  d_h = x[2] # weight for secBF
  d_b = x[3] # weight for the bonus
  thres = x[4] # threshold
  nDT = x[5] # non-decision time
  tHin = x[6] # relative starting time of secBF vs. primBF
  bias = x[7] # starting point bias
  
  probs = NULL
  for (i in 1:length(vd)) {
    rts = ddm2_parallel(d_v,d_h,d_b,thres,nDT,tHin,bias,vd[i],hd[i],bd[i],3000) # simulates rts
    rts = rts[rts!=0] #remove all rts that are 0
    xdens = density(rts, from=-5, to=5, n=1024, bw=0.11)
    idx = which(dataBeh2$primBF==vd[i] & dataBeh2$secBF==hd[i] & dataBeh2$bonus==bd[i]) 
    probs = c(probs, dt*xdens$y[dataBeh2$RTddm_pos[idx]])
  }
  
  probs[probs==0] = 1e-100
  return (-sum(log(probs)))
  
}

fitSub <- function(s, dataBeh) {
  fits=matrix(0,1,11)
  idx = which(dataBeh$subject_id==s)
  dataBeh2 = dataBeh[idx,]
  
  data1 = ddply(dataBeh2, .(primBF, secBF, bonus), summarize, acc = mean(response_n)) 
  
  vd = data1$primBF # performance score for primBF
  hd = data1$secBF # performance score for secBF
  bd = data1$bonus # performance score for bonus
  
  # boundaries for parameters
  lower <- c(-2,-2,-2,0.6,0.01,-1,-1)
  upper <- c(2,2,2,3,1,1,1)
  
  fit_s = DEoptim(ll_ddm2, lower, upper, DEoptim.control(itermax = 200,NP=100), dataBeh2=dataBeh2, vd=vd, hd=hd, bd=bd) 
  fits[1,1:7] = fit_s$optim$bestmem #fitted parameters
  fits[1,8] = fit_s$optim$bestval #LL
  fits[1,9] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) #BIC
  fits[1,10] = 2*fit_s$optim$bestval + 2*length(lower) #AIC
  fits[1,11] = s
  
  return(fits)
}

#get subject ids
inputs = unique(dataBeh$subject_id)

#'mc.cores' > 1 is not supported for windows
fits = mclapply(inputs, fitSub, mc.cores = 1, dataBeh=dataBeh)

#save now in case unlist fails
fitsF=fits
write.csv(fitsF, file = "fits_tDDM_self.csv")

fits = as.data.frame(matrix(unlist(fits),ncol=11,byrow=TRUE))
names(fits)<-c("d_primBF", "d_secBF", "d_bonus", "thres", "nDT", "timesecBFIn", "bias", "LL", "BIC", "AIC","ppt")
#will overwrite previous save, but that is intended
fitsF=fits
write.csv(fitsF, file = "fits_tDDM_self.csv")
