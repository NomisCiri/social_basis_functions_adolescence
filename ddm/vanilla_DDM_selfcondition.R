###This script fits a tDDM to the self condition data

# clear workspace
rm(list=ls())

# set current wd and seed
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(007)
here::i_am("./ddm/vanilla_DDM_selfcondition.R")
# Install required packages if necessary:
want = c("DEoptim", "Rcpp", "plyr", "parallel", "BH", "tidyverse")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }
# Now load them all
lapply(want, require, character.only = TRUE)

RcppParallel::setThreadOptions(numThreads = 1) #this is critical for running on Mac OS.

### --- how to run the tDDM:
# load the c++ file containing the functions to simulate the time-varying DDM
sourceCpp(here::here("ddm","ddm_models_c++","vanilla_ddm_Rcpp.cpp"))

# read in behavioral data
dataBeh=read_csv(here::here("data","sod_data.csv"))%>%filter(!is.na(version))%>%
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
    relevant=case_when(self_partner==1 ~ (S-Or)+bonus,
                       self_partner==2 ~ (P-Or)+bonus
    ),
    irrelevant=case_when(self_partner==1 ~ P-Oi,
                         self_partner==2 ~ S-Oi
    )
  )

#filter for self condition
dataBeh <- filter(dataBeh, self_partner == 1)

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
ll_vanilla_ddm2 <- function(x, dataBeh2, vd) {
  
  d_v = x[1] # weight for primBF
  thres = x[2] # threshold
  nDT = x[3] # non-decision time
  bias = x[4] # starting point bias
  
  probs = NULL
  for (i in 1:length(vd)) {
    rts = vanilla_ddm2_parallel(d_v,thres,nDT,bias,vd[i],3000) # simulates rts
    rts = rts[rts!=0] #remove all rts that are 0
    xdens = density(rts, from=-5, to=5, n=1024, bw=0.11)
    idx = which(dataBeh2$relevant==vd[i]) 
    probs = c(probs, dt*xdens$y[dataBeh2$RTddm_pos[idx]])
  }
  
  probs[probs==0] = 1e-100
  return (-sum(log(probs)))
  
}

fitSub <- function(s, dataBeh) {
  fits=matrix(0,1,8)
  idx = which(dataBeh$subject_id==s)
  dataBeh2 = dataBeh[idx,]
  
  data1 = ddply(dataBeh2, .(relevant), summarize, acc = mean(response_n)) 
  
  vd = data1$relevant#relevant DV
  
  # boundaries for parameters
  lower <- c(-2,0.6,0.01,-1)
  upper <- c(2,3,1,1)
  
  fit_s = DEoptim(ll_vanilla_ddm2, lower, upper, DEoptim.control(itermax = 200), dataBeh2=dataBeh2, vd=vd) 
  fits[1,1:length(lower)] = fit_s$optim$bestmem #fitted parameters
  fits[1,length(lower)+1] = fit_s$optim$bestval #LL
  fits[1,length(lower)+2] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) #BIC
  fits[1,length(lower)+3] = 2*fit_s$optim$bestval + 2*length(lower) #AIC
  fits[1,length(lower)+4] = s
  
  return(fits)
}

#get subject ids
inputs = unique(dataBeh$subject_id)

#'mc.cores' > 1 is not supported for windows
fits = mclapply(inputs, fitSub, mc.cores = 1, dataBeh=dataBeh)

#save now in case unlist fails
fitsF=fits
write.csv(fitsF, file = here::here("ddm","fits","fits_vanilla_DDM_self.csv"))

fits = as.data.frame(matrix(unlist(fits),ncol=8,byrow=TRUE))
names(fits)<-c("d_rel","thres", "nDT", "bias", "LL", "BIC", "AIC","ppt")
#will overwrite previous save, but that is intended
fitsF=fits
write.csv(fitsF, file = here::here("ddm","fits","fits_vanilla_DDM_self.csv"))
