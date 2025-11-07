# clear workspace
rm(list=ls())

# set current wd 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("tDDM_selfcondition.R", echo = TRUE)
source("tDDM_partnercondition.R", echo = TRUE)
