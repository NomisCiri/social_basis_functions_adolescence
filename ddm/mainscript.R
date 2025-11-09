# clear workspace
rm(list=ls())

# set current wd 
here::i_am("./ddm/mainscript.R")

source(here::here("ddm","tDDM_selfcondition.R"))
source(here::here("ddm","tDDM_partnercondition.R"))
