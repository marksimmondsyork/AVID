# CST NVD and NVE

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Combined IPD and AD")

###################################################
# Agg data

ipd.trials <- c("CLARITY","PROTOCOL S","PROTEUS")
npdr.trials <- c("PROTOCOL W","PANORAMA")
excl.trials <- c(npdr.trials,"RECOVERY","Sao Paulo A") 

ad.data <- read_excel("data/AD other outcome data.xlsx") %>%
  filter(!is.na(Trial), !Trial %in% excl.trials) %>%
  select(Trial:N, contains("NV")) %>%
  mutate(Class=ifelse(Arm=="PRP","PRP","Anti-VEGF"))

#####################################################
# IPD

ipd.data <- read_csv("data/AVID IPD follow-up.csv") %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12,
         Class=ifelse(Arm=="PRP","PRP","Anti-VEGF")) %>%
  select(Trial,PatEye,Arm,Class,Weeks,NVD_FU,NVE_FU,CST_FU)

