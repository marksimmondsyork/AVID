# impact of additional treatment on 
# anti-VEGF effect

library(tidyverse)
library(magrittr)
library(readxl)

library(lme4)

setwd ("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/IPD")

#source("R code/support code/IPD data BCVA set-up.R")
#source("R code/support code/IPD data other outcomes set-up.R")

add.trt <- read_csv("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/IPD/Protocol S/Addtional treatment data - Prot S 1707.csv") %>%
  mutate(PatEye=paste("PROTOCOL S -", PatEye))

add.trt.ps <- read_csv("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/IPD/Protocol S/Addtional treatment data - Prot S 1707.csv") %>%
  mutate(PatEye=paste("PROTOCOL S -", PatEye))

add.trt.pr <- read_excel("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/IPD/PROTEUS/Edited PROTEUS data for AVID 250723.xlsx",sheet="Treatments")

add.trt.cl <- read_csv("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/IPD/CLARITY/AVID data/CLARITY additional treatment data.csv")

#############################################
# PROTOCOL S
# summarise additional treatments

neyes <- add.trt.ps %>%
  filter(!is.na(Treatment)) %>%
  group_by(PatEye) %>%
  slice_head(n=1) %>%
  group_by(Treatment) %>%
  summarise(Neyes=n())

summ.add.trt.byeye <- add.trt.ps %>%
  filter(!is.na(Treatment)) %>%
  group_by(PatEye,Treatment) %>%
  summarise(across(contains("Add"), ~ifelse(sum(.x)>0, 1,0))) %>%
  group_by(Treatment) %>%
  summarise(across(contains("Add"), sum))

summ.add.trt.nprocs <- add.trt.ps %>%
  filter(!is.na(Treatment)) %>%
  group_by(Treatment) %>%
  summarise(across(contains("Add"), sum))

summ.add.trt <- add.trt.ps %>%
  filter(!is.na(Treatment)) %>%
  group_by(PatEye,Treatment) %>%
  summarise(across(contains("Add"), ~ifelse(sum(.x)>0, 1,0))) %>%
  ungroup() %>%
  mutate(Add.laser.1=ifelse(Add.Laser==1 | Add.PRP==1, 1, 0),
         Add.AVEGF.2=ifelse(Add.AVEGF==1 | Add.AVEGF.1==1, 1, 0)) %>%
  mutate(Add.trt.comb=case_when((Add.laser.1==1 & Add.AVEGF.2==1) ~ "Both",
                               Add.AVEGF.2==1 ~ "Anti-VEGF",
                               Add.laser.1==1 ~ "Laser photo.",
                               TRUE ~ "Neither")) 

sum.add.trt.2 <- summ.add.trt %>%
  group_by(Treatment,Add.trt.comb) %>%
  summarise(Npats=n())

summ.add.trt.pc <- add.trt.ps %>%
  filter(!is.na(Treatment)) %>%
  mutate(Add.laser.1=ifelse(Add.Laser==1 | Add.PRP==1, 1, 0),
         Add.AVEGF.2=ifelse(Add.AVEGF==1 | Add.AVEGF.1==1, 1, 0)) %>%
  mutate(Add.trt.comb=case_when((Add.laser.1==1 & Add.AVEGF.2==1) ~ "Both",
                                Add.AVEGF.2==1 ~ "Anti-VEGF",
                                Add.laser.1==1 ~ "Laser photo.",
                                TRUE ~ "Neither")) %>%
  group_by(Treatment,Add.trt.comb) %>%
  summarise(Nprocs=n())

# by DME status

dme.data <- ipd.oc %>%
  filter(Trial=="PROTOCOL S") %>%
  select(PatEye,Treatment,TimeObs,DME_FU,VH_FU) %>%
  group_by(PatEye) %>%
  summarise(ever.DME=max(DME_FU),
            ever.VH=max(VH_FU))

add.trt.dme <- left_join(summ.add.trt,dme.data) 

summ.trt.dme.byeye <- add.trt.dme %>%
  group_by(Treatment,ever.DME,Add.trt.comb) %>%
  summarise(Npats=n())

############################################
# PROTEUS summarise treatments

pr.data <- ipd.base %>%
  filter(Trial=="PROTEUS") %>%
  select(PatientID:Treatment) %>%
  full_join(add.trt.pr)

npats.pr <- pr.data %>%
  filter(!is.na(Treatment)) %>%
  group_by(PatientID) %>%
  slice_head(n=1) %>%
  group_by(Treatment) %>%
  summarise(Neyes=n())

summ.pr.bypat <- pr.data %>%
  filter(!is.na(Treatment)) %>%#
  filter(Month>2) %>%
  group_by(PatientID,Treatment) %>%
  summarise(across(PRP:AVEGF, ~ifelse(sum(.x)>0, 1,0))) %>%
  ungroup() %>%
  mutate(Add.Trt=case_when(AVEGF==1 & PRP==1 ~ "Both",
                           AVEGF==1 ~ "Anti-VEGF",
                           PRP==1 ~ "PRP",
                           TRUE ~ "Neither")) %>%
  group_by(Treatment,Add.Trt) %>%
  summarise(npats=n())

summ.pr.nprocs <- pr.data %>%
  filter(!is.na(Treatment)) %>%#
  filter(Month>2) %>%
  group_by(Treatment) %>%
  summarise(across(PRP:AVEGF, sum))

##############################################
# CLARITY

cl.data <-  add.trt.cl %>%
  filter(!is.na(Add.trt)) %>%
  filter((Treatment==1 & Week>8) | (Treatment==0 & Week>4)) 

cl.data.2 <- cl.data %>%
  pivot_wider(id_cols=PatientID:Week,names_from=Add.trt,values_from=Add.trt) %>%
  select(-None) %>%
  mutate(across(PRP:`Aflibercept+PRP`, ~ifelse(is.na(.x),0,1)))

npats.cl <- cl.data %>%
  group_by(PatientID) %>%
  slice_head(n=1) %>%
  group_by(Treatment) %>%
  summarise(Npats=n())

summ.cl.bypat <- cl.data.2 %>%
  group_by(PatientID,Treatment) %>%
  summarise(across(PRP:`Aflibercept+PRP`, ~ifelse(sum(.x)>0, 1,0))) %>%
  group_by(Treatment) %>%
  summarise(across(PRP:`Aflibercept+PRP`, sum))
  
summ.cl.byproc <- cl.data.2 %>%
  group_by(Treatment) %>%
  summarise(across(PRP:`Aflibercept+PRP`, sum))

#############################################
# ETDRS by treatments

add.trt.2 <- add.trt %>%
  filter(!is.na(Treatment)) %>%
  group_by(PatEye,Treatment) %>%
  summarise(across(contains("Add"), ~ifelse(sum(.x)>0, 1,0))) %>%
  ungroup() %>%
  mutate(Add.laser.1=ifelse(Add.Laser==1 | Add.PRP==1, 1, 0),
         Add.AVEGF.2=ifelse(Add.AVEGF==1 | Add.AVEGF.1==1, 1, 0)) %>%
  mutate(Add.trt=case_when((Add.laser.1==1 & Add.AVEGF.2==1) ~ "Both",
                               Add.AVEGF.2==1 ~ "Anti-VEGF",
                               Add.laser.1==1 ~ "Laser",
                               TRUE ~ "Neither")) 

bcva.ps <- filter(ipd.bcva,Trial=="PROTOCOL S") %>%
  left_join(add.trt.2)

change.etdrs <- bcva.ps %>%
  select(PatEye,Arm,TimeObs,etdrs,etdrs.cfb,Add.trt) %>%
  group_by(Arm,TimeObs,Add.trt) %>%
  summarise(etdrs.mean=mean(etdrs,na.rm=T),
            etdrs.sd=sd(etdrs,na.rm=T),
            etdrs.cfb.mean=mean(etdrs.cfb,na.rm=T),
            etdrs.cfb.sd=sd(etdrs.cfb,na.rm=T),
            Narm=n()) %>%
  mutate(etdrs.sem=etdrs.sd / sqrt(Narm),
         etdrs.lower=etdrs.mean-1.96*etdrs.sem,
         etdrs.upper=pmin(etdrs.mean+1.96*etdrs.sem, 100),
         etdrs.cfb.sem=etdrs.cfb.sd / sqrt(Narm),
         etdrs.cfb.lower=etdrs.cfb.mean-1.96*etdrs.cfb.sem,
         etdrs.cfb.upper=etdrs.cfb.mean+1.96*etdrs.cfb.sem) %>%
  mutate(Year=TimeObs/12) %>%
  filter(!is.na(etdrs.sem))

ggplot(change.etdrs,aes(x=Year,y=etdrs.cfb.mean,
                        colour=Arm)) +
  geom_smooth(method="lm") +
  facet_wrap(~Add.trt) +
  scale_x_continuous("Years",breaks=seq(0,6,1)) +
  scale_y_continuous("Mean change in ETDRS letters",breaks=seq(-80,80,5))

