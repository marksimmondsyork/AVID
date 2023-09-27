# use 2-stage methods to combine IPD and AD
# for binary outcomes
# analysis at up to 1 year

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Combined IPD and AD")

source("R code/support code/binary data set-up.R")

##################################################
# summarise IPD

# remove baseline case

summ.ipd.ev <- ipd.data.ex1yr %>%
  left_join(ipd.base) %>%
  mutate(DME.1yr=ifelse(DME_FU==1 & (DMEBaseline==0 | is.na(DMEBaseline)), 1, 0),
         VH.1yr=ifelse(VH_FU==1 & (VHBaseline==0 | is.na(VHBaseline)), 1, 0)) %>%
  group_by(Trial,Class2) %>%
  summarise(across(c(DME.1yr,VH.1yr,Vitrectomy), sum))

summ.ipd.n <- ipd.data.ex1yr %>%
  group_by(Trial,Class2) %>%
  summarise(N=n())

summ.ipd <- left_join(summ.ipd.ev,summ.ipd.n)

##################################################
# DME

ag.dme <- filter(ag.data.1yr, !is.na(N.DME), Arm!="Ranibizumab + PRP", 
                 !Trial %in% c("CLARITY","PROTEUS")) %>%
  rename(N=Nevents) %>%
  select(Trial,Class2,N,N.DME)

ipd.dme <- summ.ipd %>%
  rename(N.DME=DME.1yr) %>%
  select(Trial,Class2,N,N.DME)

all.dme <- bind_rows(ipd.dme,ag.dme) %>%
  pivot_wider(id_cols=Trial,names_from=Class2,values_from=c(N,N.DME))

dme.ma <- metabin(`N.DME_Anti-VEGF`, `N_Anti-VEGF`, `N.DME_PRP`, `N_PRP`, 
                  data=all.dme,studlab=Trial)

tiff("DME forest plot - full.tiff",width=1200,height=600,compression="lzw")
forest(dme.ma,rightcols=c("effect","ci"))
dev.off()

tiff("DME forest plot - reduced.tiff",width=1200,height=600,compression="lzw")
forest(dme.ma,leftcols=c("studlab"),rightcols=c("effect","ci"))
dev.off()

###########################################
# VH

ag.vh <- filter(ag.data.1yr, !is.na(N.VH), Arm!="Ranibizumab + PRP",
                 !Trial %in% c("CLARITY","PROTEUS")) %>%
  rename(N=Nevents) %>%
  select(Trial,Class2,N,N.VH)

ipd.vh <- summ.ipd %>%
  rename(N.VH=VH.1yr) %>%
  select(Trial,Class2,N,N.VH)

all.vh <- bind_rows(ipd.vh,ag.vh) %>%
  pivot_wider(id_cols=Trial,names_from=Class2,values_from=c(N,N.VH))

vh.ma <- metabin(`N.VH_Anti-VEGF`, `N_Anti-VEGF`, `N.VH_PRP`, `N_PRP`, 
                  data=all.vh,studlab=Trial)

tiff("VH forest plot - full.tiff",width=1200,height=600,compression="lzw")
forest(vh.ma,rightcols=c("effect","ci"))
dev.off()

tiff("VH forest plot - reduced.tiff",width=1200,height=600,compression="lzw")
forest(vh.ma,leftcols=c("studlab"),rightcols=c("effect","ci"))
dev.off()

##################################################
# Vitrectomy

ag.vt <- filter(ag.data.1yr, !is.na(N.vitrectomy)) %>%
  group_by(Trial,Class2) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  rename(N=Nevents) %>%
  select(Trial,Class2,N,N.vitrectomy) %>%
  filter(Trial!="Sao Paulo A")

ipd.vt <- summ.ipd %>%
  filter(!is.na(Vitrectomy)) %>%
  rename(N.vitrectomy=Vitrectomy) %>%
  select(Trial,Class2,N,N.vitrectomy)

all.vt <- bind_rows(ipd.vt,ag.vt) %>%
  pivot_wider(id_cols=Trial,names_from=Class2,values_from=c(N,N.vitrectomy))

vt.ma <- metabin(`N.vitrectomy_Anti-VEGF`, `N_Anti-VEGF`, `N.vitrectomy_PRP`, `N_PRP`, 
                 data=all.vt,studlab=Trial)

tiff("Vitrectomy forest plot - full.tiff",width=1200,height=600,compression="lzw")
forest(vt.ma,rightcols=c("effect","ci"))
dev.off()

tiff("Vitrectomy forest plot - reduced.tiff",width=1200,height=600,compression="lzw")
forest(vt.ma,leftcols=c("studlab"),rightcols=c("effect","ci"))
dev.off()

################################################
################################################
# Other outcomes

ipd.others <- read_csv("Data/AVID IPD follow-up.csv") %>%
  filter(TimeObs>0, !is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(PatientCode,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  filter(Weeks<55) %>%
  select(PatEye,Weeks,Trial,Treatment,Arm,CST_FU,TRD_FU,NVD_FU,NVE_FU,AddTrt)

excl.trials <- c("PANORAMA","PROTOCOL W","Sao Paulo A","RECOVERY",
                 "PROTOCOL S","CLARITY","PROTEUS")

ag.others <- read_excel("Data/AD Other outcome data.xlsx") %>%
  filter(!is.na(Trial),!Trial %in% excl.trials, Arm!="Ranibizumab + PRP") %>%
  select(Trial:N,contains("NV"),N.OtherTreat)

############################################
# NVD/NVE

ipd.nv.summ <- ipd.others %>%
  select(PatEye:Arm, contains("NV")) %>%
  filter(!is.na(NVD_FU),Weeks==52) %>%
  group_by(Trial,Arm) %>%
  summarise(N=n(),
            MeanNVD=mean(NVD_FU), MeanNVD.sd=sd(NVD_FU),
            MeanNVE=mean(NVE_FU), MeanNVE.sd=sd(NVE_FU))

ag.nv <- ag.others %>%
  select(Trial,Arm,N, Weeks, contains("Mean"),contains("sd")) %>%
  filter(!is.na(MeanNVD)) %>%
  group_by(Trial,Arm) %>%
  slice_max(Weeks)

all.nv <- bind_rows(ipd.nv.summ,ag.nv) %>%
  mutate(Treat=ifelse(Arm=="PRP","PRP","Anti-VEGF")) %>%
  select(-Weeks,-Arm) %>%
  pivot_wider(id_cols=Trial,names_from=Treat,values_from=N:MeanNVE.sd)

nvd.ma <- metacont(`N_Anti-VEGF`, `MeanNVD_Anti-VEGF`,`MeanNVD.sd_Anti-VEGF`, 
                   `N_PRP`, `MeanNVD_PRP`,`MeanNVD.sd_PRP`,
                  data=filter(all.nv,Trial!="PROTEUS"),studlab=Trial)

tiff("NVD forest plot.tiff",width=1200,height=600,compression="lzw")
forest(nvd.ma,rightcols=c("effect","ci"))
dev.off()

nve.ma <- metacont(`N_Anti-VEGF`, `MeanNVE_Anti-VEGF`,`MeanNVE.sd_Anti-VEGF`, 
                   `N_PRP`, `MeanNVE_PRP`,`MeanNVE.sd_PRP`,
                   data=all.nv,studlab=Trial)

tiff("NVE forest plot.tiff",width=1200,height=600,compression="lzw")
forest(nve.ma,rightcols=c("effect","ci"))
dev.off()