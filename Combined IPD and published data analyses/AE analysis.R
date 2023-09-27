# set up AE / further outcomes analysis

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/combined IPD and AD")

###########################################################
# AE data from publications
AE.ag.data <- read_excel("data/AD AE data.xlsx") %>%
  filter(!is.na(Trial), Weeks>0) %>%
  select(-contains("Mean"),-contains("driving")) %>%
  mutate(Class=ifelse(Arm=="PRP","PRP","Anti-VEGF")) %>%
  arrange(Trial,Weeks,Class)

# trials to drop from AD set
ipd.trials <- c("CLARITY","PROTOCOL S","PROTEUS")
npdr.trials <- c("PROTOCOL W","PANORAMA")
excl.trials <- c(npdr.trials,"RECOVERY","Sao Paulo A") 

AE.ag.data.1 <-  AE.ag.data %>%
  filter(Weeks>0, !Trial %in% excl.trials)  %>%
  group_by(Trial,Class) %>%
  slice_tail(n=1) %>%
  select(-Vit.Heam,-Macular.oedema,-Vitreous.floaters) %>%
  rename(RIP=Raised.intraoc.press, RetinalNeov=Retinal.neovasuc,
         ConjHeam=Conjunctival.haem,OccDiscomf=Ocular.pain)

#######################################################
# IPD on AEs
AE.ipd <- read_csv("Data/AVID IPD follow-up.csv") %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  filter(Weeks<60) %>%
  mutate(Class=ifelse(Arm=="PRP","PRP","Anti-VEGF")) %>%
  select(PatEye:Class,Trial,SAE:AddTrt,NVD_FU:Retinopathy_FU,
         RIP,Cataracts,RetinalNeov:OccDiscomf)

AE.ipd.summ <- AE.ipd %>%
  select(-NVD_FU,-NVE_FU,-OtherAE) %>%
  mutate(across(SAE:OccDiscomf, ~replace_na(.x, 0))) %>%
  group_by(PatEye,Trial,Class) %>%
  summarise(across(SAE:OccDiscomf, ~ifelse(sum(.x)>0, 1, 0))) %>%
  ungroup() %>%
  group_by(Trial,Class) %>%
  summarise(across(SAE:OccDiscomf, sum)) 

IPD.n <- AE.ipd %>%
  group_by(PatEye,Trial,Class) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) %>%
  group_by(Trial,Class) %>%
  summarise(N=n()) 

AE.ipd.1 <- left_join(IPD.n,AE.ipd.summ)
  
all.AE <- full_join(AE.ipd.1,AE.ag.data.1)

############################################################
# convert to MA style

AE.cont <- all.AE %>%
  filter(Class=="PRP") %>%
  select(-Arm,-Weeks) %>%
  rename_with(~paste(.,".c",sep=""),N:CVD.death) %>%
  ungroup() %>%
  select(-Trial,-Class)

AE.all <- filter(all.AE, Class!="PRP") %>%
  select(-Arm,-Weeks) %>%
  bind_cols(AE.cont)

MA_binary <- function(outcome){
    data1 <- AE.all %>%
    select(Trial,Class,N,N.c,contains(outcome))
  colnames(data1)[5:6] <- c("Ev.trt","Ev.c")
  data1 <- filter(data1,!is.na(Ev.trt) & !is.na(Ev.c)) %>%
    filter(Ev.trt>0 | Ev.c>0) %>%
    group_by(Trial,Class) %>%
    slice_head(n=1)
  ma1 <- metabin(Ev.trt,N,Ev.c,N.c,data=data1,sm="RR",studlab=Trial)
  return(ma1)
}

# run
onames <- colnames(AE.all[4:19])
all.MAs <- map(onames, MA_binary)
names(all.MAs) <- onames


# save full forest plots (no MA)
forest_func <- function(outcome){
  MAi <- all.MAs[[outcome]]
  png(paste(outcome,"forest plot.png",sep=" "),width=1000)
  forest(MAi,fixed=F,random=F,
         rightcols=c("effect","ci"))
  dev.off()
}

walk(onames,forest_func)

# MAs

extract_MA_results <- function(ma.res){
  out <- data.frame(logRR=ma.res$TE.random, 
                    selogRR=ma.res$seTE.random,
                    RR=exp(ma.res$TE.random),
                    RR.low=exp(ma.res$lower.random),
                    RR.high=exp(ma.res$upper.random),
                    I2=ma.res$I2,
                    Ntrials=ma.res$k,
                    Nevents=sum(ma.res$event.e+ma.res$event.c))
  return(out)
}

onames.full <- c("SAE","Death","Myocardial infarct.","Stroke","Additional treatment",
                 "Retinopathy","Raised intraoc. pressure","Cataracts",
                 "Retinal neovascularisation","Conjuctival haemorrhage",
                 "Ocular pain","Retinal detachment","Visual disturbances","Macular fibrosis",
                 "Retinal tear","Death due to CVD")

AE.MA.res.summ <- map_dfr(all.MAs, extract_MA_results)%>%
  mutate(Outcome=onames.full) %>%
  arrange(Outcome)

# plot
pdata <- AE.MA.res.summ %>%
  filter(Ntrials>1)

tiff("Adverse event MA forest plot.tiff",width=1000, height=700)
forest(metagen(logRR,selogRR,sm="RR",data=pdata),
       studlab=Outcome,fixed=F,random=F,
       leftcols=c("studlab","Ntrials","Nevents"),
       leftlabs=c("Outcome","N. trials","N. events"),
       xlab="Relative risk",smlab="Relative risk",
       print.I2=F,print.tau2=F,print.pval.Q=F,
       label.right="Favours PRP",
       label.left="Favours Anti-VEGF")
dev.off()
