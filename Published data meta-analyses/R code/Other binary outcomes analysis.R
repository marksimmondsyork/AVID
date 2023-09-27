# binary outcomes analysis

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

# data
outcomes <- read_excel("data/Other Outcome data.xlsx") %>%
  filter(!is.na(Trial), Weeks>0) %>%
  select(-N,-contains("NV")) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "Control",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab"
  )) %>%
  arrange(Trial,Weeks,Arm)

# Control arm data
outcomes.cont <- outcomes %>%
  filter(Drug=="Control") %>%
  select(-Drug,-Arm) %>%
  rename_with(~paste(.,".c",sep=""),starts_with("N"))

outcomes.all <- filter(outcomes, Drug!="Control") %>%
  left_join(outcomes.cont)

# MA function - all trials and times
MA_binary <- function(outcome){
  
  data1 <- outcomes.all %>%
    select(Trial,Arm,Drug,Weeks,Nevents,Nevents.c,contains(outcome))
  colnames(data1)[7:8] <- c("Ev.trt","Ev.c")
  data1 <- filter(data1,!is.na(Ev.trt) & !is.na(Ev.c))
  ma1 <- metabin(Ev.trt,Nevents,Ev.c,Nevents.c,data=data1,sm="RR",studlab=Trial,subgroup=Drug)
  return(ma1)
}

# run
onames <- colnames(outcomes.all[5:11])
all.MAs <- map(onames, MA_binary)
names(all.MAs) <- onames

# save full forest plots (no MA)
forest_func <- function(outcome){
  MAi <- all.MAs[[outcome]]
  png(paste(outcome,"forest plot.png",sep=" "),width=1000)
  forest(MAi,fixed=F,random=F,
         leftcols=c("studlab","Arm","Weeks"),leftlabs=c("Trial","Interv.","Weeks"),
         hetstat=F,label.left="Favours anti-VEGF",label.right="Favours PRP")
  dev.off()
}

walk(onames,forest_func)

######################################################
# meta-analysis

# reduce data

outcomes.reduced <- outcomes.all %>%
  filter(Weeks<200) %>%
  group_by(Trial) %>%
  slice_max(order_by=Weeks) %>%
  slice_head(n=1)

outcomes.PDR <- outcomes.reduced %>%
  filter(!Trial %in% c("PANORAMA","PROTOCOL W"))

onames <- colnames(outcomes.reduced[c(5,8,9,10)])
onames.full <- c("Regression of Neovasc.","DMO","Vitreous Haemorrhage","Vitrectomy")

# ma function

MA_binary_2 <- function(outcome){
  
  data1 <- outcomes.PDR %>%
    select(Trial,Arm,Drug,Weeks,Nevents,Nevents.c,contains(outcome))
  colnames(data1)[7:8] <- c("Ev.trt","Ev.c")
  data1 <- filter(data1,!is.na(Ev.trt) & !is.na(Ev.c))
  ma1 <- metabin(Ev.trt,Nevents,Ev.c,Nevents.c,data=data1,sm="RR",studlab=Trial,subgroup=Drug)
  return(ma1)
}

extract_MA_results <- function(ma.res){
  out <- data.frame(logRR=ma.res$TE.random, 
                    selogRR=ma.res$seTE.random,
                    RR=exp(ma.res$TE.random),
                    RR.low=exp(ma.res$lower.random),
                    RR.high=exp(ma.res$upper.random),
                    I2=ma.res$I2,
                    Ntrials=ma.res$k,
                    Npats=sum(ma.res$n.e+ma.res$n.c),
                    Nevents=sum(ma.res$event.e+ma.res$event.c))
  return(out)
}

MA.res <- map(onames, MA_binary_2)

MA.res.summ <- map_dfr(MA.res, extract_MA_results) %>%
  mutate(Outcome=onames.full) %>%
  arrange(Outcome)

# plot
pdata <- filter(MA.res.summ,Ntrials>1)

tiff("Other outcomes MA forest plot.tiff",width=1000, height=700, compression="lzw")
forest(metagen(logRR,selogRR,sm="RR",data=pdata),
       studlab=Outcome,fixed=F,random=F,
       leftcols=c("studlab","Ntrials","Npats","Nevents"),
       leftlabs=c("Outcome","N. trials","N. patients","N. events"),
       xlab="Relative risk",smlab="Relative risk",
       hetstat=F,label.left="Favours anti-VEGF",label.right="Favours PRP")
dev.off()


