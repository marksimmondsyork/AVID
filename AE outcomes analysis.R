# AE outcomes analysis

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

# AE data from publications
AE.data <- read_excel("data/AE data.xlsx") %>%
  filter(!is.na(Trial), Weeks>0) %>%
  select(-contains("Mean"),-contains("driving")) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "Control",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab"
  )) %>%
  arrange(Trial,Weeks,Arm)

# Control arm data
AE.cont <- AE.data %>%
  filter(Drug=="Control") %>%
  select(-Drug,-Arm) %>%
  rename_with(~paste(.,".c",sep=""),N:SAE)

AE.all <- filter(AE.data, Drug!="Control") %>%
  left_join(AE.cont)

MA_binary <- function(outcome){
  
  data1 <- AE.all %>%
    select(Trial,Arm,Drug,Weeks,N,N.c,contains(outcome))
  colnames(data1)[7:8] <- c("Ev.trt","Ev.c")
  data1 <- filter(data1,!is.na(Ev.trt) & !is.na(Ev.c))
  ma1 <- metabin(Ev.trt,N,Ev.c,N.c,data=data1,sm="RR",studlab=Trial)
  return(ma1)
}

# run
onames <- colnames(AE.all[5:21])
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

AE.reduced <- AE.all %>%
  group_by(Trial) %>%
  slice_head(n=1)

AE.PDR <- AE.reduced %>%
  filter(!Trial %in% c("PANORAMA","PROTOCOL W","RECOVERY"))

onames <- colnames(AE.reduced[5:21])[-c(9,11)]
onames.full <- c("Raised intraocular pressure","Vitreous haemorrhage","Retinal detachment",
                 "Cataracts","Macular oedema","Visual disturbances","Ocular pain","Retinal neovasc.",
                 "Retinal tear","Conj. haemorrhage",
                 "Death","Cardiovascular death","Myocardial Infarct.","Stroke","SAE")

MA_binary_2 <- function(outcome){
  
  data1 <- AE.PDR %>%
    select(Trial,Arm,Drug,Weeks,N,N.c,contains(outcome))
  colnames(data1)[7:8] <- c("Ev.trt","Ev.c")
  data1 <- filter(data1,!is.na(Ev.trt) & !is.na(Ev.c))
  ma1 <- metabin(Ev.trt,N,Ev.c,N.c,data=data1,sm="RR",studlab=Trial)
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
                    Nevents=sum(ma.res$event.e+ma.res$event.c))
  return(out)
}

AE.MA.PDR.res <- map(onames, MA_binary_2)

AE.MA.PDR.summ <- map_dfr(AE.MA.PDR.res, extract_MA_results) %>%
  mutate(Outcome=onames.full) %>%
  arrange(Outcome)

# plot
pdata <- filter(AE.MA.PDR.summ) %>%
  filter(Ntrials>1, !grepl("Macu",Outcome), !grepl("Vitr",Outcome))

tiff("Adverse event MA forest plot.tiff",width=1000, height=700,compression="lzw")
forest(metagen(logRR,selogRR,sm="RR",data=pdata),
       studlab=Outcome,fixed=F,random=F,
       leftcols=c("studlab","Ntrials","Nevents"),
       leftlabs=c("Outcome","N. trials","N. events"),
       xlab="Relative risk",smlab="Relative risk",
       hetstat=F,label.right="Favours PRP",label.left="Favours Anti-VEGF")
dev.off()



