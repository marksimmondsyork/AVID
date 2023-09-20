# NVD/NVE analysis

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

# NVD and NVE data from publications
nvd.nve <- read_excel("data/Other Outcome data.xlsx") %>%
  filter(!is.na(Trial)) %>%
  select(Trial:N, contains("NV")) %>%
  filter(!is.na(MeanNVE) | !is.na(MeanNVD)) %>%
  filter(Trial!="RECOVERY") %>%
  arrange(Trial,Weeks,Arm)

#####################################################
# edit to calculate change from baseline
nvd.base <- nvd.nve %>%
  filter(Weeks==0) %>%
  select(-N,-Weeks,-contains("Reduc")) %>%
  rename(NVE.base=MeanNVE, NVE.base.sd=MeanNVE.sd,
         NVD.base=MeanNVD, NVD.base.sd=MeanNVD.sd) %>%
  arrange(Trial,Arm)

nvd.nve2 <- left_join(filter(nvd.nve,Weeks>0), nvd.base) %>%
  mutate(NVE.cfb=ifelse(!is.na(NVE.base), MeanNVE-NVE.base, MeanNVE),
         NVD.cfb=ifelse(!is.na(NVD.base), MeanNVD-NVD.base, MeanNVD),
         NVE.cfb.sd=ifelse(!is.na(NVE.base), sqrt(MeanNVE.sd^2+NVE.base.sd^2), MeanNVE.sd),
         NVD.cfb.sd=ifelse(!is.na(NVD.base), sqrt(MeanNVD.sd^2+NVD.base.sd^2), MeanNVD.sd)) %>%
  mutate(NVD.cfb.se=NVD.cfb.sd / sqrt(N), NVE.cfb.se=NVE.cfb.sd / sqrt(N)) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "Control",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab"
  ))

# differences between arms
nvd.cont <- nvd.nve2 %>%
  filter(Drug=="Control") %>%
  select(Trial,Weeks,N,contains("cfb")) %>%
  rename(Nc=N,
         NVE.cfb.c=NVE.cfb, NVE.cfb.c.sd=NVE.cfb.sd, NVE.cfb.c.se=NVE.cfb.se,
         NVD.cfb.c=NVD.cfb, NVD.cfb.c.sd=NVD.cfb.sd, NVD.cfb.c.se=NVD.cfb.se)

nvd.nve.cfb <- left_join(filter(nvd.nve2,Drug!="Control"), nvd.cont)

# up to 13 weeks
nvd.nve.13 <- nvd.nve.cfb %>%
  select(Trial,Arm,Drug,Weeks,N,Nc,contains("NV")) %>%
  filter(Weeks<15) %>%
  group_by(Trial) %>%
  slice_max(order_by=Weeks) %>%
  slice_head(n=1)

# max. follow-up
nvd.nve.max <- nvd.nve.cfb %>%
  select(Trial,Arm,Drug,Weeks,N,Nc,contains("NV")) %>%
  group_by(Trial) %>%
  slice_max(order_by=Weeks) %>%
  slice_head(n=1)

##############################################
# NVE forest plots and MAs

ma.nve <- metacont(N,NVE.cfb,NVE.cfb.sd,Nc,NVE.cfb.c,NVE.cfb.c.sd, 
                   data=filter(nvd.nve.cfb,!is.na(NVE.cfb)), 
                   studlab=Trial,subgroup=Drug)

png("NVE all times forest plot.png",width=1000,height=600)
forest(ma.nve,fixed=F,random=F,
       leftcols=c("studlab","Arm","Weeks"),
       leftlabs=c("Trial","Interv.","Weeks"))
dev.off()

ma.nve.13 <- metacont(N,NVE.cfb,NVE.cfb.sd,Nc,NVE.cfb.c,NVE.cfb.c.sd, 
                   data=nvd.nve.13 ,studlab=Trial,subgroup=Drug)

png("NVE 13weeks forest plot.png",width=1000,height=600)
forest(ma.nve.13,fixed=F,random=T,
       leftcols=c("studlab","Arm"),leftlabs=c("Trial","Interv."))
dev.off()

ma.nve.max <- metacont(N,NVE.cfb,NVE.cfb.sd,Nc,NVE.cfb.c,NVE.cfb.c.sd, 
                      data=nvd.nve.max ,studlab=Trial,subgroup=Drug)

png("NVE max time forest plot.png",width=1000,height=600)
forest(ma.nve.max,fixed=F,random=T,
       leftcols=c("studlab","Arm"),leftlabs=c("Trial","Interv."))
dev.off()


##################################################
# NVD 

ma.nvd <- metacont(N,NVD.cfb,NVD.cfb.sd,Nc,NVD.cfb.c,NVD.cfb.c.sd, 
                   data=filter(nvd.nve.cfb,!is.na(NVD.cfb)), 
                   studlab=Trial,subgroup=Drug)

png("NVD all times forest plot.png",width=1000,height=600)
forest(ma.nvd,fixed=F,random=F,
       leftcols=c("studlab","Arm","Weeks"),
       leftlabs=c("Trial","Interv.","Weeks"))
dev.off()

ma.nvd.13 <- metacont(N,NVD.cfb,NVD.cfb.sd,Nc,NVD.cfb.c,NVD.cfb.c.sd, 
                      data=nvd.nve.13 ,studlab=Trial)

png("NVD 13weeks forest plot.png",width=1000,height=600)
forest(ma.nvd.13,fixed=F,random=T,
       leftcols=c("studlab","Arm"),leftlabs=c("Trial","Interv."))
dev.off()

ma.nvd.max <- metacont(N,NVD.cfb,NVD.cfb.sd,Nc,NVD.cfb.c,NVD.cfb.c.sd, 
                      data=nvd.nve.max ,studlab=Trial)

png("NVD max time forest plot.png",width=1000,height=600)
forest(ma.nvd.max,fixed=F,random=T,
       leftcols=c("studlab","Arm"),leftlabs=c("Trial","Interv."))
dev.off()
