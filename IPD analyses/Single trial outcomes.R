# Analysis of outcomes in single trials

library(tidyverse)
library(magrittr)
library(lme4)

setwd ("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/IPD")
options(dplyr.summarise.inform = FALSE)

###################################################
# data

ipd.bin <- read_csv("Data/AVID IPD binary.csv")
ipd.fu.1 <- read_csv("Data/AVID IPD follow-up.csv")
ipd.base <- read_csv("Data/AVID IPD baseline.csv")

# edit data
ipd.fu <- ipd.fu.1 %>%
  select(Trial,PatientID:TimeObs,
         CST_FU:Employ_FU, TRD_FU:AddTrt, NVD_FU:Retinopathy_FU, RIP:OccDiscomf) %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) 

ipd.base <- ipd.base %>%
  filter(!is.na(Treatment)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) 

################################################
# amount of data
data_count <- function(var){
ipd.count <- tibble(trial=ipd.fu$Trial, 
                    var=ipd.fu[[var]],
                    ID=ipd.fu$PatEye) %>%
  group_by(trial,ID) %>%
  summarise(perc.var1=100*length(var[!is.na(var)])/length(var)) %>%
  ungroup() %>%
  group_by(trial) %>%
  summarise(perc.var=100*length(perc.var1[perc.var1>0])/length(perc.var1))

return(ipd.count)
}

outcome.names <- names(select(ipd.fu,-c(Trial:TimeObs,PatEye)))

data.summ <- map_dfr(outcome.names, data_count) %>%
  mutate(outcome=rep(outcome.names,each=3))

oc.onetrial <- data.summ %>%
  group_by(outcome) %>%
  summarise(max.perc=max(perc.var),
            n.tr=length(perc.var[perc.var>0])) %>%
  filter(max.perc>5, n.tr==1)


#############################################################
# DRSS

drss.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, DRSS_FU) %>%
  left_join(select(ipd.base,c(PatEye,DRSSBaseline))) %>%
  filter(!is.na(DRSS_FU),DRSS_FU<11,DRSSBaseline<11, TimeObs %in% c(seq(12,60,12))) %>%
  mutate(DRSS.change=DRSSBaseline-DRSS_FU,
         Arm=ifelse(Treatment==1,"Ranibizumab","PRP"))

drss.data %$% t.test(DRSS_FU~Treatment)
drss.data %$% t.test(DRSS.change~Treatment)
drss.data %$% wilcox.test(DRSS_FU~Treatment)
drss.data %$% wilcox.test(DRSS.change~Treatment)

ggplot(drss.data,aes(x=DRSS_FU,fill=Arm)) +
  geom_bar() +
  scale_x_continuous("DRSS at 1 year",breaks=1:10)
ggsave("DRSS 1 year PrS.tiff",width=6,height=4,compression="lzw")

ggplot(drss.data,aes(x=DRSS.change,fill=Arm)) +
  geom_bar() +
  scale_x_continuous("DRSS improvement at 1 year",breaks=-10:10)
ggsave("DRSS improvement 1 year PrS.tiff",width=6,height=4,compression="lzw")


################################################################
# Driving/Employment/Reading

read.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, Reading_FU) %>%
  left_join(select(ipd.base,c(PatEye,ReadingBaseline))) %>%
  filter(!is.na(Reading_FU) & Reading_FU<9 & ReadingBaseline<9) %>%
  filter(TimeObs %in% c(seq(12,60,12))) %>%
  mutate(Read.change=ReadingBaseline-Reading_FU)

read.summ <- read.data %>%
  group_by(TimeObs,Treatment,Read.change) %>%
  summarise(Npats=n())

drive.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, Driving_FU) %>%
  left_join(select(ipd.base,c(PatEye,DrivingBaseline))) %>%
  filter(!is.na(Driving_FU) & Driving_FU<9 & DrivingBaseline<9) %>%
  filter(TimeObs %in% c(seq(12,60,12))) %>%
  mutate(Drive.change=Driving_FU-DrivingBaseline)

drive.summ <- drive.data %>%
  group_by(TimeObs,Treatment,Drive.change) %>%
  summarise(Npats=n())

employ.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, Employ_FU) %>%
  left_join(select(ipd.base,c(PatEye,EmployBaseline))) %>%
  filter(!is.na(Employ_FU) & Employ_FU<9 & EmployBaseline<9) %>%
  filter(TimeObs %in% c(seq(12,60,12))) %>%
  mutate(Employ.change=EmployBaseline-Employ_FU)

employ.summ <- employ.data %>%
  group_by(TimeObs,Treatment,Employ.change) %>%
  summarise(Npats=n())

########################################################
# QoL

qol.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs,EQ5D_FU,NEI_FU) %>%
  left_join(select(ipd.base,c(PatEye,EQ5DBase,NEIBase))) %>%
  mutate(across(EQ5D_FU:NEI_FU, ~ifelse(.x>100, NA, .x))) %>%
  filter(!(is.na(EQ5D_FU) & is.na(NEI_FU))) %>%
  filter(TimeObs==12) %>%
  mutate(EQ5D.change=EQ5D_FU-EQ5DBase,
         NEI.change=NEI_FU-NEIBase,
         Arm=ifelse(Treatment==1,"Aflibercept","PRP"))

qol.data %$% wilcox.test(EQ5D_FU~Treatment)
qol.data %$% wilcox.test(EQ5D.change~Treatment)

qol.data %$% wilcox.test(NEI_FU~Treatment)
qol.data %$% wilcox.test(NEI.change~Treatment)

qol.data %$% summary(lm(EQ5D_FU ~ Treatment))
qol.data %$% summary(lm(EQ5D.change ~ Treatment))
qol.data %$% summary(lm(NEI_FU ~ Treatment))
qol.data %$% summary(lm(NEI.change ~ Treatment))

ggplot(qol.data,aes(x=EQ5D.change,fill=Arm)) +
  geom_bar() +
  scale_x_continuous("Change in EQ5D after 1 year",breaks=-10:10)
ggsave("EQ5D 1 year Clarity.tiff",width=6,height=4,compression="lzw")

ggplot(qol.data,aes(x=NEI.change,fill=Arm)) +
  geom_bar() +
  scale_x_continuous("Change in NEI after 1 year",breaks=seq(-30,30,5))
ggsave("NEI 1 year Clarity.tiff",width=6,height=4,compression="lzw")

####################################################
# NVD/NVE

nvd.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs,NVD_FU,NVE_FU) %>%
  filter(!(is.na(NVD_FU) & is.na(NVE_FU)))

filter(nvd.data, TimeObs==12) %$% summary(lm(NVD_FU ~ Treatment))
filter(nvd.data, TimeObs==12) %$% summary(lm(NVE_FU ~ Treatment))

###################################################
# AEs

# Death (no data)
death.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, Death) %>%
  filter(Death==1)

# Conj. Heam. 
ch.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, ConjHeam) %>%
  filter(ConjHeam==1)

# Ocular discomfort
occdisc.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs,OccDiscomf) %>%
  filter(OccDiscomf==1)

# retinal neovasc. (no data)
rn.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs,RetinalNeov) %>%
  filter(!is.na(RetinalNeov))

# stroke
stroke.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs,Stroke) %>%
  filter(Stroke==1)

################################################
# Retinopathy grade (coding in PROTEUS?)

retin.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, Retinopathy_FU) %>%
  left_join(select(ipd.base,c(PatEye,RetinopathyBaseline))) %>%
  filter(!is.na(Retinopathy_FU))

filter(retin.data, TimeObs==12) %$% summary(glm(Retinopathy_FU ~ Treatment,family=binomial))

###############################################
# VA (only in AVEGF arm)

va.data <- ipd.fu %>%
  select(PatEye,Trial:TimeObs, VA_FU) %>%
  filter(!is.na(VA_FU))


############################################################
############################################################
# Modifiers of BCVA (DRSS)

# BCVA data
ipd.base.bcva <- ipd.base %>%
  select(Trial,PatEye,Arm,ETDRSBaseline,BCVABaseline) %>%
  rename(etdrs.base=ETDRSBaseline,
         logMAR.base=BCVABaseline) 

ipd.bcva <- ipd.fu.1 %>%
  select(Trial,PatientID:TimeObs,ETDRS_FU,BCVA_FU) %>% 
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) %>%
  filter(TimeObs>0) %>%
  left_join(ipd.base.bcva) %>%
  mutate(etdrs.cfb=ETDRS_FU-etdrs.base,
         logMAR.cfb=BCVA_FU-logMAR.base) %>%
  rename(etdrs=ETDRS_FU,
         logMAR=BCVA_FU) %>%
  filter(!is.na(logMAR.cfb))

base.drss <- ipd.base %>%
  select(PatEye,DRSSBaseline)

drss.bcva <- left_join(ipd.bcva,base.drss) %>%
  filter(!is.na(DRSSBaseline))

# regression
# 1 year

reg1y <- filter(drss.bcva, TimeObs==12) %$% lm(etdrs.cfb ~ factor(Treatment)*DRSSBaseline)
summary(reg1y)

reg2y <- filter(drss.bcva, TimeObs==24) %$% lm(etdrs.cfb ~ factor(Treatment)*DRSSBaseline)
summary(reg2y)

reg5y <- filter(drss.bcva, TimeObs==60) %$% lm(etdrs.cfb ~ factor(Treatment)*DRSSBaseline)
summary(reg5y)
