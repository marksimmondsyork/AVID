library(tidyverse)
library(magrittr)
library(lme4)

setwd ("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/IPD")

source("R code/support code/IPD data BCVA set-up.R")
source("R code/support code/IPD Analysis functions.R")

####################################################
# repeated measures and time analyses
# 1, 2 years and 5 years follow-up

ipd.1yr <- ipd.bcva %>%
  filter(TimeObs %in% c(3,7,9,12)) %>%
  mutate(Year=TimeObs/12-1)

ipd.2yr <- ipd.bcva %>%
  filter(TimeObs %in% c(3,7,9,12,15,18,21,24)) %>%
  mutate(Year=TimeObs/12-1)

ipd.5yr <- ipd.bcva %>%
  filter(TimeObs %in% c(3,7,9,12*(1:5))) %>%
  mutate(Year=TimeObs/12-1)

ipd.1yr.2 <- left_join(ipd.1yr,ipd.base.2)
ipd.2yr.2 <- left_join(ipd.2yr,ipd.base.2)
ipd.5yr.2 <- left_join(ipd.5yr,ipd.base.2)

covars.1y <- ipd.1yr.2 %$% list(Age.c,Sex.c,etdrs.base.c,logMAR.base.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

covars.2y <- ipd.2yr.2 %$% list(Age.c,Sex.c,etdrs.base.c,logMAR.base.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

covars.5y <- ipd.5yr.2 %$% list(Age.c,Sex.c,etdrs.base.c,logMAR.base.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

###################################################
# rep measures with meta-regression


repmeas_interaction <- function(outcome,treat,trial,time,covar,IDVar){
  res <- lmer(outcome ~  factor(treat)*(time+covar) + (1|IDVar/trial) + (1+treat|trial))
}

repmeas_interaction_2way <- function(outcome,treat,trial,time,covar,IDVar){
  res <- lmer(outcome ~  factor(treat)*time*covar - factor(treat):time:covar + (1|IDVar/trial) + (1+treat|trial))
}

##################################################
# models without 2-way interactions

mreg.etdrs.1y <- map(covars.1y,
                     ~ipd.1yr.2 %$%
                       repmeas_interaction(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.etdrs.2y <- map(covars.2y,
                     ~ipd.2yr.2 %$%
                       repmeas_interaction(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.etdrs.5y <- map(covars.5y,
                     ~ipd.5yr.2 %$%
                       repmeas_interaction(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.1y <- map(covars.1y,
                    ~ipd.1yr.2 %$%
                      repmeas_interaction(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.2y <- map(covars.2y,
                    ~ipd.2yr.2 %$%
                      repmeas_interaction(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.5y <- map(covars.5y,
                    ~ipd.5yr.2 %$%
                      repmeas_interaction(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))


##################################################
# models with 2-way interactions

mreg.etdrs.1y.2w <- map(covars.1y,
                     ~ipd.1yr.2 %$%
                       repmeas_interaction_2way(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.etdrs.2y.2w <- map(covars.2y,
                     ~ipd.2yr.2 %$%
                       repmeas_interaction_2way(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.etdrs.5y.2w <- map(covars.5y,
                     ~ipd.5yr.2 %$%
                       repmeas_interaction_2way(etdrs.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.1y.2w <- map(covars.1y,
                    ~ipd.1yr.2 %$%
                      repmeas_interaction_2way(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.2y.2w <- map(covars.2y,
                    ~ipd.2yr.2 %$%
                      repmeas_interaction_2way(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))

mreg.lmar.5y.2w <- map(covars.5y,
                    ~ipd.5yr.2 %$%
                      repmeas_interaction_2way(logMAR.cfb,Treatment,Trial,Year,.x,PatEye))

##########################################################
# summary

covnames <- c("Age","Sex","Base ETDRS","Base logMAR","Diabetes","PrevAVEGF","PrevPRP",
              "VHBaseline","DMEBaseline","HbA1c","CSTBaseline")

all.res <- list(mreg.etdrs.1y,mreg.etdrs.2y,mreg.etdrs.5y,
                mreg.lmar.1y,mreg.lmar.2y,mreg.lmar.5y)

all.res.2w <- list(mreg.etdrs.1y.2w,mreg.etdrs.2y.2w,mreg.etdrs.5y.2w,
                mreg.lmar.1y.2w,mreg.lmar.2y.2w,mreg.lmar.5y.2w)

all.res.conv <- map_dfr(all.res, ~map_dfr(., ~cont_convert(summary(.x))))

all.res.conv.2w <- map_dfr(all.res.2w, ~map_dfr(., ~cont_convert(summary(.x))))

summ.mreg.rm <- all.res.conv %>%
  mutate(Parameter=rep(c("Intercept","Trt effect","Time","Covar","TrtxTime","TrtxCovar"),66),
         Covariate=rep(rep(covnames,each=6),6),
         Time=rep(rep(c("1 year","2 years","5 years"),each=66), 2),
         Outcome=rep(c("ETDRS","logMAR"),each=198))

summ.mreg.rm.2w <- all.res.conv.2w %>%
  mutate(Parameter=rep(c("Intercept","Trt effect","Time","Covar","TrtxTime","TrtxCovar","TimexCovar"),66),
         Covariate=rep(rep(covnames,each=7),6),
         Time=rep(rep(c("1 year","2 years","5 years"),each=77), 2),
         Outcome=rep(c("ETDRS","logMAR"),each=231))


write_csv(summ.mreg.rm,"results/IPD repeated measures regression models A -no 2-way.csv")

write_csv(summ.mreg.rm.2w,"results/IPD repeated measures regression models B -no 3-way.csv")
