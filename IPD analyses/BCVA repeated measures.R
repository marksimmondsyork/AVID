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

# repeated measures (linear time effect)
repmeas.lmar.1yr <- ipd.1yr %$%
  rep_measures_MA(logMAR.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

repmeas.lmar.2yr <- ipd.2yr %$%
  rep_measures_MA(logMAR.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

repmeas.lmar.5yr <- ipd.5yr %$%
  rep_measures_MA(logMAR.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

repmeas.etdrs.1yr <- ipd.1yr %$%
  rep_measures_MA(etdrs.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

repmeas.etdrs.2yr <- ipd.2yr %$%
  rep_measures_MA(etdrs.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

repmeas.etdrs.5yr <- ipd.5yr %$%
  rep_measures_MA(etdrs.cfb,Treatment,Trial,Year,PatEye,
                  model="linear",measure="MD")

# summarise results
res.repmeas <- map_dfr(list(repmeas.lmar.1yr,repmeas.lmar.2yr,repmeas.lmar.5yr,
                            repmeas.etdrs.1yr,repmeas.etdrs.2yr,repmeas.etdrs.5yr),
                       ~cont_convert(summary(.x))) %>%
  mutate(Parameter=rep(c("Intercept","Treatment","Time","Interaction"),6),
         Model="Rep. measures",
         Outcome=rep(c("LogMAR","ETDRS"),each=12),
         Time=rep(c("1 year","2 years","5 years","1 year","2 years","5 years"),each=4))


write_csv(res.repmeas,"results/IPD repeated measures models B.csv")


###################################################
# rep measures with meta-regression

repmeas_interaction <- function(outcome,treat,trial,time,covar,IDVar){
  res <- lmer(outcome ~  treat*time*covar - treat:time:covar + (1|IDVar/trial) + (1+treat|trial))
}

ipd.1yr.2 <- left_join(ipd.1yr,ipd.base.2)
ipd.2yr.2 <- left_join(ipd.2yr,ipd.base.2)
ipd.5yr.2 <- left_join(ipd.5yr,ipd.base.2)

covars.1y <- ipd.1yr.2 %$% list(Age.c,Sex.c,BCVABaseline.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

covars.2y <- ipd.2yr.2 %$% list(Age.c,Sex.c,BCVABaseline.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

covars.5y <- ipd.5yr.2 %$% list(Age.c,Sex.c,BCVABaseline.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

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

# summary

covnames <- c("Age","Sex","BCVABaseline","Diabetes","PrevAVEGF","PrevPRP",
              "VHBaseline","DMEBaseline","HbA1c","CSTBaseline")

all.res <- list(mreg.etdrs.1y,mreg.etdrs.2y,mreg.etdrs.5y,
                mreg.lmar.1y,mreg.lmar.2y,mreg.lmar.5y)

all.res.conv <- map_dfr(all.res, ~map_dfr(., ~cont_convert(summary(.x))))

summ.mreg.rm <- all.res.conv %>%
  mutate(Parameter=rep(c("Intercept","Trt effect","Time","Covar","TrtxTime","TrtxCovar","TimexCovar"),60),
         Covariate=rep(rep(covnames,each=7),6),
         Time=rep(rep(c("1 year","2 years","5 years"),each=70), 2),
         Outcome=rep(c("ETDRS","logMAR"),each=210))

write_csv(summ.mreg.rm,"results/IPD repeated measures regression models C -no 3-way.csv")
