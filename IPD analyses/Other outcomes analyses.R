# TEST IPD regression analyses

library(tidyverse)
library(magrittr)
library(lme4)

setwd ("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/IPD")

source("R code/support code/IPD Analysis functions.R")

# data set-up
source("R code/support code/non-BCVA data set-up.R")

##########################################
# DME and VH


dmevh.data <- dmevh.data.ex1yr %>%
  left_join(ipd.base.2) %>%
  mutate(DME.1yr=ifelse(DME_FU==1 & (DMEBaseline==0 | is.na(DMEBaseline)), 1, 0),
         VH.1yr=ifelse(VH_FU==1 & (VHBaseline==0 | is.na(VHBaseline)), 1, 0)) 

#outcome.list <- ipd.ex1yr %$% list(DME.1yr,RIP,SAE,TRD_FU,VH.1yr)
outcome.list <- dmevh.data %$% list(DME.1yr,VH.1yr)

dmevh.1y <- map(outcome.list, 
                    ~dmevh.data %$% one_stage_MA(.x,Trial,Treatment,"RR","RE_both"))

############################################################
# meta-regressions

covars.list <- dmevh.data %$% list(Age.c,Sex.c,BCVABaseline.c,Diabetes,
                               PrevAVEGF,PrevPRP,VHBaseline,
                               DMEBaseline,HbA1c.c,CSTBaseline.c)

# for one outcome
p_mreg <- function(outcome,trial,treatment,covar){
  tryCatch(one_stage_MA_inter(outcome,trial,treatment,covar,"OR","RE_both"),
           error=function(e) NULL)
}

mreg_one <- function(outcome,covars){
  out <- map(covars, 
             ~dmevh.data %$% p_mreg(outcome,Trial,Treatment,.x))
  return(out)
}

# all models
dmevh.1y.mreg <- map(outcome.list, ~mreg_one(.x,covars.list))

#####################################################
# CST
ipd.1y.cst <- ipd.other %>%
  filter(TimeObs==12) %>%
  select(PatEye:TimeObs,contains("CST")) %>%
  left_join(ipd.base.2) %>%
  mutate(CST.cfb=CST_FU-CSTBaseline)

cst.1y.cfb <- ipd.1y.cst %$% one_stage_MA(CST.cfb,Trial,Treatment,"MD","RE_both")
cst.1y.anova <- ipd.1y.cst %$% one_stage_MA(CST_FU,Trial,Treatment,"anova","RE_both",CSTBaseline)

# meta-regressions
covars.list <- ipd.1y.cst %$% list(Age.c,Sex.c,BCVABaseline.c,Diabetes,
                                   PrevAVEGF,PrevPRP,VHBaseline,
                                   DMEBaseline,HbA1c.c,CSTBaseline.c)

cst.1y.mreg <- map(covars.list,
                   ~ipd.1y.cst %$% 
                     one_stage_MA_inter(CST.cfb,Trial,Treatment,.x,"MD","RE_both"))

#########################################
# other binary

ipd.1yr.other <- ipd.other %>%
  filter(TimeObs==12) %>%
  select(-CST_FU) %>%
  mutate(across(RIP:TRD_FU, ~replace_na(.x,0)))

outcome.list <- ipd.1yr.other %$% list(RIP,SAE,TRD_FU)

other.1y <- map(outcome.list, 
                ~ipd.1yr.other %$% one_stage_MA(.x,Trial,Treatment,"RR","RE_both"))

##########################################
# extract results
onames1 <- c("CST","DME","VH","RIP","SAE","TRD")
onames2 <- c("CST","DME","VH")
covnames <- c("Age","Sex","BCVABaseline","Diabetes",
              "PrevAVEGF","PrevPRP","VHBaseline","DMEBaseline",
              "HbA1c","CSTBaseline")

# main analyses 
summ.cst.1y <- cont_convert(summary(cst.1y.cfb))
summ.dmevh.1y <- map_dfr(dmevh.1y, ~log_convert(summary(.x)))
summ.other.1y <- map_dfr(other.1y, ~log_convert(summary(.x)))

all.res.main <- rbind(summ.cst.1y,summ.dmevh.1y,summ.other.1y) %>%
  mutate(Outcome=rep(onames1,each=2),
         Param=rep(c("Int","OR"),6))

# m-regs
summ.mreg.cst.1y <- map_dfr(cst.1y.mreg, ~cont_convert(summary(.x))) 

# function to extract results from m-regs
res_extract <- function(res){
  if(!is.null(res))
    out <- log_convert(summary(res))
  if(is.null(res))
    out <- data.frame(Estimate=rep(NA,4),CIlow=rep(NA,4),CIhigh=rep(NA,4))
  return(out)
}

summ.mreg.dmevh.1y <- map_dfr(dmevh.1y.mreg, 
                            ~map_dfr(.x, res_extract)) 

all.mreg.1yr <- bind_rows(summ.mreg.cst.1y,summ.mreg.dmevh.1y) %>%
  mutate(Outcome=rep(onames2,each=10*4),
         Covariate=rep(rep(covnames,each=4),3),
         Parameter=rep(c("Intercept","Trt effect","Covar","Interaction"),10*3))

# save output
write_csv(all.res.main,"IPD non-BCVA outcomes 1 yr results.csv")
write_csv(all.mreg.1yr,"IPD non-BCVA outcomes 1 yr mreg results.csv")




##########################################################
##########################################################
# repeated measures analyses (3 monthly intervals)

rep_measures_model <- function(outcome,treat,trial,time,IDVar,oc="bin"){
  if(oc=="cont")
    out <- lmer(outcome ~ treat*time + (1|IDVar) + (1+treat|trial))
  if(oc=="bin")
    out <- glmer(outcome ~ treat*time + (1|IDVar) + (1+treat|trial),family=binomial)
  return(out)
}

ipd.1yr <- ipd.ma %>%
  filter(TimeObs %in% c(3,7,9,12)) %>%
  mutate(Year=TimeObs/12-1)

ipd.2yr <- ipd.ma %>%
  filter(TimeObs %in% seq(3,24,3)) %>%
  mutate(Year=TimeObs/12-1)
  
ipd.5yr <- ipd.ma %>%
   filter(TimeObs %in% c(3,7,9,12*(1:5))) %>%
  mutate(Year=TimeObs/12-1)

outcome.list.1yr <- ipd.1yr %$% list(DME_FU,RIP,SAE,TRD_FU,VH_FU)
outcome.list.2yr <- ipd.2yr %$% list(DME_FU,RIP,SAE,TRD_FU,VH_FU)
outcome.list.5yr <- ipd.5yr %$% list(DME_FU,RIP,SAE,TRD_FU,VH_FU)

cst.1yr <- ipd.1yr %$% rep_measures_MA(CST_FU,Treatment,Trial,TimeObs,PatEye,"linear","MD")
bin.1yr <- map(outcome.list.1yr, 
               ~ipd.1yr %$% rep_measures_MA(.x,Treatment,Trial,TimeObs,PatEye,"linear","OR"))

cst.2yr <- ipd.2yr %$% rep_measures_MA(CST_FU,Treatment,Trial,TimeObs,PatEye,"linear","MD")
bin.2yr <- map(outcome.list.2yr, 
               ~ipd.2yr %$% rep_measures_MA(.x,Treatment,Trial,TimeObs,PatEye,"linear","OR"))

cst.5yr <- ipd.5yr %$% rep_measures_MA(CST_FU,Treatment,Trial,TimeObs,PatEye,"linear","MD")
bin.5yr <- map(outcome.list.5yr, 
               ~ipd.5yr %$% rep_measures_MA(.x,Treatment,Trial,TimeObs,PatEye,"linear","OR"))

# extract results

onames <- c("CST","DME","RIP","SAE","TRD","VH")

sum.cst.1yr <- cont_convert(summary(cst.1yr))
sum.cst.2yr <- cont_convert(summary(cst.2yr))
sum.cst.5yr <- cont_convert(summary(cst.5yr))

sum.bin.1yr <- map_dfr(bin.1yr, ~log_convert(summary(.x)))
sum.bin.2yr <- map_dfr(bin.2yr, ~log_convert(summary(.x)))
sum.bin.5yr <- map_dfr(bin.5yr, ~log_convert(summary(.x)))

all.rep.res <- bind_rows(sum.cst.1yr,sum.bin.1yr,
                         sum.cst.2yr,sum.bin.2yr,
                         sum.cst.5yr,sum.bin.5yr) %>%
  mutate(Outcome=rep(rep(onames,each=4),3),
         Parameter=rep(c("Intercept","Trt effect","Time","Interaction"),6*3),
         Year=rep(c(1,2,5),each=4*6))


write_csv(all.rep.res,"IPD non-BCVA outcomes repeated measures results.csv")
