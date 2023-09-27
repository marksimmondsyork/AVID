# IPD regression analyses

library(tidyverse)
library(magrittr)
library(lme4)

setwd ("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/IPD")

source("R code/support code/IPD data BCVA set-up.R")
source("R code/support code/IPD Analysis functions.R")

#######################################################
# BCVA at 1 yr (meta-analysis)

bcva.1y <- ipd.bcva %>%
  filter(TimeObs==12) %>%
  left_join(ipd.base.2)

# treatment only (assume common trt effect)
lmar.1yr.anova <- bcva.1y %$%
  one_stage_cont_MA(logMAR.base,logMAR,Trial,Treatment,"ancova","MD","RE_both")
lmar.1yr.cfb <- bcva.1y %$%
  one_stage_cont_MA(logMAR.base,logMAR,Trial,Treatment,"change","MD","RE_both")

etdrs.1yr.anova <- bcva.1y %$%
  one_stage_cont_MA(etdrs.base,etdrs,Trial,Treatment,"ancova","MD","RE_both")
etdrs.1yr.cfb <- bcva.1y %$%
  one_stage_cont_MA(etdrs.base,etdrs,Trial,Treatment,"change","MD","RE_both")

# summarise results
res.1yr <- map_dfr(list(lmar.1yr.anova,lmar.1yr.cfb,etdrs.1yr.anova,etdrs.1yr.cfb),
                   ~cont_convert(summary(.x))) %>%
  mutate(Parameter=rep(c("Intercept","Baseline","Treatment","Intercept","Treatment"),2),
         Model=rep(c("Anova","CFB","Anova","CFB"),c(3,2,3,2)),
         Outcome=rep(c("LogMAR","ETDRS"),each=5))

write_csv(res.1yr,"results/IPD 1 year treatment effect models.csv")


###################################################################
# regress against covariates
covars.list <- bcva.1y %$% list(Age.c,Sex.c,etdrs.base.c,logMAR.base.c,Diabetes,
                                PrevAVEGF,PrevPRP,VHBaseline,
                                DMEBaseline,HbA1c.c,CSTBaseline.c)

lmar.1y.mreg <- map(covars.list,
                       ~bcva.1y %$% one_stage_MA_inter(logMAR.cfb,Trial,Treatment,.x,"MD","RE_both"))

etdrs.1y.mreg <- map(covars.list,
                    ~bcva.1y %$% one_stage_MA_inter(etdrs.cfb,Trial,Treatment,.x,"MD","RE_both"))


covnames <- c("Age","Sex","Base ETDRS","Base logMAR","Diabetes","PrevAVEGF","PrevPRP",
              "VHBaseline","DMEBaseline","HbA1c","CSTBaseline")

# summary
summ.mreg.lmar.1y <- map_dfr(lmar.1y.mreg, ~cont_convert(summary(.x))) %>%
  mutate(Parameter=rep(c("Intercept","Trt effect","Covar","Interaction"),11),
         Covariate=rep(covnames,each=4),
         Time="1 year - no RM",
         Outcome="LogMAR")

summ.mreg.etdrs.1y <- map_dfr(etdrs.1y.mreg, ~cont_convert(summary(.x))) %>%
  mutate(Parameter=rep(c("Intercept","Trt effect","Covar","Interaction"),11),
         Covariate=rep(covnames,each=4),
         Time="1 year - no RM",
         Outcome="ETDRS")

mreg.results <- bind_rows(summ.mreg.etdrs.1y,summ.mreg.lmar.1y) %>%
  arrange(Outcome,Time)

write_csv(mreg.results,"results/IPD 1 year meta-regression models.csv")

##################################################
# additional checks

mreg.2.1yr <- ipd.1yr.2 %$% 
  lmer(logMAR.cfb ~  Treatment*(Year+BCVABaseline.c+Sex.c+VHBaseline) + (1|PatEye/Trial) + (1+Treatment|Trial))

mreg.etdrs.1yr <- ipd.1yr.2 %$% 
  lmer(etdrs.cfb ~  Treatment*(Year+etdrs.base) + (1|PatEye/Trial) + (1+Treatment|Trial))

