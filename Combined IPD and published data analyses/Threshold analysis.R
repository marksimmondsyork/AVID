###############################################
# Run Threshold analyses of combined AD and IPD
###############################################

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)
library(lme4)
library(multinma)
library(nmathresh)

options(mc.cores = parallel::detectCores())

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/combined IPD and AD")

source("R code/support code/threshold support functions.R")
source("R code/support code/run threshold analysis code.R")
source("R code/support code/AD NMA data set-up.R")
source("R code/support code/IPD data set-up.R")

fct2num <- function(x) {as.numeric(levels(x))}

############################################################
# tidy up

ad.bcva.1yr <- ad.bcva.tb %>%
  group_by(Trial,Arm.orig) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) %>%
  ungroup()

ad.bcva.adj <- ad.bcva.1yr %>%
  filter(!is.na(logMAR.base))

ipd.bcva.1yr <- ipd.bcva.1yr %>%
  mutate(Year=Weeks/52 - 1)

######################################################
# add numeric Trial/Treat

ad.bcva.1yr$Trialn <- as.numeric(factor(ad.bcva.1yr$Trial)) + 3
ad.bcva.1yr$Trtn <- as.numeric(fct_infreq(ad.bcva.1yr$Arm,ordered=T))

ad.bcva.adj$Trialn <- as.numeric(factor(ad.bcva.adj$Trial)) + 3
ad.bcva.adj$Trtn <- as.numeric(fct_infreq(ad.bcva.adj$Arm,ordered=T))
ad.bcva.adj$Drugn <- as.numeric(fct_infreq(ad.bcva.adj$Drug,ordered=T))
ad.bcva.adj$Classn <- case_when(
  ad.bcva.adj$Class=="PRP" ~ 1,
  ad.bcva.adj$Class=="Anti-VEGF + PRP" ~ 3,
  ad.bcva.adj$Class=="Anti-VEGF" ~ 2)

ipd.bcva.1yr$Trialn <- as.numeric(factor(ipd.bcva.1yr$Trial))
ipd.bcva.1yr <- ipd.bcva.1yr %>%
  mutate(Trtn=case_when(
    Arm=="PRP" ~ 1,
    Arm=="Aflibercept" ~ 6,
    Arm=="Ranibizumab" ~ 4,
    Arm=="Ranibizumab + PRP" ~ 5),
    Trtn.adj=case_when(
      Arm=="PRP" ~ 1,
      Arm=="Aflibercept" ~ 5,
      Arm=="Ranibizumab" ~ 4,
      Arm=="Ranibizumab + PRP" ~ 3),
    Drugn=case_when(
      Drug=="PRP" ~ 1,
      Drug=="Aflibercept" ~ 4,
      Drug=="Ranibizumab" ~ 2),
    Classn=case_when(
      Class=="PRP" ~ 1,
      Class=="Anti-VEGF + PRP" ~ 3,
      Class=="Anti-VEGF" ~ 2))

trt.names.1yr <- c(levels(fct_infreq(ad.bcva.1yr$Arm)),"Aflibercept")
trt.names.adj <- c(levels(fct_infreq(ad.bcva.adj$Arm)),"Aflibercept")
dr.names <- c(levels(fct_infreq(ad.bcva.adj$Drug)),"Aflibercept")
cl.names <- c("PRP","Anti-VEGF","Anti-VEGF + PRP")

#################################################
# data structures

comb.data.1yr <- ipd.bcva.1yr %>%
  select(Trialn,Trtn) %>%
  group_by(Trialn,Trtn) %>%
  slice_head(n=1) %>%
  bind_rows(select(ad.bcva.1yr,Trialn,Trtn)) %>%
  arrange(Trialn,Trtn) %>%
  select(-PatEye)

comb.data.adj <- ipd.bcva.1yr %>%
  select(Trialn,Trtn,Drugn,Classn) %>%
  group_by(Trialn,Trtn) %>%
  slice_head(n=1) %>%
  bind_rows(select(ad.bcva.adj,Trialn,Trtn,Drugn,Classn)) %>%
  arrange(Trialn,Trtn) %>%
  select(-PatEye) 

########################################################
# Unadjusted 1 year

AD.setup <- set_agd_arm(ad.bcva.1yr,
                        study=Trialn,
                        trt=Trtn,
                        y=logMAR.mcfb,
                        se=logMAR.mcfb.sem,
                        sample_size=N)

IPD.setup <- set_ipd(ipd.bcva.1yr,
                     study=Trialn,
                     trt=Trtn,
                     y=logMAR.cfb)

NMA.1yr <- combine_network(AD.setup, IPD.setup,trt_ref="1")

NMA.1yr.run <- nma(NMA.1yr,
                   trt_effects="random",
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10),
                   prior_het=half_normal(scale=5),
                   chains=4,iter=10000)

thresh.res.1yr <- threshold_run(NMA.1yr.run,comb.data.1yr,trt.names,
                                "Threshold plot 1yr unadjusted.png")


##################################################
# Analysis adjusting for time/baseline

NMA.adj <- combine_network(set_agd_arm(ad.bcva.adj,
                                       study=Trialn,
                                       trt=Trtn,
                                       y=logMAR.mcfb,
                                       se=logMAR.mcfb.sem,
                                       sample_size=N,
                                       trt_class=Classn),
                           set_ipd(ipd.bcva.1yr,
                                   study=Trialn,
                                   trt=Trtn,
                                   y=logMAR.cfb,
                                   trt_class=Classn),
                           trt_ref="1")

NMA.adj.run <- nma(NMA.adj,
                   trt_effects="random",
                   regression= ~(logMAR.base + Year):.trtclass,
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10),
                   prior_reg=normal(scale=2),
                   prior_het=half_normal(scale=5),
                   prior_aux=half_normal(scale=5) ,
                   center=T,
                   iter=5000)

logMAR.ref <- mean(ipd.bcva.1yr$logMAR.base[ipd.bcva.1yr$Trial=="PROTOCOL S"])
ref.data <- data.frame(logMAR.base=logMAR.ref, Year=0)

thresh.res.adj <- threshold_run(NMA.adj.run,comb.data.adj,trt.names,
                                "Threshold plot adjusted.png",
                                adjusted=T,ref.data=ref.data)


#############################################################
# class model

NMA.class.adj <- combine_network(set_agd_arm(ad.bcva.adj,
                                             study=Trialn,
                                             trt=Classn,
                                             y=logMAR.mcfb,
                                             se=logMAR.mcfb.sem,
                                             sample_size=N,
                                             trt_class=Class2),
                                 set_ipd(ipd.bcva.1yr,
                                         study=Trialn,
                                         trt=Classn,
                                         y=logMAR.cfb,
                                         trt_class=Class2),
                                 trt_ref="1")

NMA.class.adj.run <- nma(NMA.class.adj,
                         trt_effects="random",
                         regression= ~(logMAR.base + Year):.trtclass,
                         prior_intercept=normal(scale=10),
                         prior_trt=normal(scale=10),
                         prior_reg=normal(scale=2),
                         prior_het=half_normal(scale=5),
                         prior_aux=half_normal(scale=5) ,
                         center=T,
                         iter=5000)

class.data <- comb.data.adj %>%
  ungroup() %>%
  select(Trialn,Classn) %>%
  rename(Trtn=Classn) 

thresh.res.adj.class <- threshold_run(NMA.class.adj.run,class.data,cl.names,
                                "Threshold plot adjusted by class.png",
                                adjusted=T,ref.data=ref.data)


#############################################################
# by drug model

NMA.drug.adj <- combine_network(set_agd_arm(ad.bcva.adj,
                                            study=Trialn,
                                            trt=Drugn,
                                            y=logMAR.mcfb,
                                            se=logMAR.mcfb.sem,
                                            sample_size=N,
                                            trt_class=Class2),
                                set_ipd(ipd.bcva.1yr,
                                        study=Trialn,
                                        trt=Drugn,
                                        y=logMAR.cfb,
                                        trt_class=Class2),
                                trt_ref="1")

NMA.drug.adj.run <- nma(NMA.drug.adj,
                        trt_effects="random",
                        regression= ~(logMAR.base + Year):.trtclass,
                        prior_intercept=normal(scale=10),
                        prior_trt=normal(scale=10),
                        prior_reg=normal(scale=2),
                        prior_het=half_normal(scale=5),
                        prior_aux=half_normal(scale=5) ,
                        center=T,
                        iter=5000)

drug.data <- comb.data.adj %>%
  ungroup() %>%
  select(Trialn,Drugn) %>%
  rename(Trtn=Drugn)

thresh.res.adj.drug <- threshold_run(NMA.drug.adj.run,drug.data,dr.names,
                                "Threshold plot adjusted by drug.png",
                                adjusted=T,ref.data=ref.data)

