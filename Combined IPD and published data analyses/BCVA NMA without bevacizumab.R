# NMA of bcva data - removing bevacizumab

library(tidyverse)
library(magrittr)
library(readxl)
library(multinma)

options(mc.cores = parallel::detectCores())

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Combined IPD and AD")

source("R code/support code/AD NMA data set-up.R")
source("R code/support code/IPD data set-up.R")

############################################################
# remove bevacizumab
ad.bcva.tb1 <- ad.bcva.tb %>%
  filter(!Trial %in% c("Ahmad","Ali","Marashi","Rebecca","Roohipoor")) %>%
  group_by(Trial,Arm.orig) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) %>%
  ungroup()

ipd.bcva.1yr <- ipd.bcva.1yr %>%
  mutate(Year=Weeks/52 - 1)

#############################################################

AD.setup <- set_agd_arm(ad.bcva.tb1,
                        study=Trial,
                        trt=Arm,
                        y=logMAR.mcfb,
                        se=logMAR.mcfb.sem,
                        sample_size=N,
                        trt_ref="PRP")

IPD.setup <- set_ipd(ipd.bcva.1yr,
                     study=Trial,
                     trt=Arm,
                     y=logMAR.cfb,
                     trt_ref="PRP"
)

NMA.setup <- combine_network(AD.setup, IPD.setup,trt_ref="PRP")

#plot(NMA.setup)

NMA.re <- nma(NMA.setup,
              trt_effects="random",
              prior_intercept=normal(scale=10),
              prior_trt=normal(scale=10),
              prior_het=half_normal(scale=5),
              chains=4,iter=10000)

NMA.res <- relative_effects(NMA.re,all_contrasts=TRUE)$summary
#plot(relative_effects(NMA.re,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95))

#######################################################
# class model

NMA.class <- combine_network(set_agd_arm(ad.bcva.tb1,
                                         study=Trial,
                                         trt=Class,
                                         y=logMAR.mcfb,
                                         se=logMAR.mcfb.sem,
                                         sample_size=N),
                             set_ipd(ipd.bcva.1yr,
                                     study=Trial,
                                     trt=Class,
                                     y=logMAR.cfb),
                             trt_ref="PRP")

#plot(NMA.class)

NMA.class.re <- nma(NMA.class,
                    trt_effects="random",
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_het=half_normal(scale=5))

NMA.class.res <- relative_effects(NMA.class.re,all_contrasts=TRUE)$summary
#plot(relative_effects(NMA.class.re,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95))

#####################################################################
# time/baseline adjusted

NMA.tb1 <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(logMAR.base)),
                                      study=Trial,
                                      trt=Arm,
                                      y=logMAR.mcfb,
                                      se=logMAR.mcfb.sem,
                                      sample_size=N,
                                      trt_class=Class),
                          set_ipd(ipd.bcva.1yr,
                                  study=Trial,
                                  trt=Arm,
                                  y=logMAR.cfb,
                                  trt_class=Class),
                          trt_ref="PRP")

NMA.time.base <- nma(NMA.tb1,
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

NMA.timebase.res <- relative_effects(NMA.time.base,newdata=ref.data,all_contrasts=TRUE)$summary
#plot(NMA.time.base, ref_line=0, .width=c(0,0.95))

#############################################################
# class model adjusted for time/baseline

NMA.tb2 <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(logMAR.base)),
                                       study=Trial,
                                       trt=Class,
                                       y=logMAR.mcfb,
                                       se=logMAR.mcfb.sem,
                                       sample_size=N,
                                       trt_class=Class2),
                           set_ipd(ipd.bcva.1yr,
                                   study=Trial,
                                   trt=Class,
                                   y=logMAR.cfb,
                                   trt_class=Class2),
                           trt_ref="PRP")

NMA.tb.class <- nma(NMA.tb2,
                     trt_effects="random",
                     regression= ~(logMAR.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     prior_het=half_normal(scale=5),
                     prior_aux=half_normal(scale=5) ,
                     center=T,
                     iter=5000)

NMA.tb.class.res <- relative_effects(NMA.tb.class,newdata=ref.data,all_contrasts=TRUE)$summary
#plot(NMA.tb.class, ref_line=0, .width=c(0,0.95))


###########################################################
###########################################################
# REPEATED USING ETDRS

AD.setup <- set_agd_arm(ad.bcva.tb1,
                        study=Trial,
                        trt=Arm,
                        y=etdrs.mcfb,
                        se=etdrs.mcfb.sem,
                        sample_size=N,
                        trt_ref="PRP")

IPD.setup <- set_ipd(ipd.bcva.1yr,
                     study=Trial,
                     trt=Arm,
                     y=etdrs.cfb,
                     trt_ref="PRP"
)

NMA.ETDRS.setup <- combine_network(AD.setup, IPD.setup,trt_ref="PRP")

#plot(NMA.ETDRS.setup)

NMA.ETDRS.re <- nma(NMA.ETDRS.setup,
              trt_effects="random",
              prior_intercept=normal(scale=10),
              prior_trt=normal(scale=10),
              prior_het=half_normal(scale=5),
              chains=4,iter=10000)

NMA.ETDRS.res <- relative_effects(NMA.ETDRS.re,all_contrasts=TRUE)$summary
#plot(relative_effects(NMA.ETDRS.re,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95))

#######################################################
# class model

NMA.ETDRS.class <- combine_network(set_agd_arm(ad.bcva.tb1,
                                         study=Trial,
                                         trt=Class,
                                         y=etdrs.mcfb,
                                         se=etdrs.mcfb.sem,
                                         sample_size=N),
                             set_ipd(ipd.bcva.1yr,
                                     study=Trial,
                                     trt=Class,
                                     y=etdrs.cfb),
                             trt_ref="PRP")

#plot(NMA.ETDRS.class)

NMA.ETDRS.class.re <- nma(NMA.ETDRS.class,
                    trt_effects="random",
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_het=half_normal(scale=5))

NMA.ETDRS.class.res <- relative_effects(NMA.ETDRS.class.re,all_contrasts=TRUE)$summary
#plot(relative_effects(NMA.ETDRS.class.re,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95))

#####################################################################
# time/baseline adjusted

NMA.ETDRS.tb1 <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(etdrs.base)),
                                       study=Trial,
                                       trt=Arm,
                                       y=etdrs.mcfb,
                                       se=etdrs.mcfb.sem,
                                       sample_size=N,
                                       trt_class=Class),
                           set_ipd(ipd.bcva.1yr,
                                   study=Trial,
                                   trt=Arm,
                                   y=etdrs.cfb,
                                   trt_class=Class),
                           trt_ref="PRP")

NMA.ETDRS.time.base <- nma(NMA.ETDRS.tb1,
                     trt_effects="random",
                     regression= ~(etdrs.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     prior_het=half_normal(scale=5),
                     prior_aux=half_normal(scale=5) ,
                     center=T,
                     iter=5000)

etdrs.ref <- mean(ipd.bcva.1yr$etdrs.base[ipd.bcva.1yr$Trial=="PROTOCOL S"])
ref.data <- data.frame(etdrs.base=etdrs.ref, Year=0)

NMA.ETDRS.timebase.res <- relative_effects(NMA.ETDRS.time.base,newdata=ref.data,all_contrasts=TRUE)$summary
#plot(NMA.ETDRS.timebase.res, ref_line=0, .width=c(0,0.95))

#############################################################
# class model adjusted for time/baseline

NMA.ETDRS.tb2 <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(etdrs.base)),
                                       study=Trial,
                                       trt=Class,
                                       y=etdrs.mcfb,
                                       se=etdrs.mcfb.sem,
                                       sample_size=N,
                                       trt_class=Class2),
                           set_ipd(ipd.bcva.1yr,
                                   study=Trial,
                                   trt=Class,
                                   y=etdrs.cfb,
                                   trt_class=Class2),
                           trt_ref="PRP")

NMA.ETDRS.tb.class <- nma(NMA.ETDRS.tb2,
                    trt_effects="random",
                    regression= ~(etdrs.base + Year):.trtclass,
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_reg=normal(scale=2),
                    prior_het=half_normal(scale=5),
                    prior_aux=half_normal(scale=5) ,
                    center=T,
                    iter=5000)

NMA.ETDRS.tb.class.res <- relative_effects(NMA.ETDRS.tb.class,newdata=ref.data,all_contrasts=TRUE)$summary
#plot(NMA.ETDRS.tb.class.res, ref_line=0, .width=c(0,0.95))

###########################################################
# save all results

all.res <- mget(ls(pattern=".res"))
l.res <- map_int(all.res, nrow)

sel_cols <- function(res){
  res %>% select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
}

all.res.df <- map_dfr(all.res,sel_cols) %>%
  mutate(Model=rep(names(all.res),l.res))
  
write_csv(all.res.df,"NMA results without bevacizumab.csv")
