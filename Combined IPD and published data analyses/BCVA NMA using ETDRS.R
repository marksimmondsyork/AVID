# NMA of bcva data - using ETDRS

library(tidyverse)
library(magrittr)
library(readxl)
library(multinma)

options(mc.cores = parallel::detectCores())

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Combined IPD and AD")

source("R code/support code/AD NMA data set-up.R")
source("R code/support code/IPD data set-up.R")

############################################################
# tidy up

ad.bcva.tb1 <- ad.bcva.tb %>%
  group_by(Trial,Arm.orig) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) %>%
  ungroup()

ad.bcva.tb2 <- ad.bcva.tb1 %>%
  filter(Weeks>44 & Weeks<60)
  
ipd.bcva.1yr <- ipd.bcva.1yr %>%
  mutate(Year=Weeks/52 - 1)

#############################################################
#############################################################
# up to 1 year

# full model

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

NMA.1yr <- combine_network(AD.setup, IPD.setup,trt_ref="PRP")

NMA.1yr.run <- nma(NMA.1yr,
              trt_effects="random",
              prior_intercept=normal(scale=10),
              prior_trt=normal(scale=10),
              prior_het=half_normal(scale=5),
              chains=4,iter=10000)


#######################################################
# class model

NMA.class.1yr <- combine_network(set_agd_arm(ad.bcva.tb1,
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

NMA.class.1yr.run <- nma(NMA.class.1yr,
                    trt_effects="random",
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_het=half_normal(scale=5))


#######################################################
# drug comparison model

NMA.drug.1yr <- combine_network(set_agd_arm(ad.bcva.tb1,
                                         study=Trial,
                                         trt=Drug,
                                         y=etdrs.mcfb,
                                         se=etdrs.mcfb.sem,
                                         sample_size=N),
                             set_ipd(ipd.bcva.1yr,
                                     study=Trial,
                                     trt=Drug,
                                     y=etdrs.cfb),
                             trt_ref="PRP")

NMA.drug.1yr.run <- nma(NMA.drug.1yr,
                    trt_effects="random",
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_het=half_normal(scale=5))

#############################################################
#############################################################
# exactly 1 year

# full model

NMA.ex1yr <- combine_network(set_agd_arm(ad.bcva.tb2,
                                         study=Trial,
                                         trt=Arm,
                                         y=etdrs.mcfb,
                                         se=etdrs.mcfb.sem,
                                         sample_size=N,
                                         trt_ref="PRP"),
                             set_ipd(ipd.bcva.1yr,
                                     study=Trial,
                                     trt=Arm,
                                     y=etdrs.cfb,
                                     trt_ref="PRP"),
                             trt_ref="PRP")

NMA.ex1yr.run <- nma(NMA.ex1yr,
                   trt_effects="random",
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10),
                   prior_het=half_normal(scale=5),
                   chains=4,iter=10000)


#######################################################
# class model

NMA.class.ex1yr <- combine_network(set_agd_arm(ad.bcva.tb2,
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

NMA.class.ex1yr.run <- nma(NMA.class.ex1yr,
                         trt_effects="random",
                         prior_intercept=normal(scale=10),
                         prior_trt=normal(scale=10),
                         prior_het=half_normal(scale=5))


#######################################################
# drug comparison model

NMA.drug.ex1yr <- combine_network(set_agd_arm(ad.bcva.tb2,
                                            study=Trial,
                                            trt=Drug,
                                            y=etdrs.mcfb,
                                            se=etdrs.mcfb.sem,
                                            sample_size=N),
                                set_ipd(ipd.bcva.1yr,
                                        study=Trial,
                                        trt=Drug,
                                        y=etdrs.cfb),
                                trt_ref="PRP")

NMA.drug.ex1yr.run <- nma(NMA.drug.ex1yr,
                        trt_effects="random",
                        prior_intercept=normal(scale=10),
                        prior_trt=normal(scale=10),
                        prior_het=half_normal(scale=5))

###########################################################
# save all results

get_summ_res <- function(nma.run){
  NMA.res <- relative_effects(nma.run,all_contrasts=TRUE)$summary %>%
    select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
  return(NMA.res)
}

get_ranks <- function(nma.run)
  print(as.matrix(posterior_rank_probs(nma.run)))

get_ranks.b <- function(nma.run)
  print(posterior_ranks(nma.run))


NMA_res_plot <- function(nma.run,filename){
plot(relative_effects(nma.run,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95)) +
    scale_x_continuous("Mean difference in ETDRS (1st named vs 2nd)", breaks=seq(-20,20,2))
  ggsave(filename=filename,device="tiff",width=6,height=4,dpi=200,compression="lzw")
}


all.models <- mget(ls(pattern=".run"))

all.res <- map(all.models,get_summ_res)
l.res <- map_int(all.res, nrow)
all.res.df <- do.call(rbind,all.res) %>%
  mutate(Model=rep(names(all.res),l.res))

write_csv(all.res.df,"NMA ETDRS unadjusted results.csv")

sink("NMA ETDRS unadjusted prob rankings.txt")
walk(all.models, get_ranks)
walk(all.models, get_ranks.b)
sink()

sink("NMA ETDRS unadjusted all results.txt")
walk(all.models, print)
sink()

filename.set <- c("NMA ETDRS 1-year main.tiff", 
                  "NMA ETDRS 1-year by class.tiff",
                  "NMA ETDRS 1-year by drug.tiff",
                  "NMA ETDRS exact 1-year main.tiff", 
                  "NMA ETDRS exact 1-year by class.tiff",
                  "NMA ETDRS exact 1-year by drug.tiff")

walk2(all.models,filename.set,NMA_res_plot)

  


###################################################################
###################################################################
# Adjusted models

# full model

NMA.adj <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(etdrs.base)),
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

NMA.adj.run <- nma(NMA.adj,
                     trt_effects="random",
                     regression= ~(etdrs.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     prior_het=half_normal(scale=5),
                     prior_aux=half_normal(scale=5) ,
                     center=T,
                     iter=5000)

#############################################################
# class model

NMA.class.adj <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(etdrs.base)),
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

NMA.class.adj.run <- nma(NMA.class.adj,
                    trt_effects="random",
                    regression= ~(etdrs.base + Year):.trtclass,
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_reg=normal(scale=2),
                    prior_het=half_normal(scale=5),
                    prior_aux=half_normal(scale=5) ,
                    center=T,
                    iter=5000)

#############################################################
# by drug model

NMA.drug.adj <- combine_network(set_agd_arm(filter(ad.bcva.tb1,!is.na(etdrs.base)),
                                            study=Trial,
                                            trt=Drug,
                                            y=etdrs.mcfb,
                                            se=etdrs.mcfb.sem,
                                            sample_size=N,
                                            trt_class=Class2),
                                set_ipd(ipd.bcva.1yr,
                                        study=Trial,
                                        trt=Drug,
                                        y=etdrs.cfb,
                                        trt_class=Class2),
                                trt_ref="PRP")

NMA.drug.adj.run <- nma(NMA.drug.adj,
                        trt_effects="random",
                        regression= ~(etdrs.base + Year):.trtclass,
                        prior_intercept=normal(scale=10),
                        prior_trt=normal(scale=10),
                        prior_reg=normal(scale=2),
                        prior_het=half_normal(scale=5),
                        prior_aux=half_normal(scale=5) ,
                        center=T,
                        iter=5000)

#########################################################################
# save results

etdrs.ref <- mean(ipd.bcva.1yr$etdrs.base[ipd.bcva.1yr$Trial=="PROTOCOL S"])
ref.data <- data.frame(etdrs.base=etdrs.ref, Year=0)

get_summ_res2 <- function(nma.run){
  NMA.res <- relative_effects(nma.run,newdata=ref.data,all_contrasts=TRUE)$summary %>%
    select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
  return(NMA.res)
}

NMA_res_plot2 <- function(nma.run,filename){
  plot(relative_effects(nma.run,newdata=ref.data,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95)) +
    scale_x_continuous("Mean difference in ETDRS (1st named vs 2nd)", breaks=seq(-20,20,2))
  ggsave(filename=filename,device="tiff",width=6,height=4,dpi=200,compression="lzw")
}

get_ranks2 <- function(nma.run)
  print(posterior_rank_probs(nma.run,newdata=ref.data))

get_ranks.b2 <- function(nma.run)
  print(posterior_ranks(nma.run,newdata=ref.data))

all.models2 <- mget(ls(pattern="adj.run"))

all.res2 <- map(all.models2,get_summ_res2)
l.res <- map_int(all.res2, nrow)
all.res.df2 <- do.call(rbind,all.res2) %>%
  mutate(Model=rep(names(all.res2),l.res))

write_csv(all.res.df2,"NMA ETDRS adjusted results.csv")

sink("NMA ETDRS adjusted prob rankings.txt")
walk(all.models2, get_ranks2)
walk(all.models2, get_ranks.b2)
sink()

sink("NMA ETDRS adjusted all results.txt")
walk(all.models2, print)
sink()

filename.set <- c("NMA ETDRS adjusted main.tiff", 
                  "NMA ETDRS adjusted by class.tiff",
                  "NMA ETDRS adjusted by drug.tiff")

walk2(all.models2,filename.set,NMA_res_plot2)
