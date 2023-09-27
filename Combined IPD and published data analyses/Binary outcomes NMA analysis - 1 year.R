# use NMA methods to combine IPD and AD
# for binary outcomes
# analysis at up to 1 year

library(tidyverse)
library(magrittr)
library(readxl)
library(multinma)

options(mc.cores = parallel::detectCores())

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Combined IPD and AD")

source("R code/support code/binary data set-up.R")

#################################################
# ipd remove baseline cases

ipd.data.a <- ipd.data.ex1yr %>%
  left_join(ipd.base) %>%
  mutate(DME.1yr=ifelse(DME_FU==1 & (DMEBaseline==0 | is.na(DMEBaseline)), 1, 0),
         VH.1yr=ifelse(VH_FU==1 & VHBaseline==0, 1, 0))

##################################################
# DME

ag.dme.data <- filter(ag.data.1yr, !is.na(N.DME), !Trial %in% c("CLARITY","PROTEUS")) 
ipd.dme.data <- filter(ipd.data.a,!is.na(DME.1yr))

dme.net.d <- combine_network(
  set_agd_arm(ag.dme.data,
              study=Trial,
              trt=Arm,
              r=N.DME,
              n=Nevents),
  set_ipd(ipd.dme.data,
          study=Trial,
          trt=Arm,
          r=DME.1yr),
  trt_ref="PRP")

dme.nma.d <- nma(dme.net.d,
                 trt_effects="random",
                 prior_intercept=normal(scale=10),
                 prior_trt=normal(scale=10),
                 prior_het=half_normal(scale=5),
                 iter=5000)

# combining anti-VEGFs

dme.net.t <- combine_network(
  set_agd_arm(ag.dme.data,
              study=Trial,
              trt=Class2,
              r=N.DME,
              n=Nevents),
  set_ipd(ipd.dme.data,
          study=Trial,
          trt=Class2,
          r=DME.1yr),
  trt_ref="PRP")

dme.nma.t <- nma(dme.net.t,
                 trt_effects="random",
                 prior_intercept=normal(scale=10),
                 prior_trt=normal(scale=10),
                 prior_het=half_normal(scale=5),
                 iter=5000)

# time adjusted

dme.nma.adj <- nma(dme.net.t,
                 regression= ~Weeks:.trt,
                 trt_effects="random",
                 prior_intercept=normal(scale=10),
                 prior_trt=normal(scale=10),
                 prior_het=half_normal(scale=5),
                 prior_reg=normal(scale=2),
                 prior_aux=half_normal(scale=5) ,
                 center=F,
                 iter=5000)



###########################################
# VH

ag.vh.data <- filter(ag.data.1yr, !is.na(N.VH), !Trial %in% c("CLARITY","PROTEUS"))
ipd.vh.data <- filter(ipd.data.a,!is.na(VH.1yr))

# At 1 year by drug
vh.net.d <- combine_network(
  set_agd_arm(ag.vh.data,
                     study=Trial,
                     trt=Arm,
                     r=N.VH,
                     n=Nevents),
  set_ipd(ipd.vh.data,
          study=Trial,
          trt=Arm,
          r=VH.1yr),
  trt_ref="PRP")

vh.nma.d <- nma(vh.net.d,
                  trt_effects="random",
                  prior_intercept=normal(scale=10),
                  prior_trt=normal(scale=10),
                  prior_het=half_normal(scale=5),
                  iter=5000)

# At 1 year combining all anti-VEGFs

vh.net.t <- combine_network(
  set_agd_arm(ag.vh.data,
              study=Trial,
              trt=Class2,
              r=N.VH,
              n=Nevents),
  set_ipd(ipd.vh.data,
          study=Trial,
          trt=Class2,
          r=VH.1yr),
  trt_ref="PRP")

vh.nma.t <- nma(vh.net.t,
                trt_effects="random",
                prior_intercept=normal(scale=10),
                prior_trt=normal(scale=10),
                prior_het=half_normal(scale=5),
                iter=5000)

# time adjusted

vh.nma.adj <- nma(vh.net.t,
                   regression= ~Weeks:.trt,
                   trt_effects="random",
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10),
                   prior_het=half_normal(scale=5),
                   prior_reg=normal(scale=2),
                   prior_aux=half_normal(scale=5) ,
                   center=F,
                   iter=5000)

###################################################
# save

ref.vals <- data.frame(Weeks=52)

dme.d.res <- relative_effects(dme.nma.d,all_contrasts=TRUE)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
dme.t.res <- relative_effects(dme.nma.t)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
dme.adj.res <- relative_effects(dme.nma.adj,newdata=ref.vals)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)

vh.d.res <- relative_effects(vh.nma.d,all_contrasts=TRUE)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
vh.t.res <- relative_effects(vh.nma.t)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)
vh.adj.res <- relative_effects(vh.nma.adj,newdata=ref.vals)$summary %>%
  select(`.trtb`,`.trta`,`mean`,`sd`,`2.5%`,`97.5%`)

all.res <- rbind(dme.d.res,dme.t.res,dme.adj.res,
                 vh.d.res,vh.t.res,vh.adj.res) %>%
  mutate(OR=exp(mean),OR.low=exp(`2.5%`),OR.high=exp(`97.5%`)) %>%
  mutate(Outcome=rep(c("DME","VH"),each=10),
         Model=rep(c("Drug","Combined","Adjusted","Drug","Combined","Adjusted"),c(6,1,1,10,1,1)))

write_csv(all.res,"Binary outcome NMA results 2.csv")

plot(relative_effects(dme.nma.d,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA DME 1-year by drug.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")

plot(relative_effects(dme.nma.t), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA DME 1-year main.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")

plot(relative_effects(dme.nma.adj,newdata=ref.vals), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA DME 1-year adjusted.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")

plot(relative_effects(vh.nma.d,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA VH 1-year by drug.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")

plot(relative_effects(vh.nma.t), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA VH 1-year main.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")

plot(relative_effects(vh.nma.adj,newdata=ref.vals), ref_line=0, .width=c(0,0.95)) +
  scale_x_continuous("Odds ratio (1st named vs 2nd)")
ggsave(filename="NMA VH 1-year adjusted.tiff",device="tiff",
       width=6,height=4,dpi=200,compression="lzw")


all.models <- mget(ls(pattern=".nma"))
sink("NMA binary outcomes all results.txt")
walk(all.models, print)
sink()




