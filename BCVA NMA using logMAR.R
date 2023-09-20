# NMA of BCVA data

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

source("R code/support code/NMA data set-up.R")
source("R code/support code/run NMA code.R")

#################################################
# without NPDR trials

bcva.1yr <- filter(bcva.1yr, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.2yr <- filter(bcva.2yr, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.max <- filter(bcva.max, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.tb <- filter(bcva.tb, !Trial %in% c("PROTOCOL W","PANORAMA"))

# sizes
bcva.tb.size <- bcva.tb %>%
  select(Trial,Arm,N) %>%
  group_by(Trial) %>%
  summarise(Npats=sum(N))

##################################################
# 1 year

NMA.1yr <- set_agd_arm(bcva.1yr,
                       study=Trial,
                       trt=Arm2,
                       y=logMAR.mcfb,
                       se=logMAR.mcfb.sem,
                       sample_size=N,
                       trt_ref="PRP")

NMA.1yr.res <- run_NMA(NMA.1yr,"BCVA logMAR NMA 1 year")

##################################################
# 1 - 2 years

NMA.2yr <- set_agd_arm(bcva.2yr,
                       study=Trial,
                       trt=Arm2,
                       y=logMAR.mcfb,
                       se=logMAR.mcfb.sem,
                       sample_size=N,
                       trt_ref="PRP")

NMA.2yr.res <- run_NMA(NMA.2yr,"BCVA logMAR NMA 1-2 years")

##################################################
# Max FU (up to 2 years)

NMA.maxFU <- set_agd_arm(bcva.max,
                         study=Trial,
                         trt=Arm2,
                         y=logMAR.mcfb,
                         se=logMAR.mcfb.sem,
                         sample_size=N,
                         trt_ref="PRP")

NMA.maxFU.res <- run_NMA(NMA.maxFU,"BCVA logMAR NMA max FU")


###########################################
#  drug-based model 
# max FU only

NMA.reduced <- set_agd_arm(bcva.max,
                       study=Trial,
                       trt=Drug,
                       y=logMAR.mcfb,
                       se=logMAR.mcfb.sem,
                       sample_size=N,
                       trt_ref="PRP")

NMA.reduced.res <- run_NMA(NMA.reduced,"BCVA logMAR NMA by drug")

###########################################################
# class based model
# max FU only

NMA.class <- set_agd_arm(bcva.max,
                         study=Trial,
                         trt=Class,
                         y=logMAR.mcfb,
                         se=logMAR.mcfb.sem,
                         sample_size=N,
                         trt_ref="PRP")

NMA.class.res <- run_NMA(NMA.class,"BCVA logMAR NMA by class")
