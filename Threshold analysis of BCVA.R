###############################################
# Run Threshold analyses using PDR trials only
###############################################

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

source("R code/support code/threshold support functions.R")
source("R code/support code/run threshold analysis code.R")
source("R code/support code/NMA data set-up.R")

#################################################
# without NPDR trials

bcva.1yr.pdr <- filter(bcva.1yr, !Trial %in% c("PROTOCOL W","PANORAMA")) 
bcva.2yr.pdr <- filter(bcva.2yr, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.max.pdr <- filter(bcva.max, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.tb.pdr <- filter(bcva.tb, !Trial %in% c("PROTOCOL W","PANORAMA"))

#################################################
# Analyses at 1 and 2 years

# add numeric Trial/treat
bcva.1yr.pdr$Trialn <- as.numeric(factor(bcva.1yr.pdr$Trial))
bcva.1yr.pdr$Trtn <- as.numeric(fct_infreq(bcva.1yr.pdr$Arm2))
bcva.2yr.pdr$Trialn <- as.numeric(factor(bcva.2yr.pdr$Trial))
bcva.2yr.pdr$Trtn <- as.numeric(fct_infreq(bcva.2yr.pdr$Arm2))

trt.names <- levels(fct_infreq(bcva.1yr.pdr$Arm2))

NMA.1yr.setup <- set_agd_arm(bcva.1yr.pdr,
                         study=Trialn,
                         trt=Trtn,
                         y=logMAR.mcfb,
                         se=logMAR.mcfb.sem,
                         sample_size=N,
                         trt_ref="1")

NMA.1yr <- nma(NMA.1yr.setup,
              trt_effects="fixed",
              prior_intercept=normal(scale=10),
              prior_trt=normal(scale=10))

thresh.res.1yr <- threshold_run(NMA.1yr,bcva.1yr.pdr,trt.names,
                                "Threshold plot 1yr.png")

trt.names <- levels(fct_infreq(bcva.2yr.pdr$Arm2))

NMA.2yr.setup <- set_agd_arm(bcva.2yr.pdr,
                             study=Trialn,
                             trt=Trtn,
                             y=logMAR.mcfb,
                             se=logMAR.mcfb.sem,
                             sample_size=N,
                             trt_ref="1")

NMA.2yr <- nma(NMA.2yr.setup,
               trt_effects="fixed",
               prior_intercept=normal(scale=10),
               prior_trt=normal(scale=10))

thresh.res.2yr <- threshold_run(NMA.2yr,bcva.2yr.pdr,trt.names,
                                "Threshold plot 2yr.png")


##################################################
# Analysis at max FU time

bcva.max.pdr$Trialn <- as.numeric(factor(bcva.max.pdr$Trial))
bcva.max.pdr$Trtn <- as.numeric(fct_infreq(bcva.max.pdr$Arm2))

trt.names <- levels(fct_infreq(bcva.max.pdr$Arm2))

NMA.max.setup <- set_agd_arm(bcva.max.pdr,
                             study=Trialn,
                             trt=Trtn,
                             y=logMAR.mcfb,
                             se=logMAR.mcfb.sem,
                             sample_size=N,
                             trt_ref="1")

NMA.max <- nma(NMA.max.setup,
               trt_effects="fixed",
               prior_intercept=normal(scale=10),
               prior_trt=normal(scale=10))

thresh.res.max <- threshold_run(NMA.max,bcva.max.pdr,trt.names,
                                "Threshold plot max FU.png")


##################################################
# Analysis adjusting for time/baseline

# add numeric Trial/treat
bcva.tb.pdr$Trialn <- as.numeric(factor(bcva.tb.pdr$Trial))
bcva.tb.pdr$Trtn <- as.numeric(fct_infreq(bcva.tb.pdr$Arm2))
bcva.tb.pdr <- filter(bcva.tb.pdr, !is.na(logMAR.base))
trt.names <- levels(fct_infreq(bcva.tb.pdr$Arm2))

NMA.adj.setup <- set_agd_arm(bcva.tb.pdr,
                          study=Trialn,
                          trt=Trtn,
                          y=logMAR.mcfb,
                          se=logMAR.mcfb.sem,
                          sample_size=N,
                          trt_ref="1",
                          trt_class=Class)

NMA.adj <- nma(NMA.adj.setup,
                     trt_effects="fixed",
                     regression= ~(logMAR.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     center=T)

ref.data <- data.frame(logMAR.base=mean(bcva.tb.pdr$logMAR.base,na.rm=T),
                       Year=2)

thresh.res.adj <- threshold_run(NMA.adj,bcva.tb.pdr,trt.names,
                                "Threshold plot time and base adjusted.png",
                                adjusted=T,ref.data=ref.data)


###########################################
# reduced network to compare anti-VEGFs (anti-VEGF = anti-VEGF+PRP)

bcva.tb.pdr$Trtn <- as.numeric(fct_infreq(bcva.tb.pdr$Drug))
trt.names <- levels(fct_infreq(bcva.tb.pdr$Drug))

# unadjusted

NMA.red.setup <- set_agd_arm(bcva.tb.pdr,
                                 study=Trialn,
                                 trt=Trtn,
                                 y=logMAR.mcfb,
                                 se=logMAR.mcfb.sem,
                                 sample_size=N,
                                 trt_ref="1",
                                 trt_class=Class2)

NMA.red <- nma(NMA.red.setup,
                   trt_effects="fixed",
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10))

thresh.res.red <- threshold_run(NMA.red,bcva.tb.pdr,trt.names,
                                "Threshold plot drug model.png")

# adjusted

NMA.red.adj <- nma(NMA.red.setup,
                     trt_effects="fixed",
                     regression= ~(logMAR.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     center=T)

ref.data <- data.frame(logMAR.base=mean(bcva.tb.pdr$logMAR.base,na.rm=T),
                       Year=2)

thresh.res.red.adj <- threshold_run(NMA.red.adj,bcva.tb.pdr,trt.names,
                                "Threshold plot adjusted drug model.png",
                                adjusted=T,ref.data=ref.data)


###########################################
# class network (all anti-VEGFs equivalent)

bcva.tb.pdr$Trtn <- as.numeric(fct_infreq(bcva.tb.pdr$Class))
trt.names <- levels(fct_infreq(bcva.tb.pdr$Class))

# unadjusted

NMA.class.setup <- set_agd_arm(bcva.tb.pdr,
                                   study=Trialn,
                                   trt=Trtn,
                                   y=logMAR.mcfb,
                                   se=logMAR.mcfb.sem,
                                   sample_size=N,
                                   trt_ref="1")

NMA.class <- nma(NMA.class.setup,
                     trt_effects="fixed",
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2))

thresh.res.class <- threshold_run(NMA.class,bcva.tb.pdr,trt.names,
                                      "Threshold plot class model.png")

# adjusted

NMA.class.adj <- nma(NMA.class.setup,
                   trt_effects="fixed",
                   regression= ~(logMAR.base + Year),
                   prior_intercept=normal(scale=10),
                   prior_trt=normal(scale=10),
                   prior_reg=normal(scale=2),
                   center=T)

ref.data <- data.frame(logMAR.base=mean(bcva.tb.pdr$logMAR.base,na.rm=T),
                       Year=2)

thresh.res.class.adj <- threshold_run(NMA.class.adj,bcva.tb.pdr,trt.names,
                                    "Threshold plot adjusted class model.png",
                                    adjusted=T,ref.data=ref.data)
