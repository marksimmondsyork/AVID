# NMA of BCVA data - time and baseline adjustments

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

source("R code/support code/NMA data set-up.R")

#################################################
# without NPDR trials

bcva.max.pdr <- filter(bcva.max, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.tb.pdr <- filter(bcva.tb, !Trial %in% c("PROTOCOL W","PANORAMA")) %>%
  filter(Weeks<150) %>%
  group_by(Trial,Arm) %>%
  slice_max(Weeks)

#################################################
# Treatment class interacting with time

NMA.t1 <- set_agd_arm(bcva.tb.pdr,
                         study=Trial,
                         trt=Arm2,
                         y=etdrs.mcfb,
                         se=etdrs.mcfb.sem,
                         sample_size=N,
                         trt_ref="PRP",
                         trt_class=Class2)

NMA.time <- nma(NMA.t1,
                      trt_effects="random",
                      regression= ~Year:.trtclass,
                      prior_intercept=normal(scale=10),
                      prior_trt=normal(scale=10),
                      prior_reg=normal(scale=2),
                      center=F)

ref.data <- data.frame(etdrs.base=mean(bcva.tb.pdr$etdrs.base,na.rm=T),
                       Year=0)

NMA.time.res <- relative_effects(NMA.time,newdata=ref.data,all_contrasts=TRUE)
NMA.time.ranks <- posterior_rank_probs(NMA.time,newdata=ref.data)

tiff("BCVA ETDRS NMA with time - relative effects forest.tiff",width=1200,height=800,res=200,compression="lzw")
plot(NMA.time.res, ref_line=0, .width=c(0,0.95))
dev.off()

png("BCVA ETDRS NMA with time - prob rankings.png",width=600,height=400)
plot(NMA.time.ranks)
dev.off()

sink("BCVA ETDRS NMA with time - full results.txt")
print(NMA.time)
sink()


###########################
# Treatment class interaction with baseline etdrs and time

NMA.b1 <- set_agd_arm(filter(bcva.tb.pdr,!is.na(etdrs.base)),
                          study=Trial,
                          trt=Arm2,
                          y=etdrs.mcfb,
                          se=etdrs.mcfb.sem,
                          sample_size=N,
                          trt_ref="PRP",
                          trt_class=Class2)

NMA.base <- nma(NMA.b1,
                     trt_effects="random",
                     regression= ~(etdrs.base + Year):.trtclass,
                     prior_intercept=normal(scale=10),
                     prior_trt=normal(scale=10),
                     prior_reg=normal(scale=2),
                     center=T)

NMA.base.res <- relative_effects(NMA.base,newdata=ref.data,all_contrasts=TRUE)
NMA.base.ranks <- posterior_rank_probs(NMA.base,newdata=ref.data)

tiff("BCVA ETDRS NMA time and baseline - relative effects forest.tiff",width=1200,height=800,res=200,compression="lzw")
plot(NMA.base.res, ref_line=0, .width=c(0,0.95))
dev.off()

png("BCVA ETDRS NMA time and baseline - prob rankings.png",width=600,height=400)
plot(NMA.base.ranks)
dev.off()

sink("BCVA ETDRS NMA time and baseline - full results.txt")
print(NMA.base)
sink()

#########################################################
#########################################################
# reduced network

NMA.red <- set_agd_arm(filter(bcva.tb.pdr,!is.na(etdrs.base)),
                          study=Trial,
                          trt=Drug,
                          y=etdrs.mcfb,
                          se=etdrs.mcfb.sem,
                          sample_size=N,
                          trt_ref="PRP",
                          trt_class=Class2)

NMA.bt.red <- nma(NMA.red,
                    trt_effects="random",
                    regression= ~(etdrs.base + Year):.trtclass,
                    prior_intercept=normal(scale=10),
                    prior_trt=normal(scale=10),
                    prior_reg=normal(scale=2),
                    center=T)

NMA.bt.red.res <- relative_effects(NMA.bt.red,newdata=ref.data,all_contrasts=TRUE)
NMA.bt.red.ranks <- posterior_rank_probs(NMA.bt.red,newdata=ref.data)

tiff("BCVA ETDRS adjusted drug model - relative effects forest.tiff",width=1200,height=800,res=200,compression="lzw")
plot(NMA.bt.red.res, ref_line=0, .width=c(0,0.95))
dev.off()

png("BCVA ETDRS adjusted drug model - prob rankings.png",width=600,height=400)
plot(NMA.bt.red.ranks)
dev.off()

sink("BCVA ETDRS adjusted drug model - full results.txt")
print(NMA.bt.red)
sink()

###############################################
# class model 

NMA.red2 <- set_agd_arm(filter(bcva.tb.pdr,!is.na(etdrs.base)),
                       study=Trial,
                       trt=Class,
                       y=etdrs.mcfb,
                       se=etdrs.mcfb.sem,
                       sample_size=N,
                       trt_ref="PRP")

NMA.bt.red2 <- nma(NMA.red2,
                  trt_effects="random",
                  regression= ~(etdrs.base + Year),
                  prior_intercept=normal(scale=10),
                  prior_trt=normal(scale=10),
                  prior_reg=normal(scale=2),
                  center=T)

NMA.bt.red2.res <- relative_effects(NMA.bt.red2,newdata=ref.data,all_contrasts=TRUE)
NMA.bt.red2.ranks <- posterior_rank_probs(NMA.bt.red2,newdata=ref.data)

tiff("BCVA ETDRS adjusted class model - relative effects forest.tiff",width=1200,height=800,res=200,compression="lzw")
plot(NMA.bt.red2.res, ref_line=0, .width=c(0,0.95))
dev.off()

png("BCVA ETDRS adjusted class model - prob rankings.png",width=600,height=400)
plot(NMA.bt.red2.ranks)
dev.off()

sink("BCVA ETDRS adjusted class model - full results.txt")
print(NMA.bt.red2)
sink()

#################################################
# tabulate all results

model.set <- list(NMA.time.res,NMA.base.res,NMA.bt.red.res,NMA.bt.red2.res)
ranks.set <- list(NMA.time.ranks,NMA.base.ranks,NMA.bt.red.ranks,NMA.bt.red2.ranks)
model.names1 <- c("Time by class","Time and baseline by class",
                  "Time and baseline - drug","Time and baseline - class")

m.len <- map_dbl(model.set, ~length(.$summary$parameter))
m.len2 <- map_dbl(ranks.set, ~length(.$summary$parameter))

all.res <- map_dfr(model.set, ~.$summary) %>%
  mutate(Model=rep(model.names1,m.len))
all.ranks <- map_dfr(ranks.set, ~.$summary) %>%
  mutate(Model=rep(model.names1,m.len2))

write_csv(all.res,"BCVA ETDRS NMA adjusted - all effects.csv")
write_csv(all.ranks,"BCVA ETDRS NMA adjusted - all prob ranks.csv")


