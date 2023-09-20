# Published data MA for BCVA
# various plots and examinations of data

library(tidyverse)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

source("R code/support code/NMA data set-up.R")
source("R code/support code/ggplot slides.R")

# plot all data
png("All ETDRS data by arm.png",width=1000,height=600)
ggplot(bcva, aes(x=Weeks,y=etdrs.mcfb,colour=Drug,shape=Class,size=N,group=Trial)) +
  geom_point() +
  scale_x_continuous("Follow-up (Weeks)") +
  scale_y_continuous("ETDRS (Mean change from baseline)")
dev.off()

png("All ETDRS data by trial.png",width=1000,height=600)
ggplot(bcva, aes(x=Weeks,y=etdrs.mcfb,colour=Trial,shape=Drug,size=N)) +
  geom_point() +
  scale_x_continuous("Follow-up (Weeks)") +
  scale_y_continuous("ETDRS (Mean change from baseline)")
dev.off()

# Control arm data
bcva.cont <- bcva.2 %>%
  filter(Drug=="PRP") %>%
  select(-Drug) %>%
  rename_with(~paste(.,".c",sep=""),contains("mcfb")) %>%
  rename(N.c=N,Arm.c=Arm,Arm2.c=Arm2,Class.c=Class,Class2.c=Class2)

# expanded full data for analysis
bcva.all <- filter(bcva.2, Drug!="PRP", Trial %in% bcva.cont$Trial) %>%
  left_join(bcva.cont) %>%
  arrange(Arm,Weeks,Trial)

# difference data
ma.etdrs.all <- metacont(N, etdrs.mcfb, etdrs.mcfb.sd,
                         N.c, etdrs.mcfb.c, etdrs.mcfb.sd.c,
                         sm="MD", studlab=Trial, subgroup=Drug, data=bcva.all)

bcva.all$MD <- ma.etdrs.all$TE
bcva.all$MD.low <- ma.etdrs.all$lower
bcva.all$MD.high <- ma.etdrs.all$lupper

base.data <- bcva %>%
  filter(Weeks==0) %>%
  select(Trial,Arm,etdrs.mean) %>%
  rename(etdrs.base=etdrs.mean) %>%
  filter(!Trial %in% c("RECOVERY","Sao Paulo A"))

bcva.all <- left_join(bcva.all,base.data) %>%
  mutate(Months=Weeks/52*12)

png("ETDRS MD vs Time.png",width=1000,height=600)
ggplot(bcva.all, aes(x=Months,y=MD,colour=Drug,shape=Drug,size=N)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  scale_x_continuous("Follow-up (Months)") +
  scale_y_continuous("Mean difference in ETDRS beween Anti-VEGF and Control") +
  scale_size_area(max_size=5)
dev.off()

png("ETDRS MD vs baseline ETDRS.png",width=1000,height=600)
ggplot(bcva.all, aes(x=etdrs.base,y=MD,colour=Drug,shape=Drug,size=N)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  scale_x_continuous("Mean ETDRS at baseline") +
  scale_y_continuous("Mean difference in ETDRS beween Anti-VEGF and Control") +
  scale_size_area(max_size=10)
dev.off()

#################################################
# as above but PDR only

bcva.pdr <- filter(bcva, !Trial %in% c("PROTOCOL W","PANORAMA"))
bcva.all.pdr <- filter(bcva.all, !Trial %in% c("PROTOCOL W","PANORAMA"))

# plot all data
png("All ETDRS data by arm - PDR only.png",width=1000,height=600)
ggplot(bcva.pdr, aes(x=Weeks,y=etdrs.mcfb,colour=Drug,shape=Class,size=N,group=Trial)) +
  geom_point() +
  scale_x_continuous("Follow-up (Weeks)") +
  scale_y_continuous("ETDRS (Mean change from baseline)")
dev.off()

png("All ETDRS data by trial - PDR only.png",width=1000,height=600)
ggplot(bcva.pdr, aes(x=Weeks,y=etdrs.mcfb,colour=Trial,shape=Drug,size=N)) +
  geom_point() +
  scale_x_continuous("Follow-up (Weeks)") +
  scale_y_continuous("ETDRS (Mean change from baseline)")
dev.off()


png("ETDRS MD vs Time - PDR only.png",width=1000,height=600)
ggplot(bcva.all.pdr, aes(x=Months,y=MD,colour=Drug,shape=Drug,size=N)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  scale_x_continuous("Follow-up (Months)") +
  scale_y_continuous("Mean difference in ETDRS beween Anti-VEGF and Control") +
  scale_size_area(max_size=5)
dev.off()

png("ETDRS MD vs baseline ETDRS - PDR only.png",width=1000,height=600)
ggplot(bcva.all.pdr, aes(x=etdrs.base,y=MD,colour=Drug,shape=Drug,size=N)) +
  geom_point() +
  scale_x_continuous("Mean ETDRS at baseline") +
  scale_y_continuous("Mean difference in ETDRS beween Anti-VEGF and Control") +
  scale_size_area(max_size=5)
dev.off()
