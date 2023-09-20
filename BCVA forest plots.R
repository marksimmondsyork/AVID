library(tidyverse)
library(readxl)
library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

# bcva data from publications
bcva <- read_excel("Data/BCVA data edited.xlsx") %>%
  filter(!is.na(Trial)) %>%
  select(-contains("..")) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "Control",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab"
  ))

####################################################
# just required data for CFB analysis
bcva.2 <- bcva %>%
  select(Trial:Weeks,Drug,contains("mcfb")) %>%
  filter(Weeks>0, !is.na(etdrs.mcfb)) %>%
  filter(!Trial %in% c("RECOVERY","Sao Paulo A"))

# Control arm data
bcva.cont <- bcva.2 %>%
  filter(Drug=="Control") %>%
  select(-Drug) %>%
  rename_with(~paste(.,".c",sep=""),contains("mcfb")) %>%
  rename(N.c=N,Arm.c=Arm)

# expanded full data for analysis
bcva.all <- filter(bcva.2, Drug!="Control") %>%
  left_join(bcva.cont) %>%
  arrange(Arm,Weeks,Trial)

bcva.pdr <- filter(bcva.all, !Trial %in% c("PROTOCOL W","PANORAMA")) %>%
  rename(Intervention=Arm) %>%
  mutate(Time=case_when(
    Weeks<15 ~ "3 months or less",
    Weeks<45 ~ "3 months to 1 year",
    Weeks<=60 ~ "1 year",
    Weeks<200 ~ "2 years",
    Weeks>200 ~ "5 years")) %>%
  mutate(Group=paste(Intervention,"-",Time))

bcva.pdr.short <- bcva.pdr %>%
  group_by(Group,Trial) %>%
  slice_tail(n=1) %>%
  arrange(Intervention,Weeks)

ma.MAR.all <- metacont(N, logMAR.mcfb, logMAR.mcfb.sd,
                       N.c, logMAR.mcfb.c, logMAR.mcfb.sd.c,
                       sm="MD", studlab=Trial, subgroup=Group, data=bcva.pdr.short)

tiff("All logMAR forest plot.tiff",height=2000,width=1200,res=200,compression="lzw")
forest(ma.MAR.all,fixed=F,random=F,print.subgroup.name = F,
       leftcols=c("studlab","Weeks"),subgroup.hetstat=F,overall.hetstat=F,
       xlab="MD in logMAR",label.left="Favours Anti-VEGF",label.right="Favours PRP")
dev.off()
 
####################################################
# Up to 1 year plot
bcva.1yr <- bcva.pdr.short %>%
  filter(Weeks<60) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

ma.lmar.1yr <- metacont(N, logMAR.mcfb, logMAR.mcfb.sd,
                       N.c, logMAR.mcfb.c, logMAR.mcfb.sd.c,
                       sm="MD", studlab=Trial, subgroup=Intervention, 
                       data=bcva.1yr)

ma.etdrs.1yr <- metacont(N, etdrs.mcfb, etdrs.mcfb.sd,
                        N.c, etdrs.mcfb.c, etdrs.mcfb.sd.c,
                        sm="MD", studlab=Trial, subgroup=Intervention, 
                        data=bcva.1yr)

tiff("BCVA logMAR 1 year forest plot.tiff",height=2000,width=1200,res=200,compression="lzw")
forest(ma.lmar.1yr,fixed=F,random=T,print.subgroup.name = F,
       leftcols=c("studlab","Weeks"),leftlabs=c("Trial","Weeks"),
       subgroup.hetstat=F,overall.hetstat=F, sort.subgroup=T,test.subgroup=F,
       xlab="MD in logMAR",label.left="Favours Anti-VEGF",label.right="Favours PRP")
dev.off()

tiff("BCVA ETDRS 1 year forest plot.tiff",height=2000,width=1200,res=200,compression="lzw")
forest(ma.etdrs.1yr,fixed=F,random=F,print.subgroup.name = F,
       leftcols=c("studlab","Weeks"),leftlabs=c("Trial","Weeks"),
       subgroup.hetstat=F,overall.hetstat=F, sort.subgroup=T,test.subgroup=F,
       xlab="MD in logMAR",label.right="Favours Anti-VEGF",label.left="Favours PRP")
dev.off()

#########################################
# 1 to 2 years

bcva.2yr <- bcva.pdr.short %>%
  filter(Weeks>45 & Weeks<120) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

ma.lmar.2yr <- metacont(N, logMAR.mcfb, logMAR.mcfb.sd,
                        N.c, logMAR.mcfb.c, logMAR.mcfb.sd.c,
                        sm="MD", studlab=Trial, subgroup=Intervention, 
                        data=bcva.2yr)

ma.etdrs.2yr <- metacont(N, etdrs.mcfb, etdrs.mcfb.sd,
                         N.c, etdrs.mcfb.c, etdrs.mcfb.sd.c,
                         sm="MD", studlab=Trial, subgroup=Intervention, 
                         data=bcva.2yr)

tiff("BCVA logMAR 2 year forest plot.tiff",height=2000,width=1200,res=200,compression="lzw")
forest(ma.lmar.2yr,fixed=F,random=F,print.subgroup.name = F,
       leftcols=c("studlab","Weeks"),leftlabs=c("Trial","Weeks"),
       subgroup.hetstat=F,overall.hetstat=F, sort.subgroup=T,test.subgroup=F,
       xlab="MD in logMAR",label.left="Favours Anti-VEGF",label.right="Favours PRP")
dev.off()

tiff("BCVA ETDRS 2 year forest plot.tiff",height=2000,width=1200,res=200,compression="lzw")
forest(ma.etdrs.2yr,fixed=F,random=F,print.subgroup.name = F,
       leftcols=c("studlab","Weeks"),leftlabs=c("Trial","Weeks"),
       subgroup.hetstat=F,overall.hetstat=F, sort.subgroup=T,test.subgroup=F,
       xlab="MD in logMAR",label.right="Favours Anti-VEGF",label.left="Favours PRP")
dev.off()

#############################################
# all data

ma.etdrs.all <- metacont(N, etdrs.mcfb, etdrs.mcfb.sd,
                         N.c, etdrs.mcfb.c, etdrs.mcfb.sd.c,
                         sm="MD", studlab=Trial, subgroup=Drug, data=bcva.all)

png("All ETDRS forest plot with NPDR.png",width=600, height=700)
forest(ma.etdrs.all,fixed=F,random=F,
       leftcols=c("studlab","Arm","Weeks"))
dev.off()

ma.MAR.all <- metacont(N, logMAR.mcfb, logMAR.mcfb.sd,
                       N.c, logMAR.mcfb.c, logMAR.mcfb.sd.c,
                       sm="MD", studlab=Trial, subgroup=Drug, data=bcva.all)

png("All logMAR forest plot with NPDR.png",width=600, height=700)
forest(ma.MAR.all,fixed=F,random=F,
       leftcols=c("studlab","Arm","Weeks"))
dev.off()