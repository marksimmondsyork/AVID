# Analyses of NPDR trials

library(meta)

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

source("R code/support code/NMA data set-up.R")

########################################
# NPDR only

bcva.max.npdr <- filter(bcva.max, Trial %in% c("PROTOCOL W","PANORAMA")) %>%
  filter(Arm!="Aflibercept 2q8") %>%
  mutate(across(c(Drug,Class2), ~ifelse(.x=="PRP", "Control",.x)))
bcva.tb.npdr <- filter(bcva.tb, Trial %in% c("PROTOCOL W","PANORAMA")) %>%
  filter(Arm!="Aflibercept 2q8") %>%
  mutate(across(c(Drug,Class2), ~ifelse(.x=="PRP", "Control",.x)))

##########################################
# BCVA

bcva.data <- bcva.max.npdr %>%
  select(Trial,N,Class2,etdrs.mcfb,etdrs.mcfb.sd,logMAR.mcfb,logMAR.mcfb.sd) %>%
  pivot_wider(id_cols=Trial,names_from=Class2,values_from=c(N,contains("etdrs"),contains("MAR")))

npdr.etdrs.ma <- metacont(`N_Anti-VEGF`,`etdrs.mcfb_Anti-VEGF`,`etdrs.mcfb.sd_Anti-VEGF`,
                          N_Control,etdrs.mcfb_Control,etdrs.mcfb.sd_Control,
                          data=bcva.data,sm="MD",studlab=Trial)


npdr.MAR.ma <- metacont(`N_Anti-VEGF`,`logMAR.mcfb_Anti-VEGF`,`logMAR.mcfb.sd_Anti-VEGF`,
                          N_Control,logMAR.mcfb_Control,logMAR.mcfb.sd_Control,
                          data=bcva.data,sm="MD",studlab=Trial)

png("NPDR BCVA LogMAR forest plot.png",width=1000, height=700)
#forest(npdr.MAR.ma,leftcols="studlab",rightcols=c("effect","ci"))
forest(npdr.MAR.ma,rightcols=c("effect","ci"), 
       label.right="Favours sham inj.",label.left="Favours Anti-VEGF")
dev.off()

png("NPDR BCVA ETDRS forest plot.png",width=1000, height=700)
# forest(npdr.etdrs.ma,leftcols="studlab",rightcols=c("effect","ci"))
forest(npdr.etdrs.ma,rightcols=c("effect","ci"),
       label.left="Favours sham inj.",label.right="Favours Anti-VEGF")
dev.off()

############################################
# other outcomes

other.oc <- read_excel("Data/other outcome data.xlsx") %>%
  filter(Trial %in% c("PROTOCOL W","PANORAMA"),
         Arm!="Aflibercept 2q8",
         Weeks>60) %>%
  select_if(~!(all(is.na(.))))

var.names <- names(other.oc)[6:11]

onames.full <- c("DME","Improved DRSS","Other treatment","Proliferative DR",
                 "Vitrectomy","Vitreous Haemorrhage")

# functions to calculate RRs

RR_func <- function(data,var){
   rr.res <- metabin(data[[var]][1],data$N[1],data[[var]][2],data$N[2],sm="RR")
   out <- data.frame(RR=rr.res$TE,RRlow=rr.res$lower,RRhigh=rr.res$upper)
   out <- data.frame(exp(out),logRR=rr.res$TE,selogRR=rr.res$seTE)
   return(out)
 } 

RR_func2 <- function(data,var.names){
  rr.res <- map_dfr(var.names,~RR_func(data,.x)) %>%
    mutate(Outcome=var.names)
  return(rr.res)
}
 
# calculate RRs for all outcomes
RR.res <- other.oc %>%
  group_by(Trial) %>%
  nest() %>%
  mutate(rr.res=map(data,~RR_func2(.x,var.names))) %>%
  select(-data) %>%
  unnest(cols=everything()) %>%
  arrange(Outcome,Trial)

RR.res$Outcome2 <- rep(onames.full,each=2)

png("NPDR other outcomes forest plot.png",width=1000, height=700)
forest(metagen(logRR,selogRR,sm="RR",data=RR.res),
       studlab=Outcome2,fixed=F,random=F,
       leftcols=c("studlab","Trial"),leftlabs=c("Outcome","Trial"),
       xlab="Relative risk",smlab="Relative risk")
dev.off()

#################################
# DME

DME.data <- other.oc %>%
  select(Trial,Arm,N,N.DME) %>%
  mutate(Arm=ifelse(grepl("Afl",Arm), "AV","C")) %>%
  pivot_wider(id_cols=Trial,names_from=Arm,values_from=contains("N"))

DME.ma <- DME.data %$% metabin(N.DME_AV,N_AV,N.DME_C,N_C,sm="RR",studlab=Trial)

png("NPDR DME forest plot.png",width=1000, height=700)
forest(DME.ma,xlab="Relative risk",smlab="Relative risk",
       rightcols=c("effect","ci"),random=F,
       label.left="Favours Anti-VEGF.",label.right="Favours sham inj.")
dev.off()

##########################################
# adverse events

ae.data <- read_excel("Data/AE data.xlsx") %>%
  filter(Trial %in% c("PROTOCOL W","PANORAMA"),
         Arm!="Aflibercept 2q8",
         Weeks>60) %>%
  select_if(~!(all(is.na(.))))

vars.aes <- names(ae.data)[5:13]

aenames.full <- c("Cataracts","Conj. haemorrhage","Death","Myocardial Infarct.",
                  "Ocular pain","Retinal detachment","Stroke","Vitreous haemorrhage",
                  "Vitreous floaters")

# calculate RRs for all outcomes
AE.RR.res <- ae.data %>%
  group_by(Trial) %>%
  nest() %>%
  mutate(rr.res=map(data,~RR_func2(.x,vars.aes))) %>%
  select(-data) %>%
  unnest(cols=everything()) %>%
  arrange(Outcome,Trial)

AE.RR.res$Outcome2 <- rep(aenames.full,each=2)

png("NPDR AE forest plot.png",width=1000, height=700)
forest(metagen(logRR,selogRR,sm="RR",data=AE.RR.res),
       studlab=Outcome2,fixed=F,random=F,
       leftcols=c("studlab","Trial"),leftlabs=c("Outcome","Trial"),
       xlab="Relative risk",smlab="Relative risk",
       label.left="Favours Anti-VEGF.",label.right="Favours sham inj.")
 dev.off()
