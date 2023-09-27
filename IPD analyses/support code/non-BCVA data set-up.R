# Set up binary outcomes data

#############################################
# baseline

ipd.base <- read_csv("Data/AVID IPD baseline.csv") %>% 
  filter(!is.na(Treatment)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye))

ipd.base.2 <- ipd.base %>%
  select(Trial,Arm,PatientID:Treatment,
         Age,Sex,BCVABaseline,Diabetes,PrevAVEGF,PrevPRP,
         VHBaseline,DMEBaseline,HbA1c,HBA1c,CSTBaseline,DRSSBaseline) %>%
  mutate(Diabetes=ifelse(Diabetes==1, 1, 0),
         PrevAVEGF=ifelse(PrevAVEGF>0, 1, 0),
         HbA1c=ifelse(is.na(HbA1c), HBA1c, HbA1c)) %>%
  mutate(HbA1c=ifelse(grepl("PROTOCOL",Trial), 10.93*HbA1c-23.5, HbA1c)) %>%
  #mutate(across(contains("Prev"), ~replace_na(.x,0))) %>%
  mutate(across(c(Age,Sex,BCVABaseline,HbA1c,CSTBaseline,DRSSBaseline), ~.x - mean(.x,na.rm=T), .names="{.col}.c")) %>%
  select(-HBA1c) %>%
  filter(!is.na(Treatment)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) %>%
  relocate(PatEye)


#######################################################
# DME and VH

bin.data <- read_csv("Data/AVID IPD binary.csv") %>%
  filter(TimeObs>0)

dmevh.data <- bin.data %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12)

dmevh.data.ex1yr <- dmevh.data %>%
  filter(Weeks>45 & Weeks<55) %>%
  group_by(PatEye,Trial,Treatment) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) 

###############################################################
 # other outcomes

ipd.other <- read_csv("Data/AVID IPD follow-up.csv") %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  filter(Weeks<60) %>%
  select(Trial,PatEye,Weeks,TimeObs,Treatment,CST_FU,RIP,SAE,TRD_FU)

