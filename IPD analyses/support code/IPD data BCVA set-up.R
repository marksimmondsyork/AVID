# Edit AVID IPD data for analysis


# all data

ipd.fu <- read_csv("Data/AVID IPD follow-up.csv")
ipd.base <- read_csv("Data/AVID IPD baseline.csv")

# edit data
ipd.fu1 <- ipd.fu %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) 

ipd.base1 <- ipd.base %>% 
  filter(!is.na(Treatment)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) 

# edit baseline variables

ipd.base.2 <- ipd.base1 %>%
  select(PatEye,Trial,Treatment,Arm,Age,Sex,Diabetes,
         PrevAVEGF,PrevPRP,VHBaseline,DMEBaseline,HbA1c,HBA1c,
         CSTBaseline,DRSSBaseline,NVDBaseline) %>%
  mutate(Diabetes=ifelse(Diabetes==1, 1, 0),
         PrevAVEGF=ifelse(PrevAVEGF>0, 1, 0),
         HbA1c=ifelse(is.na(HbA1c), HBA1c, HbA1c)) %>%
  mutate(HbA1c=ifelse(grepl("PROTOCOL",PatEye), 10.93*HbA1c-23.5, HbA1c)) %>%
  #mutate(across(contains("Prev"), ~replace_na(.x,0))) %>%
  mutate(across(c(Age,Sex,HbA1c,CSTBaseline,DRSSBaseline,NVDBaseline), ~.x - mean(.x,na.rm=T), .names="{.col}.c")) %>%
  select(-HBA1c)
  
# get mean bcva and change from baseline
ipd.base.bcva <- ipd.base1 %>%
  select(PatEye,Trial,Treatment,Arm,ETDRSBaseline,BCVABaseline) %>%
  rename(etdrs.base=ETDRSBaseline,
         logMAR.base=BCVABaseline) %>%
  mutate(across(c(etdrs.base,logMAR.base), ~.x - mean(.x,na.rm=T), .names="{.col}.c"))
  

ipd.bcva <- ipd.fu1 %>%
  select(PatEye,Trial,Treatment,Arm,TimeObs,ETDRS_FU,BCVA_FU) %>%
  filter(TimeObs>0) %>%
  left_join(ipd.base.bcva) %>%
  mutate(etdrs.cfb=ETDRS_FU-etdrs.base,
         logMAR.cfb=BCVA_FU-logMAR.base) %>%
  rename(etdrs=ETDRS_FU,
         logMAR=BCVA_FU) %>%
  filter(!is.na(logMAR.cfb))


# Data with multiple time points and baseline bcva,  
# restricted to annual data

ipd.bcva.yrs <- ipd.bcva %>%
  mutate(Year=ceiling(TimeObs/12)) %>%
  filter(!is.na(Year)) %>%
  group_by(PatEye,Year) %>%
  slice_max(TimeObs) 

rm(ipd.fu1,ipd.base1)
