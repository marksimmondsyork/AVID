# Edit IPD for analysis

# all data
ipd.data.base <- read_csv("Data/AVID IPD baseline.csv")
ipd.data.fu <- read_csv("Data/AVID IPD follow-up.csv") %>%
  filter(TimeObs>0)

# edit data
ipd.data.fu1 <- ipd.data.fu %>%
  filter(!is.na(Treatment),!is.na(TimeObs),!is.na(ETDRS_FU)) %>%
  mutate(PatEye=paste(PatientCode,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "PRP",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibizumab")) %>%
  mutate(Class=case_when(
    grepl("+ PRP", Arm) ~ "Anti-VEGF + PRP",
    Arm=="Sham injection" ~ "Sham injection",
    Arm=="PRP" ~ "PRP",
    TRUE ~ "Anti-VEGF"))  %>%
  mutate(Class2=ifelse(Drug=="PRP","PRP","Anti-VEGF"))

# get mean BCVA and change from baseline
ipd.base.bcva <- ipd.data.base %>%
  select(PatientCode,Trial,Eye,Treatment,Arm,ETDRSBaseline,BCVABaseline) %>%
  mutate(PatEye=paste(PatientCode,"-",Eye)) %>%
  rename(etdrs.base=ETDRSBaseline,
         logMAR.base=BCVABaseline) 

ipd.bcva <- ipd.data.fu1 %>%
  select(PatEye,Trial,Arm,Drug,Class,Class2,Weeks,ETDRS_FU,BCVA_FU) %>%
  left_join(ipd.base.bcva) %>%
  mutate(etdrs.cfb=ETDRS_FU-etdrs.base,
         logMAR.cfb=BCVA_FU-logMAR.base) %>%
  filter(!is.na(logMAR.cfb))

# Data, up to 1 yr / 2yrs

ipd.bcva.1yr <- ipd.bcva %>%
  filter(Weeks<53) %>%
  group_by(PatEye) %>%
  slice_max(Weeks)

ipd.bcva.ex1yr <- ipd.bcva %>%
  filter(Weeks>45 & Weeks<60) %>%
  group_by(PatEye) %>%
  slice_max(Weeks)

ipd.bcva.2yr <- ipd.bcva %>%
  filter(Weeks>45 & Weeks<110) %>%
  group_by(PatEye) %>%
  slice_max(Weeks)

# Data at max FU time, up to 2 yrs

ipd.bcva.max <- ipd.bcva %>%
  filter(Weeks<110) %>%
  group_by(PatEye) %>%
  slice_max(Weeks)

# Data with multiple time points and baseline BCVA,  
# restricted to 3mo, 6mo, 12mo, 18mo, 24mo

ipd.bcva.multi <- ipd.bcva %>%
  mutate(Month=case_when
         (Weeks<14 ~ 3,
          Weeks<27 ~ 6,
          Weeks<53 ~ 12,
          Weeks<79 ~ 18,
          Weeks<105 ~ 24)) %>%
  filter(!is.na(Month)) %>%
  group_by(PatEye,Month) %>%
  slice_max(Weeks) %>%
  mutate(Weeks2=52*Month/12) %>%
  mutate(Year=Weeks2/52 - 1) %>%
  left_join(ipd.base.bcva)
