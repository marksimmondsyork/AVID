# Set up binary outcomes data

###########################################################
# Aggregate data set-up
ag.data <- read_excel("data/AD Other Outcome data.xlsx") %>%
  filter(!is.na(Trial), Weeks>0) %>%
  select(-N,-contains("NV")) %>%
  rename(Arm.orig=Arm) %>%
  mutate(Arm=case_when(
    grepl("Aflib",Arm.orig) ~ "Aflibercept",
    grepl("ETDRS",Arm.orig) ~ "Ranibizumab + PRP",
    grepl("PASCAL",Arm.orig) ~ "Ranibizumab + PRP",
    TRUE ~ Arm.orig)) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "PRP",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab")) %>%
  mutate(Class=case_when(
    grepl("+ PRP", Arm) ~ "Anti-VEGF + PRP",
    Arm=="Sham injection" ~ "Sham injection",
    Arm=="PRP" ~ "PRP",
    TRUE ~ "Anti-VEGF"))  %>%
  mutate(Class2=ifelse(Drug=="PRP","PRP","Anti-VEGF"))

# trials to drop from AD set
ipd.trials <- c("CLARITY","PROTOCOL S","PROTEUS")
npdr.trials <- c("PROTOCOL W","PANORAMA")
excl.trials <- c(npdr.trials,"RECOVERY") 

ag.data <-  ag.data %>%
  filter(Weeks>0, !Trial %in% excl.trials) %>%
  arrange(Trial,Drug)

ag.data.1yr <- ag.data %>%
  filter(Weeks<60) %>%
  group_by(Trial) %>%
  slice_max(Weeks) %>%
  ungroup()

ag.data.ex1yr <- ag.data %>%
  filter(Weeks>45 & Weeks<60) %>%
  group_by(Trial) %>%
  slice_max(Weeks) %>%
  ungroup()

#######################################################
# IPD set-up
ipd.data <- read_csv("Data/AVID IPD binary.csv") %>%
  filter(TimeObs>0)

ipd.data.fu <- ipd.data %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  mutate(Arm=case_when(
    Treatment==0 ~ "PRP",
    Treatment==1 & Trial=="CLARITY" ~ "Aflibercept",
    Treatment==1 & Trial=="PROTOCOL S" ~ "Ranibizumab",
    Treatment==1 & Trial=="PROTEUS" ~ "Ranibizumab + PRP",
  )) %>%
  mutate(Drug=case_when(
    Arm=="PRP" ~ "PRP",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Rani",Arm) ~ "Ranibuzimab")) %>%
  mutate(Class=case_when(
    grepl("+ PRP", Arm) ~ "Anti-VEGF + PRP",
    Arm=="PRP" ~ "PRP",
    TRUE ~ "Anti-VEGF"))  %>%
  mutate(Class2=ifelse(Drug=="PRP","PRP","Anti-VEGF")) 

ipd.base <- read_csv("Data/AVID IPD baseline.csv") %>% 
  filter(!is.na(Treatment)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye)) %>%
  mutate(Arm=case_when(
    Treatment==0 ~ "PRP",
    Treatment==1 & Trial=="CLARITY" ~ "Aflibercept",
    Treatment==1 & Trial=="PROTOCOL S" ~ "Ranibizumab",
    Treatment==1 & Trial=="PROTEUS" ~ "Ranibizumab + PRP",
  )) %>%
  select(PatEye,Trial,Arm,VHBaseline,DMEBaseline)

ipd.data.1yr <- ipd.data.fu %>%
  filter(Weeks<53) %>%
  group_by(PatEye,Trial,Arm) %>%
  summarise(across(DME_FU:Vitrectomy, ~ifelse(sum(.x)>0, 1, 0))) %>%
  mutate(Class2=ifelse(Arm=="PRP","PRP","Anti-VEGF"),
         Weeks=52)

ipd.data.ex1yr <- ipd.data.fu %>%
  filter(Weeks>45 & Weeks<55) %>%
  group_by(PatEye,Trial,Arm) %>%
  slice_max(Weeks) %>%
  slice_head(n=1) %>%
  mutate(Class2=ifelse(Arm=="PRP","PRP","Anti-VEGF"),
         Weeks=52)

ipd.data.multi <- ipd.data.fu %>%
  mutate(Month=case_when
         (Weeks<14 ~ 3,
           Weeks<27 ~ 6,
           Weeks<53 ~ 12,
           Weeks<79 ~ 18,
           Weeks<105 ~ 24)) %>%
  filter(!is.na(Month)) %>%
  group_by(PatEye,Month) %>%
  slice_max(Weeks) %>%
  mutate(Year=Weeks/52 - 1)

###############################################################
 # other outcomes

ipd.other <- read_csv("Data/AVID IPD follow-up.csv") %>%
  filter(!is.na(Treatment),!is.na(TimeObs)) %>%
  mutate(PatEye=paste(Trial,"-",PatientID,"-",Eye),
         Weeks=52*TimeObs/12) %>%
  filter(Weeks<60) %>%
  mutate(Arm=case_when(
    Treatment==0 ~ "PRP",
    Treatment==1 & Trial=="CLARITY" ~ "Aflibercept",
    Treatment==1 & Trial=="PROTOCOL S" ~ "Ranibizumab",
    Treatment==1 & Trial=="PROTEUS" ~ "Ranibizumab + PRP",
  )) %>%
  mutate(Drug=case_when(
    Arm=="PRP" ~ "PRP",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Rani",Arm) ~ "Ranibuzimab")) %>%
  mutate(Class=case_when(
    grepl("+ PRP", Arm) ~ "Anti-VEGF + PRP",
    Arm=="PRP" ~ "PRP",
    TRUE ~ "Anti-VEGF"))  %>%
  mutate(Class2=ifelse(Drug=="PRP","PRP","Anti-VEGF")) %>%
  select(PatEye:Class2,)
