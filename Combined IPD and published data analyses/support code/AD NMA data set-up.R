# general data set-up and editing for MA and NMA


# BCVA data from publications
ad.bcva <- read_excel("data/AD BCVA data.xlsx") %>%
  filter(!is.na(Trial)) %>%
  select(-contains("..")) %>%
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
excl.trials <- c(ipd.trials,npdr.trials,"RECOVERY")

ad.bcva.2 <-  ad.bcva %>%
  select(Trial:Weeks,Arm,Drug,Class,Class2,contains("mcfb")) %>%
  filter(Weeks>0, !is.na(etdrs.mcfb), !Trial %in% excl.trials) %>%
  arrange(Trial,Drug)

# Data, up to 1 yr / 2yrs

ad.bcva.1yr <- ad.bcva.2 %>%
  filter(Weeks<53) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

ad.bcva.ex1yr <- ad.bcva.2 %>%
  filter(Weeks>45 & Weeks<60) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

ad.bcva.2yr <- ad.bcva.2 %>%
  filter(Weeks>45 & Weeks<120) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

# Data at max FU time

ad.bcva.max <-  ad.bcva.2 %>%
  filter(Weeks<120) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

# data adjusted for time and baseline logMAR

ad.base.data <- ad.bcva %>%
  filter(Weeks==0) %>%
  select(Trial,Arm,logMAR,etdrs.mean) %>%
  rename(logMAR.base=logMAR,etdrs.base=etdrs.mean) %>%
  filter(!Trial %in% excl.trials)

ad.bcva.tb <- ad.bcva.2 %>%
  mutate(Year=Weeks/52 - 1) %>%
  left_join(ad.base.data)

rm(excl.trials,ipd.trials,npdr.trials)
