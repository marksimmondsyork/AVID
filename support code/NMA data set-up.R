# general data set-up and editing for MA and NMA

library(tidyverse)
library(magrittr)
library(readxl)
library(meta)
library(lme4)
library(multinma)
library(nmathresh)

options(mc.cores = parallel::detectCores())

setwd("S:/CRD/NIHR HTA Anti-VEGF Drugs Diabetic Retinopathy/Meta-analysis/Published data")

# BCVA data from publications
bcva <- read_excel("data/BCVA data edited.xlsx") %>%
  filter(!is.na(Trial)) %>%
  select(-contains("..")) %>%
  mutate(Drug=case_when(
    Arm %in% c("PRP","Sham injection") ~ "PRP",
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("Beva",Arm) ~ "Bevacizumab",
    grepl("Rani",Arm) ~ "Ranibuzimab")) %>%
  mutate(Arm2=case_when(
    grepl("Aflib",Arm) ~ "Aflibercept",
    grepl("ETDRS",Arm) ~ "Ranibizumab + PRP",
    grepl("PASCAL",Arm) ~ "Ranibizumab + PRP",
    TRUE ~ Arm)) %>%
  mutate(Class=case_when(
    grepl("+ PRP", Arm) ~ "Anti-VEGF + PRP",
    Arm=="Sham injection" ~ "Sham injection",
    Arm=="PRP" ~ "PRP",
    TRUE ~ "Anti-VEGF"))  %>%
  mutate(Class2=ifelse(Drug=="PRP","PRP","Anti-VEGF"))

bcva.2 <-  bcva %>%
  select(Trial:Weeks,Arm2,Drug,Class,Class2,contains("mcfb")) %>%
  filter(Weeks>0, !is.na(etdrs.mcfb),Trial!="RECOVERY") %>%
  arrange(Trial,Drug)

# Data, up to 1 yr / 2yrs

bcva.1yr <- bcva.2 %>%
  filter(Weeks<53) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

bcva.2yr <- bcva.2 %>%
  filter(Weeks>45 & Weeks<120) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

# Data at max FU time

bcva.max <-  bcva.2 %>%
  filter(Weeks<120) %>%
  group_by(Trial) %>%
  slice_max(Weeks)

# data adjusted for time and baseline logMAR

base.data <- bcva %>%
  filter(Weeks==0) %>%
  select(Trial,Arm,logMAR,etdrs.mean) %>%
  rename(logMAR.base=logMAR, etdrs.base=etdrs.mean) %>%
  filter(!Trial %in% c("RECOVERY","Sao Paulo A"))

bcva.tb <- bcva.2 %>%
  mutate(Year=Weeks/52 - 1) %>%
  left_join(base.data)
