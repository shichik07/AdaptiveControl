# Analysis on the behavioral data of healthy control participants and participants with Parkinson's disease
setwd('C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data')

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)

#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32946)

#Load the data 

filename <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/participants.tsv"
Clinic <- read.csv(file = filename, header = TRUE, sep = "")
Clinic <- as_tibble(Clinic) %>%
  mutate(age = as.numeric(age)) %>% 
  mutate(MMST = as.numeric(MMST)) %>%
  mutate(UPDRS = as.numeric(UPDRS)) %>%
  mutate(educationyears = as.numeric(educationyears)) %>%
  filter(!grepl('CY', participant_id)) %>%
  filter(handedness != "both") %>%
  filter(participant_id != "sub-PD_16") %>%
  mutate(Exp_group = ifelse(grepl("CO", participant_id), 'CO', 'PD')) 

Summary_table <- Clinic %>%
  group_by(Exp_group) %>%
  summarize(mean_age = mean(age, na.rm = TRUE), sd_age = sd(age, na.rm = TRUE),  N = n(),
            mean_education = mean(educationyears, na.rm = TRUE), sd_education = sd(educationyears, na.rm = TRUE),
            mean_MMST = mean(MMST, na.rm = TRUE), sd_MMST = sd(MMST, na.rm = TRUE),
            mean_UPDRS = mean(UPDRS, na.rm = TRUE), sd_UPDRS = sd(UPDRS, na.rm = TRUE))
            n_gender = n(gender))

Summary_table <- Clinic %>%
  group_by(Exp_group, gender) %>%
  summarize(n_gender = n())

Summary_ttests <- Clinic %>%
  group_by(Exp_group) %>% 
  mutate(Exp_group = as.factor(Exp_group)) %>%
  summarise(value = list(age)) %>%
  spread(Exp_group, value) %>%
  mutate(p_value = t.test(unlist(CO), unlist(PD))$p.value,
         t_value = t.test(unlist(CO), unlist(PD))$statistic)
