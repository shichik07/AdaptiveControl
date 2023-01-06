setwd('E:/AdaptiveControl/Data/BehaviorResults')

# Load packages
library(BayesFactor)
library(dplyr)
library(tidyr)
library(rstudioapi)

#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32936)

#Load the data 

filename <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/participants.tsv"
Clinic <- read.csv(file = filename, header = TRUE)
Clinic <- as_tibble(Clinic)

path1 = "Z:/JuliusKricheldorff/AdaptiveControl/BehavioralData"
files = list.files(path = path1, pattern = "AdapCon_PD|AdapCon_CO", full.names = TRUE, recursive = FALSE)
df = tibble()

for (x in files) {
  
  t <- read.csv(x, header = TRUE)
  t <- as_tibble(t)
  name <- substr(x, 62,63)
  num <- regmatches(x, gregexpr("[[:digit:]]+", x))
  subject_id <- paste0('sub-', name, '_', num)
  t$Subject <- subject_id
  if (substr(subject_id, 5,6) == "PD"){
    updrs_ind <- which(Clinic$participant_id == subject_id)
    if (length(updrs_ind) == 0){
      t$UPDRS = NA
    } else{ 
      if (Clinic$UPDRS[updrs_ind] == "n/a"){
        t$UPDRS = NA
      } else {
        t$UPDRS <- as.numeric(Clinic$UPDRS[updrs_ind])
      }
    }
  } else {
    t$UPDRS = NA
  }
  
  df <- bind_rows(df, t)
}



Data <- df %>% 
  mutate(number_1 = as.character(number_1))%>%
  mutate(number_2 = as.character(number_2))%>%
  unite('Item', number_1:number_2, remove = FALSE) %>%
  mutate(Item = as_factor(Item)) %>%
  mutate(Item_specific = ifelse(Analysis_type == "main_con" | Analysis_type == "main_incon", 'Item_spec', 'List_wide' )) %>%
  mutate( Congruency_main = ifelse(Analysis_type == "main_con"| Analysis_type == "MC", "M_Congruent", "M_Incongruent")) %>% 
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main, UPDRS) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response) %>%
  filter(Subject != "sub-PD_c(\"11\", \"2\")" &  Subject != "sub-PD_16") %>%#exclude participant 16 and the partial dataset by participant
  filter(Subject != "sub-CO_0" &  Subject != "sub-CO_3" & Subject != "sub-CO_7" &  Subject != "sub-CO_9") %>% # exclude because they used both hands
  mutate(Error = case_when( # we define accuracy in terms of errors, easier to interpret
    Correct_A == 1 ~ 0,
    Correct_A == 0 ~1
  )) %>%
  filter( RT > 200) %>%
  filter(Block != "practice")# Remove data with small RT, otherwise messes with our non-decision parameter

# Calculate the BF
TbySub <- Data %>%
  group_by(Exp_group, Subject) %>%
  summarise(mean_RT = mean(RT)) %>%
  ungroup()

TbySub <- as.data.frame(TbySub)
TbySub$Exp_group <- unlist(TbySub$Exp_group)
TbySub$mean_RT <- unlist(TbySub$mean_RT)
TbySub$Subject <- unlist(TbySub$Subject)

# Get the Bayesfactor 
BF_group <- ttestBF(formula = mean_RT ~ Exp_group, data = TbySub)

CO_RT <- TbySub %>%
  filter(Exp_group == "CO")
PD_RT <- TbySub %>%
  filter(Exp_group == "PD")
# Get the corresponding p-value
t.test(CO_RT$mean_RT, PD_RT$mean_RT, alternative = "two.sided", var.equal = FALSE)

# Now calculate the congruency effects
TbySub_con <- Data %>%
  group_by(Exp_group, Subject) %>%
  summarize(av.Con = mean(RT[Congruency == "congruent"]),
            av.Incon = mean(RT[Congruency == "incongruent"]), 
            av.RT = (av.Con+ av.Incon)/2,
            Effect = av.Incon - av.Con) %>%
  ungroup()

# Get the Bayesfactor 
BF_group <- ttestBF(formula = Effect ~ Exp_group, data = TbySub_con)

CO_RT <- TbySub_con %>%
  filter(Exp_group == "CO")
PD_RT <- TbySub_con %>%
  filter(Exp_group == "PD")
# Get the corresponding p-value
t.test(CO_RT$Effect, PD_RT$Effect, alternative = "two.sided", var.equal = FALSE)


