#####
# Author: Julius Kricheldorff
# Exploring time series effects of adaptive control in the LWPC effect for PD and CO seperately

####
setwd('E:/AdaptiveControl/Data/BehaviorResults')

# Load packages
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstudioapi)
library(xtable)
library(ggplot2)

# Set a seed for sake of reproducibility
set.seed(32946)

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
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main, UPDRS, total_trial) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response) %>%
  filter(Subject != "sub-PD_c(\"11\", \"2\")" &  Subject != "sub-PD_16") %>%#exclude participant 16 and the partial dataset by participant
  filter(Subject != "sub-CO_0" &  Subject != "sub-CO_3" & Subject != "sub-CO_7" &  Subject != "sub-CO_9") %>% # exclude because they used both hands
  mutate(Error = case_when( # we define accuracy in terms of errors, easier to interpret
    Correct_A == 1 ~ 0,
    Correct_A == 0 ~1
  )) 

Data_one_P <- Data %>%
  filter(Subject == "sub-CO_1",
         Block != "practice",
         Item_specific == "List_wide") %>%
  group_by(Analysis_type, Block, Trl_type, Congruency) %>%
  summarize(N = n())
  

Data_correct <- Data %>%
  filter( RT > 200) %>% # Remove data with small RT, otherwise messes with our non-decision parameter 
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group)) %>% # rank items by group
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency)%>%
  mutate(rank = dense_rank(total_trial)) %>% # get rank of individual trials
  mutate(quant_rank = ntile(rank,4)) %>% # get quantile rank of items
  ungroup() %>%
  mutate(Analysis_type = case_when(
    Analysis_type == "main_con" | Analysis_type == "MC" ~ "MC",
    Analysis_type == "main_incon" | Analysis_type == "MI" ~ "MI"
  ))


#### Next we want to calculate the average reaction time by percentile of items:
# Here we go into 

SUM_RT  <- Data_correct %>% 
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  mutate(Av_RT = mean(RT)) %>%
  ungroup() %>% 
  group_by(Exp_group, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  mutate(grand_RT = mean(Av_RT)) %>% # calculate grand average 
  ungroup() %>%
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  mutate(corrected_RT = RT - Av_RT + grand_RT) %>% # correct RT for within subject variance
  ungroup() %>%
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>% # now calculate the average again with the corrected RT
  ungroup() %>%
  mutate(Subject  = as_factor(Subject)) %>%
  group_by(Exp_group, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  summarize(SE = sd(corrected_RT)/n_distinct(Subject), Group_RT=  mean(RT)) # average RT by group including within subject corrected SEs
  
## Next let us plot the results for the LWPC effect first
save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"
LWPC_SUM <- SUM_RT %>%
  filter(Item_specific == "List_wide")

ggplot(data = LWPC_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
  facet_grid(Exp_group ~ Analysis_type) +
  geom_line(size = 1) + 
  geom_point() +
  coord_cartesian(ylim = c(600, 820)) +
  labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item Type") +
  ggtitle("RT over Time LWPC (proactive C.)") + 
  theme(
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  scale_linetype_manual(values = c(4,1))
  geom_errorbar(aes(ymin = - SE, ymax = SE), position = position_dodge(width = 0.9))
  save_n <- "RTbyTime_LWPC.png"
  
  ggsave(path = save_path, filename = save_n,  dpi=300)
  
  
  ISPC_SUM <- SUM_RT %>%
    filter(Item_specific != "List_wide")
  
  ggplot(data = ISPC_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
    facet_grid(Exp_group ~ Analysis_type) +
    geom_line(size = 1) + 
    geom_point() +
    coord_cartesian(ylim = c(600, 820)) +
    labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item Type") +
    ggtitle("RT over Time ISPC (reactive C.)") + 
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title.position  =  "panel"
    ) +
    scale_color_manual(values = c("#74C476","#00441B")) +
    scale_linetype_manual(values = c(4,1))
  geom_errorbar(aes(ymin = - SE, ymax = SE), position = position_dodge(width = 0.9))

  save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"
  save_n <- "RTbyTime_ISPC.png"
  
  ggsave(path = save_path, filename = save_n,  dpi=300)
  