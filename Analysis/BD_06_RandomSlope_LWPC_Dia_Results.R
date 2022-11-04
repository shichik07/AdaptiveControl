# Plot Model Outputs

# Analysis on the behavioral data of healthy control participants and participants with Parkinson's disease
setwd('C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data')

# Load packages
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(bayesplot)
library(tidybayes)
library(stringr)
library(RColorBrewer)
library(plotrix) # to calculate standard error
library(haven)
library(ggpubr)
library(forcats) # so we can simply reorder the variables with fct_inorder

# Set a seed for sake of reproducibility
set.seed(32944)

#load("fit_ShiftedLog_PDOnly_inducer_LW.rda") # Load data Listwise Inducer
load("fit_ShiftedLog_PDOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic
#load("fit_ShiftedLog_COOnly_inducer_LW.rda") # Load data Listwise Inducer
load("fit_ShiftedLog_COOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic


# posterior samples aus dem Model entnehmen
#post_LW_Ind_PD <- posterior_samples(fit_ShiftedLog_PDOnly_inducer_LW)
#post_LW_Dia_PD <- posterior_samples(fit_ShiftedLog_PDOnly_diagnostic_LW)
#post_LW_Ind_CO <- posterior_samples(fit_ShiftedLog_COOnly_inducer_LW)
#post_LW_Dia_CO <- posterior_samples(fit_ShiftedLog_COOnly_diagnostic_LW)


#load the clinical data first
Clinic <- read.table(file = "participants.tsv", header = TRUE)
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
  
  df <- bind_rows(df, t)
}




## Clean up the data a little bit
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
  filter( RT > 200) # Remove data with small RT, otherwise messes with our non-decision parameter


#### Effekte Inducer LW ####
# find unique numbers of particpants
part_nr_PD <- Data %>% 
  filter(Exp_group == "PD") %>%
  distinct(Subject)

part_nr_CO <- Data %>% 
  filter(Exp_group == "CO") %>%
  distinct(Subject)
  




# calculate effect estimates in ms
conditional_effect_calc_shift <- function(model, group, part_nr){
  posterior_cor <- tibble(
    Subject = character(),
    MC_Stroop = numeric(),
    MI_Stroop = numeric(),
    ProactiveC = numeric(),
    Item = character(),
    group = character()
  )
  
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  #transpose the participant vector for the for loop
  part_nr <- t(part_nr)
  
  if (group == "PD"){
    for (part in part_nr) {
      # First get the parameternames
      var_intercept <- paste("m_post$`r_Subject[",part, ",Intercept]`", sep = "")
      var_Congruency <- paste("m_post$`r_Subject[",part, ",Contrast_LWPD_Congruency]`", sep = "")
      var_Listwide <- paste("m_post$`r_Subject[",part, ",Contrast_LWPD_Listwide]`", sep = "")
      var_Interaction <- paste("m_post$`r_Subject[",part, ",Contrast_LWPD_LW_Block]`", sep = "")
      
      # calculate marginal in the diagnostic items  
      LW_MC_C <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) +
                       0.5* (m_post$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) +
                       0.5 * (m_post$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  +
                       0.25* (m_post$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MC_I <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) - 
                       0.5* (m_post$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) + 
                       0.5 * (m_post$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                       0.25* (m_post$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MI_C <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) + 
                       0.5* (m_post$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                       0.5 * (m_post$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                       0.25* (m_post$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MI_I <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) - 
                       0.5* (m_post$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                       0.5 * (m_post$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  + 
                       0.25* (m_post$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      
      # Calculate difference Interaction
      Stroop_MC <- (LW_MC_I - LW_MC_C)
      Stroop_MI <- (LW_MI_I - LW_MI_C)
      LW_Block <- Stroop_MC - Stroop_MI
      
      # save diagnostic variables
      posterior_cor <- posterior_cor %>% add_row(Subject = part,
                                                 MC_Stroop = mean(Stroop_MC),
                                                 MI_Stroop = mean(Stroop_MI),
                                                 ProactiveC = mean(LW_Block),
                                                 Item = "Diagnostic",
                                                 group = "PD")
      remove(Stroop_MC)
      
    }
    return(posterior_cor)
  } else if(group == "CO"){
    for (part in part_nr) {
      # First get the parameternames
      var_intercept <- paste("m_post$`r_Subject[",part, ",Intercept]`", sep = "")
      var_Congruency <- paste("m_post$`r_Subject[",part, ",Contrast_LWCO_Congruency]`", sep = "")
      var_Listwide <- paste("m_post$`r_Subject[",part, ",Contrast_LWCO_Listwide]`", sep = "")
      var_Interaction <- paste("m_post$`r_Subject[",part, ",Contrast_LWCO_LW_Block]`", sep = "")
      
      # calculate marginal in the diagnostic items  
      LW_MC_C <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) +
                       0.5* (m_post$b_Contrast_LWCO_Congruency + eval(parse(text = var_Congruency))) +
                       0.5 * (m_post$b_Contrast_LWCO_Listwide + eval(parse(text = var_Listwide)))  +
                       0.25* (m_post$b_Contrast_LWCO_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MC_I <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) - 
                       0.5* (m_post$b_Contrast_LWCO_Congruency + eval(parse(text = var_Congruency))) + 
                       0.5 * (m_post$b_Contrast_LWCO_Listwide + eval(parse(text = var_Listwide)))  - 
                       0.25* (m_post$b_Contrast_LWCO_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MI_C <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) + 
                       0.5* (m_post$b_Contrast_LWCO_Congruency + eval(parse(text = var_Congruency))) - 
                       0.5 * (m_post$b_Contrast_LWCO_Listwide + eval(parse(text = var_Listwide)))  - 
                       0.25* (m_post$b_Contrast_LWCO_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      LW_MI_I <- exp((m_post$b_Intercept + eval(parse(text = var_intercept))) - 
                       0.5* (m_post$b_Contrast_LWCO_Congruency + eval(parse(text = var_Congruency))) - 
                       0.5 * (m_post$b_Contrast_LWCO_Listwide + eval(parse(text = var_Listwide)))  + 
                       0.25* (m_post$b_Contrast_LWCO_LW_Block + eval(parse(text = var_Interaction))) + 
                       m_post$sigma/2)
      
      # Calculate difference Interaction
      Stroop_MC <- (LW_MC_I - LW_MC_C)
      Stroop_MI <- (LW_MI_I - LW_MI_C)
      LW_Block <- Stroop_MC - Stroop_MI
      
      # save diagnostic variables
      posterior_cor <- posterior_cor %>% add_row(Subject = part,
                                                 MC_Stroop = mean(Stroop_MC),
                                                 MI_Stroop = mean(Stroop_MI),
                                                 ProactiveC = mean(LW_Block),
                                                 Item = "Diagnostic",
                                                 group = "CO")
      remove(Stroop_MC)
      
    }
    return(posterior_cor)
    
  }
}
      
PD_effects = conditional_effect_calc_shift(model = fit_ShiftedLog_PDOnly_diagnostic_LW, 
                              group = "PD", 
                              part_nr = part_nr_PD)
CO_effects = conditional_effect_calc_shift(model = fit_ShiftedLog_COOnly_diagnostic_LW, 
                                           group = "CO", 
                                           part_nr = part_nr_CO)
