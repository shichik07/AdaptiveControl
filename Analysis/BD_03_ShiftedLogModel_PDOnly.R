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
library(purrr)

#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32946)

#### Load the  behavioral data ####
# not we use not participant PD_16 because he slept in during the experiment 

# load the clinical data
Clinic <- read.table(file = "participants.tsv", header = TRUE)
Clinic <- as_tibble(Clinic) %>%
  matches()
  mutate(find_me = any(str_det))
 
detect_index(Clinic, "sub-PD_1")
which(Clinic$participant_id == "sub-PD_1")

path1 = "Z:/JuliusKricheldorff/AdaptiveControl/BehavioralData"
files = list.files(path = path1, pattern = "AdapCon_PD", full.names = TRUE, recursive = FALSE)
df = tibble()

for (x in files) {
  #fil <- x
  t <- read.csv(x, header = TRUE)
  t <- as_tibble(t)
  name <- substr(x, 62,63)
  num <- regmatches(x, gregexpr("[[:digit:]]+", x))
  subject_id <- paste0('sub-', name, '_', num)
  t$Subject <- subject_id
  updrs_ind <- which(Clinic$participant_id == subject_id)
  t$UPDRS <- as.numeric(Clinic$UPDRS[70])
  df <- bind_rows(df, t)
}

## Clean up the data a little bit
df <- df %>% 
  mutate(number_1 = as.character(number_1))%>%
  mutate(number_2 = as.character(number_2))%>%
  unite('Item', number_1:number_2, remove = FALSE) %>%
  mutate(Item = as_factor(Item)) %>%
  mutate(Item_specific = ifelse(Analysis_type == "main_con" | Analysis_type == "main_incon", 'Item_spec', 'List_wide' )) %>%
  mutate( Congruency_main = ifelse(Analysis_type == "main_con"| Analysis_type == "MC", "M_Congruent", "M_Incongruent")) %>% 
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main, UPDRS) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response)

# Contrasts only for the list-wide effect only
StroopCon_LW <- hypr(
  PD_Congruency = (PD_LW_incon_C + PD_LW_con_C)/2 ~ (PD_LW_con_I + PD_LW_incon_I)/2, # main effect Congruence healthy controls
  PD_Listwide = (PD_LW_con_C + PD_LW_con_I)/2 ~ (PD_LW_incon_C + PD_LW_incon_I)/2, # main effect Block Listwide control
  PD_LW_Block = (PD_LW_incon_I - PD_LW_incon_C)/2 ~ (PD_LW_con_I - PD_LW_con_C)/2, # interaction listwide effect
  levels = c("PD_LW_incon_I", "PD_LW_incon_C", "PD_LW_con_C", "PD_LW_con_I")
)
StroopCon_LW

# Contrasts for the Item-specific effects only - check that later (might be more difficult)
StroopCon_IS <- hypr(
  PD_Congruency = (PD_IS_incon_C + PD_IS_con_C)/2 ~ (PD_IS_con_I + PD_IS_incon_I)/2, # main effect Congruence healthy controls
  PD_Itemspecific = (PD_IS_con_C + PD_IS_con_I)/2 ~ (PD_IS_incon_C + PD_IS_incon_I)/2, # main effect Block Listwide control
  PD_IS_Block = (PD_IS_incon_I - PD_IS_incon_C)/2 ~ (PD_IS_con_I - PD_IS_con_C)/2, # interaction listwide effect
  levels = c("PD_IS_con_C", "PD_IS_con_I", "PD_IS_incon_I", "PD_IS_incon_C")
)
StroopCon_IS


# Now we create an extra column for the contrasts for the List wide effect
df <- df %>%
  mutate(Contrast_LW = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "PD"  ~ "PD_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_LW_incon_I",
  )) %>%
  mutate(Contrast_LW  = as_factor(Contrast_LW)) 

# assign the generated contrast matrix to the List Wide Factor
contrasts(df$Contrast_LW) <- contr.hypothesis(StroopCon_LW)

contrasts(df$Contrast_LW)    

# Now we create an extra column for the contrasts for the Item Specific effect
df <- df %>%
  mutate(Contrast_IS = case_when(
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_incon_I",
  )) %>%
  mutate(Contrast_IS  = as_factor(Contrast_IS))

# assign the generated contrast matrix to the Item Specific Factor
contrasts(df$Contrast_IS) <- contr.hypothesis(StroopCon_IS)

contrasts(df$Contrast_IS)    




## Get seperate datasets
Data_y_inducer_LWPC <- df %>%
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic") %>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

Data_y_diagnostic_LWPC <- df %>%
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic") %>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))


Data_y_inducer_ISPC <- df %>%
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

Data_y_diagnostic_ISPC <- df %>%
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

# brmsformula object List Wide
m2_PD_LW <- bf(RT ~ 1 + Contrast_LW + (Contrast_LW|Subject) + (1|Item)) 
a <- get_prior(data = Data_y_inducer_LWPC, m2_PD_LW, family = shifted_lognormal()) 

## Weakly informative prior for the listwide and item specific effect

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(student_t(3, 6.5, 0.5), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0,0.5), class = sigma), 
  prior(cauchy(0, 5), class = ndt),
  prior(lkj(1), class = cor), # use default prior for the correlation 
  prior(lkj(1), class = cor, group = Subject),
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0, 0.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 0.5), class = sd, coef = Contrast_LWPD_Congruency, group = Subject),
  prior(normal(0, 0.5), class = sd, coef = Contrast_LWPD_Listwide, group = Subject),
  prior(normal(0, 0.5), class = sd, coef = Contrast_LWPD_LW_Block, group = Subject),
  prior(normal(0, 0.5), class = sd, coef = Intercept, group = Item)
)

fit_ShiftedLog_PDOnly_inducer_LW <- brm(formula = m2_PD_LW ,
                            family = shifted_lognormal(),
                            data = Data_y_inducer_LWPC,
                            prior = prior_weakly_informed_LW,
                            warmup = 2000,
                            iter = 12000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 423,
                            control = list(adapt_delta = 0.95),
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_ShiftedLog_PDOnly_inducer_LW, file = "fit_ShiftedLog_PDOnly_inducer_LW.rda")


