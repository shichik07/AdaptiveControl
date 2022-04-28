# Mock Analysis on the behavioral data of healthy control participants

# Load packages
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(bayesplot)
library(tictoc)
library(extraDistr) # more distributions
library(haven)
library(purrr)
library(hypr)

# Set a seed for sake of reproducibility
set.seed(33946)

#Load the data 
filename <- "/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Data/Behavior/All_subjects.csv"
AdapCon <- read.csv(file = filename)

# Get Items, and name blocks correctly
Data <- AdapCon %>%
  as_tibble() %>%
  mutate(number_1 = as.character(number_1))%>%
  mutate(number_2 = as.character(number_2))%>%
  unite('Item', number_1:number_2, remove = FALSE) %>%
  mutate(Item = as_factor(Item)) %>%
  mutate(Item_specific = ifelse(Analysis_type == "main_con" | Analysis_type == "main_incon", 'Item_spec', 'List_wide' )) %>%
  mutate( Congruency_main = ifelse(Analysis_type == "main_con"| Analysis_type == "MC", "M_Congruent", "M_Incongruent")) %>% 
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject_identifier, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Subject = Subject_identifier, Correct_A = correct_Keyboard_response) 


# Next we define our contrast matrix using the hypr package
StroopCon <- hypr(
  Congruency = (LW_incon_C + LW_incon_I + IS_incon_C + IS_incon_I)/4 ~ (LW_con_C + LW_con_I + IS_con_C + IS_con_I)/4, # main effect Congruence healthy controls
  Listwide = (LW_con_C + LW_incon_C)/2 ~ (LW_con_I + LW_incon_I)/2, # main effect Congruence patients with Parkinson
  Itemspecific = (IS_con_C + IS_incon_C)/2 ~ (IS_con_I + IS_incon_I)/2, # main effect Block healthy controls
  LW_Block = (LW_incon_C - LW_con_C)/2 ~ (LW_incon_I - LW_con_I)/2, # main effect Block patients with Parkinson
  IS_Block = (IS_incon_C - IS_con_C)/2 ~ (IS_incon_I - IS_con_I)/2, # interaction Block Congruence healthy controls
  levels = c("LW_incon_C", "LW_incon_I", "LW_con_C", "LW_con_I", "IS_incon_C", "IS_incon_I", "IS_con_C", "IS_con_I")
)
StroopCon

# Now we create an extra column for the contrasts

Data <- Data %>%
  mutate(Contrast_F = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  ~ "LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent"  ~ "LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  ~ "LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  ~ "LW_incon_I",
    Analysis_type == "main_con" & Congruency == "congruent"  ~ "IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  ~ "IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  ~ "IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  ~ "IS_incon_I",
  )) %>%
  mutate(Contrast_F  = as_factor(Contrast_F)) 

# assign the generated contrast matrix to our new factor
contrasts(Data$Contrast_F) <- contr.hypothesis(StroopCon)

contrasts(Data$Contrast_F)    

# Filter only young participants, remove practice trials and only keep correct trials
Data_y_inducer <- Data %>%
  dplyr::filter(Exp_group == "CY") %>% # only Young
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>%
  dplyr::filter(Trl_type != "Diagnostic")

Data_y_diagnostic <- Data %>%
  dplyr::filter(Exp_group == "CY") %>% # only Young
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>%
  dplyr::filter(Trl_type == "Diagnostic")

# now that the data is prepared let us create our priors
# write down priors - intercept with a mean of 6 and a standard deviation of one seems reasonable exp(5) =150ms exp(7 = 1096ms), values that seem plausible
# for each of the effects we are careful and will assume a normal distribution centered around 0, with a standard deviation of 0.5. I derived at these values by means of prior predictive simulations

prior_informed <- c(
  prior(normal(6, 1), class = Intercept),
  prior(normal(0,1), class = sigma), # truncated normal
  prior(normal(0,0.5), class = b, coef = Contrast_FCongruency),
  prior(normal(0,0.5), class = b, coef = Contrast_FListwide),
  prior(normal(0,0.5), class = b, coef = Contrast_FItemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_FLW_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_FIS_Block),
  prior(normal(0, 3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 3), class = sd, coef = Intercept, group = Item)
)

inducer_young <- brm( RT ~ 1 + Contrast_F + (1|Subject) + (1|Item),
                           data = Data_y_inducer,
                           family = lognormal(),
                           prior = prior_informed,
                           warmup = 2000,
                           iter = 4000#, # 20000 is the limit necessary for bridge sampling
                           #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
)
save(informed_artif_dat, file = "fit_artif_dat.rda")
