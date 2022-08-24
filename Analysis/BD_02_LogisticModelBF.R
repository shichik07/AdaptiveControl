#####
# Author: Julius Kricheldorff
# Analysis Behavioral Data Errors - calculate Bayes_Factors for each Model Term

####
setwd('C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data')

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)
library(rstudioapi)

#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32946)

#Load the data 
filename <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/All_subjects.csv"
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
  dplyr::rename(RT = response_time_Keyboard_response, Subject = Subject_identifier, Correct_A = correct_Keyboard_response) %>%
  mutate(Error = case_when( # we define accuracy in terms of errors, easier to interpret
    Correct_A == 1 ~ 0,
    Correct_A == 0 ~1
  ))

# Define formulas so we can loop through them
LW_formulas <- c(
  formula((CO_LW_incon_C + CO_LW_con_C)/2 ~ (CO_LW_con_I + CO_LW_incon_I)/2), # main effect Congruence healthy controls
  formula((CO_LW_con_C + CO_LW_con_I)/2 ~ (CO_LW_incon_C + CO_LW_incon_I)/2), # main effect Block Listwide control
  formula((CO_LW_incon_I - CO_LW_incon_C)/2 ~ (CO_LW_con_I - CO_LW_con_C)/2), # interaction listwide effect
  formula((PD_LW_incon_C + PD_LW_con_C)/2 ~ (PD_LW_con_I + PD_LW_incon_I)/2), # main effect Congruence healthy controls
  formula((PD_LW_con_C + PD_LW_con_I)/2 ~ (PD_LW_incon_C + PD_LW_incon_I)/2), # main effect Block Listwide control
  formula((PD_LW_incon_I - PD_LW_incon_C)/2 ~ (PD_LW_con_I - PD_LW_con_C)/2)
)

# contrast names seperately
LW_contrast_names <- c(
  "CO_Congruency",
  "CO_Listwide",
  "CO_LW_Block",
  "PD_Congruency",
  "PD_Listwide",
  "PD_LW_Block"
)

# LW levels
LW_levels <- c("CO_LW_con_C", "CO_LW_con_I", "CO_LW_incon_I", "CO_LW_incon_C", 
            "PD_LW_incon_I", "PD_LW_incon_C", "PD_LW_con_C", "PD_LW_con_I")

# same for the IS effect
IS_formulas <-c(
  formula((CO_IS_incon_C + CO_IS_con_C)/2 ~ (CO_IS_con_I + CO_IS_incon_I)/2), # main effect Congruence healthy controls
  formula((CO_IS_con_C + CO_IS_con_I)/2 ~ (CO_IS_incon_C + CO_IS_incon_I)/2), # main effect Block Listwide control
  formula((CO_IS_incon_I - CO_IS_incon_C)/2 ~ (CO_IS_con_I - CO_IS_con_C)/2), # interaction listwide effect
  formula((PD_IS_incon_C + PD_IS_con_C)/2 ~ (PD_IS_con_I + PD_IS_incon_I)/2), # main effect Congruence healthy controls
  formula((PD_IS_con_C + PD_IS_con_I)/2 ~ (PD_IS_incon_C + PD_IS_incon_I)/2), # main effect Block Listwide control
  formula((PD_IS_incon_I - PD_IS_incon_C)/2 ~ (PD_IS_con_I - PD_IS_con_C)/2) # interaction listwide effect
)

# contrast names seperately
IS_contrast_names <- c(
  "CO_Congruency",
  "CO_Itemspecific",
  "CO_IS_Block",
  "PD_Congruency",
  "PD_Itemspecific",
  "PD_IS_Block"
)

# IS levels
IS_levels <- c("CO_IS_incon_I", "CO_IS_con_I", "CO_IS_con_C", "CO_IS_incon_C", 
                        "PD_IS_con_C", "PD_IS_con_I", "PD_IS_incon_I", "PD_IS_incon_C")

LW_mods <- c("LW_min_CO_Congruency", 
             "LW_min_CO_Listwide",
             "LW_min_CO_LW_Block",
             "LW_min_PD_Congruency",
             "LW_min_PD_Listwide",
             "LW_min_PD_LW_Block")

IS_mods <- c("IS_min_CO_Congruency", 
             "IS_min_CO_Listwide",
             "IS_min_CO_LW_Block",
             "IS_min_PD_Congruency",
             "IS_min_PD_Listwide",
             "IS_min_PD_LW_Block")

# First let us get the listwide model contrasts
for(mods in 1:length(LW_mods)){
  # first create the contrast matrix
  temp_mat <- hypr(LW_formulas[-mods],
                   levels = LW_levels)
  # next add the appropriate variable names
  names(temp_mat) <- LW_contrast_names[-mods]
  
  # Lastly rename the contrast matrix
  assign(LW_mods[mods], temp_mat)
}

# Now the same for the Itemspecific model contrasts
for(mods in 1:length(IS_mods)){
  # first create the contrast matrix
  temp_mat <- hypr(IS_formulas[-mods],
                   levels = IS_levels)
  # next add the appropriate variable names
  names(temp_mat) <- IS_contrast_names[-mods]
  
  # Lastly rename the contrast matrix
  assign(IS_mods[mods], temp_mat)
}

# Now we create an extra column for the contrasts for the List wide effect 
Data <- Data %>%
  mutate(Contrast_LW = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "CO"  ~ "CO_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_LW_incon_I",
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "PD"  ~ "PD_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_LW_incon_I",
  )) %>%
  mutate(Contrast_LW  = as_factor(Contrast_LW)) 

# Now we create an extra column for the contrasts for the Item Specific effect
Data <- Data %>%
  mutate(Contrast_IS = case_when(
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_incon_I",
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_incon_I",
  )) %>%
  mutate(Contrast_IS  = as_factor(Contrast_IS))


# Filter only young participants, remove practice trials and only keep correct trials
Data_y_inducer_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group)) 

Data_y_diagnostic_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_inducer_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

# Same priors as for the model before

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(normal(-1.3, 1.5), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0 ,1.5), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0 ,1.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0 ,1.5), class = sd, coef = Intercept, group = Item)
)

# Prior informed weakly Item Specific
prior_weakly_informed_IS <- c(
  prior(normal(-1.3, 1.5), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0 ,1.5), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0 ,1.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0 ,1.5), class = sd, coef = Intercept, group = Item)
)

# brmsformula object List Wide
m1_LW <- bf(Error ~ 1  + Contrast_LW + (1|Subject) + (1|Item))

# brmsformula object Item Specific
m1_IS <- bf(Error ~ 1  + Contrast_IS + (1|Subject) + (1|Item)) 

##### now let us loop though our models and save the results for the listwide models

for(mods in 1:length(LW_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_LWPC$Contrast_LW) <- contr.hypothesis(eval(parse(text = LW_mods[mods])))
  # define prior
  LW_weakly <- prior_weakly_informed_LW[-(mods+1),]
  # get save name for variable
  save_name <- paste("fit_inducer_", LW_mods[mods], ".rda", sep = "")
  # fit the model
  LW_inducer_model <- brm(formula = m1_LW,
                                    family = bernoulli(link = logit),
                                    data = Data_y_inducer_LWPC,
                                    prior = LW_weakly,
                                    warmup = 2000,
                                    iter = 22000,# 20000 is the limit necessary for bridge sampling
                                    cores = 4, seed = 423,
                                    control = list(adapt_delta = 0.95),
                                    save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                    chains =4
  )
  # save the model
  save(LW_inducer_model, file = save_name)
  # just in case remove model data
  rm(LW_inducer_model)
} 

# Now the same with the item specific models
for(mods in 1:length(IS_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_ISPC$Contrast_IS) <- contr.hypothesis(eval(parse(text = IS_mods[mods])))
  # define prior
  IS_weakly <- prior_weakly_informed_IS[-(mods+1),]
  # get save name for variable
  save_name <- paste("fit_inducer_", IS_mods[mods], ".rda", sep = "")
  # fit the model
  IS_inducer_model <- brm(formula = m1_IS,
                          family = bernoulli(link = logit),
                          data = Data_y_inducer_ISPC,
                          prior = IS_weakly,
                          warmup = 200,
                          iter = 220,# 20000 is the limit necessary for bridge sampling
                          cores = 4, seed = 423,
                          control = list(adapt_delta = 0.95),
                          save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                          chains =4
  )
  # save the model
  save(IS_inducer_model, file = save_name)
  # just in case remove model data
  rm(IS_inducer_model)
} 



