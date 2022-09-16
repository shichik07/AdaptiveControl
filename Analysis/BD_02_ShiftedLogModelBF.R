#####
# Author: Julius Kricheldorff
# Analysis Behavioral Data Errors - calculate Bayes_Factors for each Model Term

####
setwd('E:/AdaptiveControl/Data/BehaviorResults')

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
  filter( RT > 200) # Remove data with small RT, otherwise messes with our non-decision parameter



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
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group)) 

Data_y_diagnostic_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_inducer_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type != "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

# Same priors as for the model before

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = sigma, lb = 0), 
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Item)
)


# Prior informed weakly Item Specific
prior_weakly_informed_IS <- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0,0.5), class = sigma, lb = 0), 
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0,0.3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0,0.3), class = sd, coef = Intercept, group = Item)
)

# brmsformula object List Wide
m1_LW <- bf(RT ~ 1  + Contrast_LW + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group)

# brmsformula object Item Specific
m1_IS <- bf(RT ~ 1  + Contrast_IS + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group) 

# okay create a function to pass model parameter
pass_brms = function(save_name, prior, data, model) {
  fit_model <- brm(formula = model,
                   family = shifted_lognormal(),
                   data = data,
                   prior = prior,
                   warmup = 2000,
                   iter = 12000,# 20000 is the limit necessary for bridge sampling
                   cores = 4, seed = 423,
                   control = list(adapt_delta = 0.95),
                   save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                   chains =4)
  
  # save the model
  save(fit_model, file = save_name)
}


##### now let us loop though our models and save the results for the listwide models

for(mods in 1:length(LW_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_LWPC$Contrast_LW) <- contr.hypothesis(eval(parse(text = LW_mods[mods])))
  # define prior
  LW_weakly <- prior_weakly_informed_LW[-(mods+4),]
  # get save name for variable
  save_name <- paste("fit_shifted_inducer_", LW_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = LW_weakly, data = Data_y_inducer_LWPC, model = m1_LW)
} 

# Now the same with the item specific models
for(mods in 1:length(IS_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_ISPC$Contrast_IS) <- contr.hypothesis(eval(parse(text = IS_mods[mods])))
  # define prior
  IS_weakly <- prior_weakly_informed_IS[-(mods+4),]
  # get save name for variable
  save_name <- paste("fit_shifted_inducer_", IS_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = IS_weakly, data = Data_y_inducer_ISPC, model = m1_IS)
  
}

# Lastly is the addition of the ndt distinction by group relevant
temp_hyp_LW <- hypr(LW_formulas,
                 levels = LW_levels)
temp_hyp_IS <- hypr(IS_formulas,
                    levels = IS_levels)
# next add the appropriate variable names
names(temp_hyp_LW) <- LW_contrast_names
names(temp_hyp_IS) <- IS_contrast_names
# Lastly rename the contrast matrix
assign("LW_shifted_log_ndt", temp_hyp_LW)
assign("IS_shifted_log_ndt", temp_hyp_IS)
#assign contrasts to data
contrasts(Data_y_inducer_ISPC$Contrast_IS) <- contr.hypothesis(IS_shifted_log_ndt)
contrasts(Data_y_inducer_LWPC$Contrast_LW) <- contr.hypothesis(LW_shifted_log_ndt)
# models
m1_LW_ndt <- bf(RT ~ 1  + Contrast_LW + (1|Subject) + (1|Item))
m1_IS_ndt <- bf(RT ~ 1  + Contrast_IS + (1|Subject) + (1|Item)) 

# define prior
LW_weakly <- prior_weakly_informed_LW[-(3:4),]
append(LW_weakly, prior(cauchy(0, 5), class = ndt))
IS_weakly <- prior_weakly_informed_IS[-(3:4),]
append(IS_weakly, prior(cauchy(0, 5), class = ndt))

# get save_names
save_name_LW <- paste("fit_shifted_inducer_", "LW_shifted_log_ndt", ".rda", sep = "")
save_name_IS <- paste("fit_shifted_inducer_", "IS_shifted_log_ndt", ".rda", sep = "")

# Run modesl
pass_brms(save_name = save_name_LW, prior = LW_weakly, data = Data_y_inducer_LWPC, model = m1_LW_ndt)
pass_brms(save_name = save_name_IS, prior = IS_weakly, data = Data_y_inducer_ISPC, model = m1_IS_ndt)


##### Next let us get the posterior effects for the new prior

# get posterior to inform prior
load("C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/fit_inducer_model_LW.rda")
load("C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/fit_inducer_model_IS.rda")
ind_sum_LW <- summary(fit_inducer_model_LW)
ind_sum_IS <- summary(fit_inducer_model_IS)

sd_determination <- function(mean, upper95, lower95){
  bound1 <- abs(mean - upper95)
  bound2 <- abs(mean - lower95)
  if (bound1 > bound2){
    return(bound1)
  } else {
    return(bound2)
  }
}

# first we write a prior tibble for the LW and IS effect

prior_tib_LW <- tibble(
  mean = numeric(),
  sd = numeric(),
  var = character()
)

prior_tib_IS <- tibble(
  mean = numeric(),
  sd = numeric(),
  var = character()
)

# Next we add all our variables for the LW effects
prior_tib_LW <-prior_tib_LW %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[1],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[1], ind_sum_LW$fixed$`u-95% CI`[1], ind_sum_LW$fixed$`l-95% CI`[1]),
    var = "Fixed_Intercept"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[2],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[2], ind_sum_LW$fixed$`u-95% CI`[2], ind_sum_LW$fixed$`l-95% CI`[2]),
    var = "Fixed_Contrast_LWCO_Congruency"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[3],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[3], ind_sum_LW$fixed$`u-95% CI`[3], ind_sum_LW$fixed$`l-95% CI`[3]),
    var = "Fixed_Contrast_LWCO_Listwide"
  ) %>%
  add_row(
    mean = 0, # fixed effect for the interaction has to be centered around 0
    sd = sd_determination(0, ind_sum_LW$fixed$`u-95% CI`[4], ind_sum_LW$fixed$`l-95% CI`[4]),
    var = "Fixed_Contrast_LWCO_LW_Block"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[5],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[5], ind_sum_LW$fixed$`u-95% CI`[5], ind_sum_LW$fixed$`l-95% CI`[5]),
    var = "Fixed_Contrast_LWPD_Congruency"
  ) %>%
  add_row(
    mean = 0, # fixed effect for the interaction has to be centered around 0
    sd = sd_determination(0, ind_sum_LW$fixed$`u-95% CI`[7], ind_sum_LW$fixed$`l-95% CI`[7]),
    var = "Fixed_Contrast_LWPD_LW_Block"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[6],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[6], ind_sum_LW$fixed$`u-95% CI`[6], ind_sum_LW$fixed$`l-95% CI`[6]),
    var = "Fixed_Contrast_LWPD_Listwide"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[8],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[8], ind_sum_LW$fixed$`u-95% CI`[8], ind_sum_LW$fixed$`l-95% CI`[8]),
    var = "ndt_Exp_groupCO"
  ) %>%
  add_row(
    mean = ind_sum_LW$fixed$Estimate[9],
    sd = sd_determination(ind_sum_LW$fixed$Estimate[9], ind_sum_LW$fixed$`u-95% CI`[9], ind_sum_LW$fixed$`l-95% CI`[9]),
    var = "ndt_Exp_groupPD"
  ) %>%
  add_row(
    mean = ind_sum_LW$random$Item$Estimate[1],
    sd = sd_determination(ind_sum_LW$random$Item$Estimate[1], ind_sum_LW$random$Item$`u-95% CI`[1], ind_sum_LW$random$Item$`l-95% CI`[1]),
    var = "Item"
  ) %>%
  add_row(
    mean = ind_sum_LW$random$Subject$Estimate[1],
    sd = sd_determination(ind_sum_LW$random$Subject$Estimate[1], ind_sum_LW$random$Subject$`u-95% CI`[1], ind_sum_LW$random$Subject$`l-95% CI`[1]),
    var = "Subject"
  ) %>%
  add_row(
    mean = ind_sum_LW$spec_pars$Estimate[1],
    sd = sd_determination(ind_sum_LW$spec_pars$Estimate[1], ind_sum_LW$spec_pars$`u-95% CI`[1], ind_sum_LW$spec_pars$`l-95% CI`[1]),
    var = "Sigma"
  )


# Next we add all our variables for the IS effects
prior_tib_IS <-prior_tib_IS %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[1],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[1], ind_sum_IS$fixed$`u-95% CI`[1], ind_sum_IS$fixed$`l-95% CI`[1]),
    var = "Fixed_Intercept"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[2],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[2], ind_sum_IS$fixed$`u-95% CI`[2], ind_sum_IS$fixed$`l-95% CI`[2]),
    var = "Fixed_Contrast_ISCO_Congruency"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[3],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[3], ind_sum_IS$fixed$`u-95% CI`[3], ind_sum_IS$fixed$`l-95% CI`[3]),
    var = "Fixed_Contrast_ISCO_Itemspecific"
  ) %>%
  add_row(
    mean = 0, # fixed effect for the interaction has to be centered around 0
    sd = sd_determination(0, ind_sum_IS$fixed$`u-95% CI`[4], ind_sum_IS$fixed$`l-95% CI`[4]),
    var = "Fixed_Contrast_ISCO_IS_Block"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[5],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[5], ind_sum_IS$fixed$`u-95% CI`[5], ind_sum_IS$fixed$`l-95% CI`[5]),
    var = "Fixed_Contrast_ISPD_Congruency"
  ) %>%
  add_row(
    mean = 0, # fixed effect for the interaction has to be centered around 0
    sd = sd_determination(0, ind_sum_IS$fixed$`u-95% CI`[7], ind_sum_IS$fixed$`l-95% CI`[7]),
    var = "Fixed_Contrast_ISPD_IS_Block"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[6],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[6], ind_sum_IS$fixed$`u-95% CI`[6], ind_sum_IS$fixed$`l-95% CI`[6]),
    var = "Fixed_Contrast_ISPD_Itemspecific"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[8],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[8], ind_sum_IS$fixed$`u-95% CI`[8], ind_sum_IS$fixed$`l-95% CI`[8]),
    var = "ndt_Exp_groupCO"
  ) %>%
  add_row(
    mean = ind_sum_IS$fixed$Estimate[9],
    sd = sd_determination(ind_sum_IS$fixed$Estimate[9], ind_sum_IS$fixed$`u-95% CI`[9], ind_sum_IS$fixed$`l-95% CI`[9]),
    var = "ndt_Exp_groupPD"
  ) %>%
  add_row(
    mean = ind_sum_IS$random$Item$Estimate[1],
    sd = sd_determination(ind_sum_IS$random$Item$Estimate[1], ind_sum_IS$random$Item$`u-95% CI`[1], ind_sum_IS$random$Item$`l-95% CI`[1]),
    var = "Item"
  ) %>%
  add_row(
    mean = ind_sum_IS$random$Subject$Estimate[1],
    sd = sd_determination(ind_sum_IS$random$Subject$Estimate[1], ind_sum_IS$random$Subject$`u-95% CI`[1], ind_sum_IS$random$Subject$`l-95% CI`[1]),
    var = "Subject"
  ) %>%
  add_row(
    mean = ind_sum_IS$spec_pars$Estimate[1],
    sd = sd_determination(ind_sum_IS$spec_pars$Estimate[1], ind_sum_IS$spec_pars$`u-95% CI`[1], ind_sum_IS$spec_pars$`l-95% CI`[1]),
    var = "Sigma"
  )


# next we get this stuff into stan variables for the LW effect
stan_vars_LW <- stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Intercept", ]$mean ,name = 'Fixed_Intercept_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Intercept", ]$sd ,name = 'Fixed_Intercept_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_Congruency", ]$mean ,name = 'Fixed_Contrast_LWCO_Congruency_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_Congruency", ]$sd ,name = 'Fixed_Contrast_LWCO_Congruency_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_Listwide", ]$mean ,name = 'Fixed_Contrast_LWCO_Listwide_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_Listwide", ]$sd ,name = 'Fixed_Contrast_LWCO_Listwide_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_LW_Block", ]$mean ,name = 'Fixed_Contrast_LWCO_LW_Block_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWCO_LW_Block", ]$sd ,name = 'Fixed_Contrast_LWCO_LW_Block_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_Congruency", ]$mean ,name = 'Fixed_Contrast_LWPD_Congruency_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_Congruency", ]$sd ,name = 'Fixed_Contrast_LWPD_Congruency_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_Listwide", ]$mean ,name = 'Fixed_Contrast_LWPD_Listwide_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_Listwide", ]$sd ,name = 'Fixed_Contrast_LWPD_Listwide_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_LW_Block", ]$mean ,name = 'Fixed_Contrast_LWPD_LW_Block_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Fixed_Contrast_LWPD_LW_Block", ]$sd ,name = 'Fixed_Contrast_LWPD_LW_Block_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Item", ]$mean ,name = 'Item_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Item", ]$sd ,name = 'Item_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Subject", ]$mean ,name = 'Subject_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Subject", ]$sd ,name = 'Subject_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Sigma", ]$mean ,name = 'Sigma_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "Sigma", ]$sd ,name = 'Sigma_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "ndt_Exp_groupCO", ]$mean ,name = 'ndt_Exp_groupCO_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "ndt_Exp_groupCO", ]$sd ,name = 'ndt_Exp_groupCO_sd') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "ndt_Exp_groupPD", ]$mean ,name = 'ndt_Exp_groupPD_mean') +
  stanvar(prior_tib_LW[prior_tib_LW$var == "ndt_Exp_groupPD", ]$sd ,name = 'ndt_Exp_groupPD_sd')

# next we get this stuff into stan variables for the LW effect
stan_vars_IS <- stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Intercept", ]$mean ,name = 'Fixed_Intercept_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Intercept", ]$sd ,name = 'Fixed_Intercept_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_Congruency", ]$mean ,name = 'Fixed_Contrast_ISCO_Congruency_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_Congruency", ]$sd ,name = 'Fixed_Contrast_ISCO_Congruency_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_Itemspecific", ]$mean ,name = 'Fixed_Contrast_ISCO_Itemspecific_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_Itemspecific", ]$sd ,name = 'Fixed_Contrast_ISCO_Itemspecific_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_IS_Block", ]$mean ,name = 'Fixed_Contrast_ISCO_IS_Block_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISCO_IS_Block", ]$sd ,name = 'Fixed_Contrast_ISCO_IS_Block_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_Congruency", ]$mean ,name = 'Fixed_Contrast_ISPD_Congruency_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_Congruency", ]$sd ,name = 'Fixed_Contrast_ISPD_Congruency_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_Itemspecific", ]$mean ,name = 'Fixed_Contrast_ISPD_Itemspecific_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_Itemspecific", ]$sd ,name = 'Fixed_Contrast_ISPD_Itemspecific_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_IS_Block", ]$mean ,name = 'Fixed_Contrast_ISPD_IS_Block_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Fixed_Contrast_ISPD_IS_Block", ]$sd ,name = 'Fixed_Contrast_ISPD_IS_Block_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Item", ]$mean ,name = 'Item_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Item", ]$sd ,name = 'Item_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Subject", ]$mean ,name = 'Subject_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Subject", ]$sd ,name = 'Subject_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Sigma", ]$mean ,name = 'Sigma_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "Sigma", ]$sd ,name = 'Sigma_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "ndt_Exp_groupCO", ]$mean ,name = 'ndt_Exp_groupCO_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "ndt_Exp_groupCO", ]$sd ,name = 'ndt_Exp_groupCO_sd') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "ndt_Exp_groupPD", ]$mean ,name = 'ndt_Exp_groupPD_mean') +
  stanvar(prior_tib_IS[prior_tib_IS$var == "ndt_Exp_groupPD", ]$sd ,name = 'ndt_Exp_groupPD_sd')

# now put in our variables into the prior for LW
prior_informed_LW <- c(
  prior(normal(Fixed_Intercept_mean, Fixed_Intercept_sd), class = Intercept), 
  prior(normal(Fixed_Contrast_LWCO_Congruency_mean, Fixed_Contrast_LWCO_Congruency_sd), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_LWCO_Listwide_mean, Fixed_Contrast_LWCO_Listwide_sd), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(Fixed_Contrast_LWCO_LW_Block_mean, Fixed_Contrast_LWCO_LW_Block_sd), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(Fixed_Contrast_LWPD_Congruency_mean, Fixed_Contrast_LWPD_Congruency_sd), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_LWPD_Listwide_mean, Fixed_Contrast_LWPD_Listwide_sd), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(Fixed_Contrast_LWPD_LW_Block_mean, Fixed_Contrast_LWPD_LW_Block_sd), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(Item_mean, Item_sd), class = sd, coef = Intercept, group = Item),
  prior(normal(Subject_mean, Subject_sd), class = sd, coef = Intercept, group = Subject),
  prior(normal(Sigma_mean, Sigma_sd), class = sigma),
  prior(normal(ndt_Exp_groupCO_mean, ndt_Exp_groupCO_sd), class = b, coef = Exp_groupCO, dpar = ndt), 
  prior(normal(ndt_Exp_groupPD_mean, ndt_Exp_groupPD_sd), class = b, coef = Exp_groupPD, dpar = ndt)
)

# now put in our variables into the prior for IS
prior_informed_IS <- c(
  prior(normal(Fixed_Intercept_mean, Fixed_Intercept_sd), class = Intercept), 
  prior(normal(Fixed_Contrast_ISCO_Congruency_mean, Fixed_Contrast_ISCO_Congruency_sd), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_ISCO_Itemspecific_mean, Fixed_Contrast_ISCO_Itemspecific_sd), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(Fixed_Contrast_ISCO_IS_Block_mean, Fixed_Contrast_ISCO_IS_Block_sd), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(Fixed_Contrast_ISPD_Congruency_mean, Fixed_Contrast_ISPD_Congruency_sd), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_ISPD_Itemspecific_mean, Fixed_Contrast_ISPD_Itemspecific_sd), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(Fixed_Contrast_ISPD_IS_Block_mean, Fixed_Contrast_ISPD_IS_Block_sd), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(Item_mean, Item_sd), class = sd, coef = Intercept, group = Item),
  prior(normal(Subject_mean, Subject_sd), class = sd, coef = Intercept, group = Subject),
  prior(normal(Sigma_mean, Sigma_sd), class = sigma),
  prior(normal(ndt_Exp_groupCO_mean, ndt_Exp_groupCO_sd), class = b, coef = Exp_groupCO, dpar = ndt), 
  prior(normal(ndt_Exp_groupPD_mean, ndt_Exp_groupPD_sd), class = b, coef = Exp_groupPD, dpar = ndt)
)
# okay create a function to pass model parameter for the informed decisions
pass_brms_informed = function(save_name, prior, data, model, stanvariables) {
  fit_model <- brm(formula = model,
                   family = shifted_lognormal(),
                   data = data,
                   prior = prior,
                   warmup = 2000,
                   iter = 12000,# 20000 is the limit necessary for bridge sampling
                   cores = 4, seed = 423,
                   control = list(adapt_delta = 0.95),
                   save_pars = save_pars(all = TRUE),# must be set to true for bridgesampling
                   stanvars = stanvariables,
                   chains =4)
  
  # save the model
  save(fit_model, file = save_name)
}

##### now let us loop though our models and save the results for the listwide models

for(mods in 1:length(LW_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_diagnostic_LWPC$Contrast_LW) <- contr.hypothesis(eval(parse(text = LW_mods[mods])))
  # define prior
  LW_informed <- prior_informed_LW[-(mods+1),]
  # get save name for variable
  save_name <- paste("fit_shifted_diagnostic_", LW_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_informed(save_name = save_name, prior =  LW_informed, 
                     data = Data_y_diagnostic_LWPC, model = m1_LW, stanvariables = stan_vars_LW)
} 

# Now the same with the item specific models
for(mods in 1:length(IS_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_diagnostic_ISPC$Contrast_IS) <- contr.hypothesis(eval(parse(text = IS_mods[mods])))
  # define prior
  IS_informed <- prior_informed_IS[-(mods+1),]
  # get save name for variable
  save_name <- paste("fit_shifted_diagnostic_", IS_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_informed(save_name = save_name, prior = IS_informed, 
                     data = Data_y_diagnostic_ISPC, model = m1_IS, stanvariables = stan_vars_IS)
  
} 

# ANd the same for the ndt models

contrasts(Data_y_diagnostic_ISPC$Contrast_IS) <- contr.hypothesis(IS_shifted_log_ndt)
contrasts(Data_y_diagnostic_LWPC$Contrast_LW) <- contr.hypothesis(LW_shifted_log_ndt)

# define prior
LW_informed <- prior_informed_LW[-(11:12),]
append(LW_informed, prior(normal(5.3, 0.5), class = ndt))
IS_informed <- prior_informed_IS[-(11:12),]
append(IS_informed, prior(normal(5.3, 0.5), class = ndt))

# get save_names
save_name_LW <- paste("fit_shifted_diagnostic_", "LW_shifted_log_ndt", ".rda", sep = "")
save_name_IS <- paste("fit_shifted_diagnostic_", "IS_shifted_log_ndt", ".rda", sep = "")

# Run modesl
pass_brms_informed(save_name = save_name_LW, prior = LW_informed,
                   data = Data_y_diagnostic_LWPC, model = m1_LW_ndt, stanvariables = stan_vars_LW)
pass_brms_informed(save_name = save_name_IS, prior = IS_informed,
                   data = Data_y_diagnostic_ISPC, model = m1_IS_ndt, stanvariables = stan_vars_IS)
