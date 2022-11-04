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
library(broom) # to access and manipulate summary object


#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32946)

#### Load the  behavioral data ####
# not we use not participant _16 because he slept in during the experiment 

# load the clinical data
#load the clinical data first
Clinic <- read.table(file = "participants.tsv", header = TRUE)
Clinic <- as_tibble(Clinic)

path1 = "Z:/JuliusKricheldorff/AdaptiveControl/BehavioralData"
files = list.files(path = path1, pattern = "AdapCon_CO", full.names = TRUE, recursive = FALSE)
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

# Contrasts only for the list-wide effect only
StroopCon_LW <- hypr(
  CO_Congruency = (CO_LW_incon_C + CO_LW_con_C)/2 ~ (CO_LW_con_I + CO_LW_incon_I)/2, # main effect Congruence healthy controls
  CO_Listwide = (CO_LW_con_C + CO_LW_con_I)/2 ~ (CO_LW_incon_C + CO_LW_incon_I)/2, # main effect Block Listwide control
  CO_LW_Block = (CO_LW_incon_I - CO_LW_incon_C)/2 ~ (CO_LW_con_I - CO_LW_con_C)/2, # interaction listwide effect
  levels = c("CO_LW_con_C", "CO_LW_con_I","CO_LW_incon_I", "CO_LW_incon_C")
)
StroopCon_LW

# Contrasts for the Item-specific effects only - check that later (might be more difficult)
StroopCon_IS <- hypr(
  CO_Congruency = (CO_IS_incon_C + CO_IS_con_C)/2 ~ (CO_IS_con_I + CO_IS_incon_I)/2, # main effect Congruence healthy controls
  CO_Itemspecific = (CO_IS_con_C + CO_IS_con_I)/2 ~ (CO_IS_incon_C + CO_IS_incon_I)/2, # main effect Block Listwide control
  CO_IS_Block = (CO_IS_incon_I - CO_IS_incon_C)/2 ~ (CO_IS_con_I - CO_IS_con_C)/2, # interaction listwide effect
  levels = c("CO_IS_incon_I", "CO_IS_con_I", "CO_IS_con_C", "CO_IS_incon_C")
)
StroopCon_IS


# Now we create an extra column for the contrasts for the List wide effect
Data <- Data %>%
  mutate(Contrast_LW = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "CO"  ~ "CO_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_LW_incon_I",
  )) %>%
  mutate(Contrast_LW  = as_factor(Contrast_LW)) 

# assign the generated contrast matrix to the List Wide Factor
contrasts(Data$Contrast_LW) <- contr.hypothesis(StroopCon_LW)

contrasts(Data$Contrast_LW)    

# Now we create an extra column for the contrasts for the Item Specific effect
Data <- Data %>%
  mutate(Contrast_IS = case_when(
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_incon_I",
  )) %>%
  mutate(Contrast_IS  = as_factor(Contrast_IS))

# assign the generated contrast matrix to the Item Specific Factor
contrasts(Data$Contrast_IS) <- contr.hypothesis(StroopCon_IS)

contrasts(Data$Contrast_IS)    




## Get seperate datasets
Data_y_inducer_LWPC <- Data %>%
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic") %>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

Data_y_diagnostic_LWPC <- Data %>%
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic") %>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))


Data_y_inducer_ISPC <- Data %>%
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Congruency  = as_factor(Congruency)) %>%
  mutate(Analysis_type = as_factor(Analysis_type))

# brmsformula object List Wide
m2_CO_LW <- bf(RT ~ 1 + Contrast_LW + (Contrast_LW|Subject) + (1|Item)) 
a <- get_prior(data = Data_y_inducer_LWPC, m2_CO_LW, family = shifted_lognormal()) 

## Weakly informative prior for the listwide and item specific effect

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(normal(6.5, 0.5), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = sigma, lb = 0), 
  prior(normal(5.3, 0.5), class = ndt),
  prior(lkj(1), class = cor), # use default prior for the correlation 
  prior(lkj(1), class = cor, group = Subject),
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO patients
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Contrast_LWCO_Congruency, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Contrast_LWCO_Listwide, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Contrast_LWCO_LW_Block, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Item)
)


fit_ShiftedLog_COOnly_inducer_LW <- brm(formula = m2_CO_LW ,
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

save(fit_ShiftedLog_COOnly_inducer_LW, file = "fit_ShiftedLog_COOnly_inducer_LW.rda")

# get posterior to inform prior
load("fit_ShiftedLog_COOnly_inducer_LW.rda")
ind_sum <- summary(fit_ShiftedLog_COOnly_inducer_LW)

sd_determination <- function(mean, upper95, lower95){
  bound1 <- abs(mean - upper95)
  bound2 <- abs(mean - lower95)
  if (bound1 > bound2){
    return(bound1)
  } else {
    return(bound2)
  }
}

# first we write a prior tibble

prior_tib <- tibble(
  mean = numeric(),
  sd = numeric(),
  var = character()
)

# Next we add all our variables
prior_tib <-prior_tib %>%
  add_row(
    mean = ind_sum$fixed$Estimate[1],
    sd = sd_determination(ind_sum$fixed$Estimate[1], ind_sum$fixed$`u-95% CI`[1], ind_sum$fixed$`l-95% CI`[1]),
    var = "Fixed_Intercept"
  ) %>%
  add_row(
    mean = ind_sum$spec_pars$Estimate[1],
    sd = sd_determination(ind_sum$spec_pars$Estimate[1], ind_sum$spec_pars$`u-95% CI`[1], ind_sum$spec_pars$`l-95% CI`[1]),
    var = "sigma"
  ) %>%
  add_row(
    mean = ind_sum$spec_pars$Estimate[2],
    sd = sd_determination(ind_sum$spec_pars$Estimate[2], ind_sum$spec_pars$`u-95% CI`[2], ind_sum$spec_pars$`l-95% CI`[2]),
    var = "ndt"
  ) %>%
  add_row(
    mean = ind_sum$fixed$Estimate[2],
    sd = sd_determination(ind_sum$fixed$Estimate[2], ind_sum$fixed$`u-95% CI`[2], ind_sum$fixed$`l-95% CI`[2]),
    var = "Fixed_Contrast_LWCO_Congruency"
  ) %>%
  add_row(
    mean = 0, # fixed effect for the interaction has to be centered around 0
    sd = sd_determination(0, ind_sum$fixed$`u-95% CI`[4], ind_sum$fixed$`l-95% CI`[4]),
    var = "Fixed_Contrast_LWCO_LW_Block"
  ) %>%
  add_row(
    mean = ind_sum$fixed$Estimate[3],
    sd = sd_determination(ind_sum$fixed$Estimate[3], ind_sum$fixed$`u-95% CI`[3], ind_sum$fixed$`l-95% CI`[3]),
    var = "Fixed_Contrast_LWCO_Listwide"
  ) %>%
  add_row(
    mean = ind_sum$random$Subject$Estimate[1],
    sd = sd_determination(ind_sum$random$Subject$Estimate[1], ind_sum$random$Subject$`u-95% CI`[1], ind_sum$random$Subject$`l-95% CI`[1]),
    var = "Random_intercept"
  ) %>%
  add_row(
    mean = ind_sum$random$Subject$Estimate[2],
    sd = sd_determination(ind_sum$random$Subject$Estimate[2], ind_sum$random$Subject$`u-95% CI`[2], ind_sum$random$Subject$`l-95% CI`[2]),
    var = "Random_Contrast_LWCO_Congruency"
  ) %>%
  add_row(
    mean = ind_sum$random$Subject$Estimate[3],
    sd = sd_determination(ind_sum$random$Subject$Estimate[3], ind_sum$random$Subject$`u-95% CI`[3], ind_sum$random$Subject$`l-95% CI`[3]),
    var = "Random_Contrast_LWCO_Listwide"
  ) %>%
  add_row(
    mean = ind_sum$random$Subject$Estimate[4],
    sd = sd_determination(ind_sum$random$Subject$Estimate[4], ind_sum$random$Subject$`u-95% CI`[4], ind_sum$random$Subject$`l-95% CI`[4]),
    var = "Random_Contrast_LWCO_LW_Block"
  ) %>%
  add_row(
    mean = ind_sum$random$Item$Estimate[1],
    sd = sd_determination(ind_sum$random$Item$Estimate[1], ind_sum$random$Item$`u-95% CI`[1], ind_sum$random$Item$`l-95% CI`[1]),
    var = "Item"
  ) 

# next we get this stuff into stan variables

stan_vars <- stanvar(prior_tib[prior_tib$var == "Fixed_Intercept", ]$mean ,name = 'Fixed_Intercept_mean') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Intercept", ]$sd ,name = 'Fixed_Intercept_sd') +
  stanvar(prior_tib[prior_tib$var == "sigma", ]$mean ,name = 'sigma_mean') +
  stanvar(prior_tib[prior_tib$var == "sigma", ]$sd ,name = 'sigma_sd') +
  stanvar(prior_tib[prior_tib$var == "ndt", ]$mean ,name = 'ndt_mean') +
  stanvar(prior_tib[prior_tib$var == "ndt", ]$sd ,name = 'ndt_sd') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_Congruency", ]$mean ,name = 'Fixed_Contrast_LWCO_Congruency_mean') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_Congruency", ]$sd ,name = 'Fixed_Contrast_LWCO_Congruency_sd') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_Listwide", ]$mean ,name = 'Fixed_Contrast_LWCO_Listwide_mean') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_Listwide", ]$sd ,name = 'Fixed_Contrast_LWCO_Listwide_sd') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_LW_Block", ]$mean ,name = 'Fixed_Contrast_LWCO_LW_Block_mean') +
  stanvar(prior_tib[prior_tib$var == "Fixed_Contrast_LWCO_LW_Block", ]$sd ,name = 'Fixed_Contrast_LWCO_LW_Block_sd') +
  stanvar(prior_tib[prior_tib$var == "Random_intercept", ]$mean ,name = 'Random_intercept_mean') +
  stanvar(prior_tib[prior_tib$var == "Random_intercept", ]$sd ,name = 'Random_intercept_sd') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LWCO_Congruency", ]$mean ,name = 'Random_Contrast_LWCO_Congruency_mean') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LWCO_Congruency", ]$sd ,name = 'Random_Contrast_LWCO_Congruency_sd') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LWCO_Listwide", ]$mean ,name = 'Random_Contrast_LWCO_Listwide_mean') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LWCO_Listwide", ]$sd ,name = 'Random_Contrast_LWCO_Listwide_sd') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LW_Block", ]$mean ,name = 'Random_Contrast_LW_Block_mean') +
  stanvar(prior_tib[prior_tib$var == "Random_Contrast_LW_Block", ]$sd ,name = 'Random_Contrast_LW_Block_sd') +
  stanvar(prior_tib[prior_tib$var == "Item", ]$mean ,name = 'Item_mean') +
  stanvar(prior_tib[prior_tib$var == "Item", ]$sd ,name = 'Item_sd')

# now put in our variables into the prior
# now put in our variables into the prior
prior_informed_LW <- c(
  prior(normal(Fixed_Intercept_mean, Fixed_Intercept_sd), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(sigma_mean, sigma_sd), class = sigma), 
  prior(normal(ndt_mean, ndt_sd), class = ndt),
  prior(lkj(1), class = cor), # use default prior for the correlation 
  prior(lkj(1), class = cor, group = Subject),
  prior(normal(Fixed_Contrast_LWCO_Congruency_mean, Fixed_Contrast_LWCO_Congruency_sd), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO patients
  prior(normal(Fixed_Contrast_LWCO_Listwide_mean, Fixed_Contrast_LWCO_Listwide_sd), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(Fixed_Contrast_LWCO_LW_Block_mean, Fixed_Contrast_LWCO_LW_Block_sd), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(Random_intercept_mean, Random_intercept_sd), class = sd, coef = Intercept, group = Subject),
  prior(normal(Random_Contrast_LWCO_Congruency_mean, Random_Contrast_LWCO_Congruency_sd), class = sd, coef = Contrast_LWCO_Congruency, group = Subject),
  prior(normal(Random_Contrast_LWCO_Listwide_mean, Random_Contrast_LWCO_Listwide_sd), class = sd, coef = Contrast_LWCO_Listwide, group = Subject),
  prior(normal(Random_Contrast_LW_Block_mean, Random_Contrast_LW_Block_sd), class = sd, coef = Contrast_LWCO_LW_Block, group = Subject),
  prior(normal(Item_mean, Item_sd), class = sd, coef = Intercept, group = Item)
)

fit_ShiftedLog_COOnly_diagnostic_LW <- brm(formula = m2_CO_LW ,
                                           family = shifted_lognormal(),
                                           data = Data_y_diagnostic_LWPC,
                                           prior = prior_informed_LW,
                                           warmup = 2000,
                                           iter = 12000,# 20000 is the limit necessary for bridge sampling
                                           cores = 4, seed = 423,
                                           control = list(adapt_delta = 0.95),
                                           stanvars=stan_vars,
                                           #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                           chains =4
)

save(fit_ShiftedLog_COOnly_diagnostic_LW, file = "fit_ShiftedLog_COOnly_diagnostic_LW.rda")
