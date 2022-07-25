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
  filter( RT >= 200) # filter RTs smaller than 200ms 20 cases for the whole dataset





# Next we define our contrast matrix using the hypr package
StroopCon <- hypr(
  CO_Congruency = (CO_LW_incon_C + CO_LW_con_C + CO_IS_incon_C + CO_IS_con_C)/4 ~ (CO_LW_con_I + CO_LW_incon_I + CO_IS_con_I + CO_IS_incon_I)/4, # main effect Congruence healthy controls
  CO_Listwide = (CO_LW_con_C + CO_LW_con_I)/2 ~ (CO_LW_incon_C + CO_LW_incon_I)/2, # main effect Block Listwide control
  CO_Itemspecific = (CO_IS_con_C + CO_IS_con_I)/2 ~ (CO_IS_incon_C + CO_IS_incon_I)/2, # main effect Itemspecfic
  CO_LW_Block = (CO_LW_incon_I - CO_LW_incon_C)/2 ~ (CO_LW_con_I - CO_LW_con_C)/2, # interaction listwide effect
  CO_IS_Block = (CO_IS_incon_I - CO_IS_incon_C)/2 ~ (CO_IS_con_I - CO_IS_con_C)/2, # interaction itemspecfic effect
  PD_Congruency = (PD_LW_incon_C + PD_LW_con_C + PD_IS_incon_C + PD_IS_con_C)/4 ~ (PD_LW_con_I + PD_LW_incon_I + PD_IS_con_I + PD_IS_incon_I)/4, # main effect Congruence healthy controls
  PD_Listwide = (PD_LW_con_C + PD_LW_con_I)/2 ~ (PD_LW_incon_C + PD_LW_incon_I)/2, # main effect Block Listwide control
  PD_Itemspecific = (PD_IS_con_C + PD_IS_con_I)/2 ~ (PD_IS_incon_C + PD_IS_incon_I)/2, # main effect Itemspecfic
  PD_LW_Block = (PD_LW_incon_I - PD_LW_incon_C)/2 ~ (PD_LW_con_I - PD_LW_con_C)/2, # interaction listwide effect
  PD_IS_Block = (PD_IS_incon_I - PD_IS_incon_C)/2 ~ (PD_IS_con_I - PD_IS_con_C)/2, # interaction itemspecfic effect
  levels = c("CO_LW_con_C", "CO_LW_con_I", "CO_LW_incon_I", "CO_LW_incon_C", "CO_IS_incon_I", "CO_IS_con_I", "CO_IS_con_C", "CO_IS_incon_C",
             "PD_LW_con_C", "PD_LW_con_I", "PD_LW_incon_I", "PD_LW_incon_C", "PD_IS_incon_I", "PD_IS_con_I", "PD_IS_con_C", "PD_IS_incon_C")
)
StroopCon

# Contrasts only for the list-wide effect only
StroopCon_LW <- hypr(
  CO_Congruency = (CO_LW_incon_C + CO_LW_con_C)/2 ~ (CO_LW_con_I + CO_LW_incon_I)/2, # main effect Congruence healthy controls
  CO_Listwide = (CO_LW_con_C + CO_LW_con_I)/2 ~ (CO_LW_incon_C + CO_LW_incon_I)/2, # main effect Block Listwide control
  CO_LW_Block = (CO_LW_incon_I - CO_LW_incon_C)/2 ~ (CO_LW_con_I - CO_LW_con_C)/2, # interaction listwide effect
  PD_Congruency = (PD_LW_incon_C + PD_LW_con_C)/2 ~ (PD_LW_con_I + PD_LW_incon_I)/2, # main effect Congruence healthy controls
  PD_Listwide = (PD_LW_con_C + PD_LW_con_I)/2 ~ (PD_LW_incon_C + PD_LW_incon_I)/2, # main effect Block Listwide control
  PD_LW_Block = (PD_LW_incon_I - PD_LW_incon_C)/2 ~ (PD_LW_con_I - PD_LW_con_C)/2, # interaction listwide effect
  levels = c("CO_LW_con_C", "CO_LW_con_I", "CO_LW_incon_I", "CO_LW_incon_C", 
             "PD_LW_incon_I", "PD_LW_incon_C", "PD_LW_con_C", "PD_LW_con_I")
)
StroopCon_LW

# Contrasts for the Item-specific effects only - check that later (might be more difficult)
StroopCon_IS <- hypr(
  CO_Congruency = (CO_IS_incon_C + CO_IS_con_C)/2 ~ (CO_IS_con_I + CO_IS_incon_I)/2, # main effect Congruence healthy controls
  CO_Itemspecific = (CO_IS_con_C + CO_IS_con_I)/2 ~ (CO_IS_incon_C + CO_IS_incon_I)/2, # main effect Block Listwide control
  CO_IS_Block = (CO_IS_incon_I - CO_IS_incon_C)/2 ~ (CO_IS_con_I - CO_IS_con_C)/2, # interaction listwide effect
  PD_Congruency = (PD_IS_incon_C + PD_IS_con_C)/2 ~ (PD_IS_con_I + PD_IS_incon_I)/2, # main effect Congruence healthy controls
  PD_Itemspecific = (PD_IS_con_C + PD_IS_con_I)/2 ~ (PD_IS_incon_C + PD_IS_incon_I)/2, # main effect Block Listwide control
  PD_IS_Block = (PD_IS_incon_I - PD_IS_incon_C)/2 ~ (PD_IS_con_I - PD_IS_con_C)/2, # interaction listwide effect
  levels = c("CO_IS_incon_I", "CO_IS_con_I", "CO_IS_con_C", "CO_IS_incon_C", 
             "PD_IS_con_C", "PD_IS_con_I", "PD_IS_incon_I", "PD_IS_incon_C")
)
StroopCon_IS

# Now we create an extra column for the contrasts for each group separately
Data <- Data %>%
  mutate(Contrast_F = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "CO"  ~ "CO_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_LW_incon_I",
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_incon_I",
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "PD"  ~ "PD_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_LW_incon_I",
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_incon_I",
  )) %>%
  mutate(Contrast_F  = as_factor(Contrast_F))

# assign the generated contrast matrix to our new factor
contrasts(Data$Contrast_F) <- contr.hypothesis(StroopCon)

contrasts(Data$Contrast_F)    

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
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_incon_I",
  )) %>%
  mutate(Contrast_IS  = as_factor(Contrast_IS))

# assign the generated contrast matrix to the Item Specific Factor
contrasts(Data$Contrast_IS) <- contr.hypothesis(StroopCon_IS)

contrasts(Data$Contrast_IS)    


# Filter only young participants, remove practice trials and only keep correct trials
Data_y_inducer_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group)) 

contrasts(Data_y_inducer_LWPC$Contrast_LW) 

Data_y_diagnostic_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_inducer_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))


contrasts(Data_y_inducer_LWPC$Contrast_IS) 

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

# now that the data is prepared let us create our priors
# write down priors - intercept with a mean of 6 and a standard deviation of one seems reasonable exp(5) =150ms exp(7 = 1096ms), values that seem plausible
# for each of the effects we are careful and will assume a normal distribution centered around 0, with a standard deviation of 0.5. I derived at these values by means of prior predictive simulations

prior_weakly_informed_All <- c(
  prior(student_t(3, 6.5, 0.5), class = Intercept),
  prior(cauchy(0,0.5), class = sigma), 
  prior(cauchy(0, 5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0, 5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = b, coef = Contrast_FCO_Listwide),
  prior(normal(0,0.5), class = b, coef = Contrast_FCO_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_FCO_LW_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_FCO_IS_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_FPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.5), class = b, coef = Contrast_FPD_Listwide),
  prior(normal(0,0.5), class = b, coef = Contrast_FPD_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_FPD_LW_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_FPD_IS_Block),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Item)
)

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(student_t(3, 6.5, 0.5), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0,0.5), class = sigma), 
  prior(cauchy(0, 5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0, 5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.5), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0,0.5), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0,0.5), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0, 0.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 0.5), class = sd, coef = Intercept, group = Item)
)

# Prior informed weakly Item Specific
prior_weakly_informed_IS <- c(
  prior(student_t(3, 6.5, 0.5), class = Intercept),
  prior(cauchy(0,0.5), class = sigma), 
  prior(cauchy(0, 5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0, 5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Item)
)

# brmsformula object All
m1_all <- bf(RT ~ 1  + Contrast_F + (1|Subject) + (1|Item),
             ndt ~  0 + Exp_group) 

# brmsformula object List Wide
m1_LW <- bf(RT ~ 1  + Contrast_LW + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group) 

# brmsformula object Item Specific
m1_IS <- bf(RT ~ 1  + Contrast_IS + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group) 

#### Fit Inducer Models ####

# path where we save intermediate steps
path <- create_folder(folder_name  = "chkpt_folder_fit_inducer_model")

# we should consider varying non-decision times between the groups
fit_inducer_model_LW <- brm(formula = m1_LW,
                            family = shifted_lognormal(),
                            data = Data_y_inducer_LWPC,
                            prior = prior_weakly_informed_LW,
                            warmup = 2000,
                            iter = 10000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 423,
                            control = list(adapt_delta = 0.95),
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_inducer_model_LW, file = "fit_inducer_model_LW.rda")

fit_inducer_model_IS <- brm(formula = m1_IS,
                            family = shifted_lognormal(),
                            data = Data_y_inducer_ISPC,
                            prior = prior_weakly_informed_IS,
                            warmup = 2000,
                            iter = 10000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 412,
                            control = list(adapt_delta = 0.95),
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_inducer_model_IS, file = "fit_inducer_model_IS.rda")



#### Fit the Diagnostics ####

# Prior informed weakly List wide Diagnostic Items
prior_weakly_informed_LW_Diagnostics <- c(
  prior(normal(6.14, 0.13), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0.26,0.1), class = sigma), 
  prior(normal(5.38, 0.1), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(5.28, 0.1), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(-0.15,0.02), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0.01,0.02), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0,0.06), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(-0.12,0.02), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0.01,0.02), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0,0.04), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0.28, 0.06), class = sd, coef = Intercept, group = Subject),
  prior(normal(0.13, 0.12), class = sd, coef = Intercept, group = Item)
)

# Prior informed weakly Item Specific Diagnostic Items
prior_weakly_informed_IS_Diagnostics <- c(
  prior(normal(6.03, 0.13), class = Intercept),
  prior(normal(0.29,0.1), class = sigma), 
  prior(normal(5.56,0.1), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(5.45,0.1), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(-0.17,0.02), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0.04,0.01), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0,0.06), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(-0.12,0.02), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.02), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0,0.05), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0.31,0.06), class = sd, coef = Intercept, group = Subject),
  prior(normal(0.13,0.12), class = sd, coef = Intercept, group = Item)
)


# we should consider varying non-decision times between the groups
fit_inducer_model_LW_Diagnostics <- brm(formula = m1_LW,
                            family = shifted_lognormal(),
                            data = Data_y_diagnostic_LWPC,
                            prior = prior_weakly_informed_LW_Diagnostics,
                            warmup = 2000,
                            iter = 10000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 4233,
                            control = list(adapt_delta = 0.95),
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_inducer_model_LW_Diagnostics, file = "fit_inducer_model_LW_Diagnostics.rda")

fit_inducer_model_IS_Diagnostics <- brm(formula = m1_IS,
                            family = shifted_lognormal(),
                            data = Data_y_diagnostic_ISPC,
                            prior = prior_weakly_informed_IS_Diagnostics,
                            warmup = 2000,
                            iter = 10000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 4112,
                            control = list(adapt_delta = 0.95),
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_inducer_model_IS_Diagnostics, file = "fit_inducer_model_IS_Diagnostics.rda")


#### END Fit Diagnostics

# NEXT EXTRACT pARAMETER ESTIMATES OF THE INDUCER MODEL FOR THE NEXT PRIOR
load("fit_inducer_model_LW.rda") # Load data

# get variable names
get_variables(fit_inducer_model_LW)

post <- as_draws_array(fit_inducer_model_LW) # get posterior parameter samples

post$
quantile(post$b_intercept, probs = c(.025, .975))

prior_informed_LW <- prior(
  prior(student_t(3, 6.5, 0.5), class = Intercept),
  prior(cauchy(0,0.5), class = sigma), 
  prior(cauchy(0, 5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(cauchy(0, 5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Subject),
  prior(normal(0,0.5), class = sd, coef = Intercept, group = Item)
)
# # load model again
# load("fit_inducer_model.rda")
# # Check fit
# fit_inducer_model_LW 
# pp_check(fit_inducer_model_LW, ndraws = 100) 
# 
# # Get estimate of conditional effects by model for comparison
# con_eff <- conditional_effects(fit_inducer_model_LW)
# a <- con_eff$Contrast_F
# a$estimate__
# 
# 
# # doesn't work - check r installation
# #pp_check(fit_inducer_model, type = "stat", stat = "min")
# #pp_check(fit_inducer_model, type = "stat", stat = "max")
# #pp_check(fit_inducer_model, type = "stat", stat = "median")
# 
# # check effects
# mcmc_dens(fit_inducer_model_LW, pars = variables(fit_inducer_model_LW)[1:15])
# 
# # posterior samples aus dem Model entnehmen
# post <- posterior_samples(fit_inducer_model)
# 
# # Haupteffekt Kongruenz ausrechnen (ndt Parameter wird nicht exponentiert und rechnet)
# Congruency_CO_ms <- exp(post$b_Intercept - 0.5* post$b_Contrast_FCO_Congruency) - exp(post$b_Intercept + 0.5* post$b_Contrast_FCO_Congruency)
# c(mean = mean(Congruency_CO_ms), quantile(Congruency_CO_ms, probs = c(.025, .975)))
# 
# Congruency_PD_ms <- exp(post$b_Intercept - 0.5* post$b_Contrast_FPD_Congruency) - exp(post$b_Intercept + 0.5* post$b_Contrast_FPD_Congruency)
# c(mean = mean(Congruency_PD_ms), quantile(Congruency_PD_ms, probs = c(.025, .975)))

