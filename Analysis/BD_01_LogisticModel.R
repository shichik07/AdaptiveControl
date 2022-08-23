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
  dplyr::filter(Trl_type != "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group)) 

contrasts(Data_y_inducer_LWPC$Contrast_LW) 

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


contrasts(Data_y_inducer_LWPC$Contrast_IS) 

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  mutate(Exp_group  = as_factor(Exp_group))

# For the priors we assume a base error rate for the intercept about 6% and for the conditions we will say 3% change is reasonable to be reflected as one std
# Because that does not really work in log odds, we use a normal with mean -1.3 and scale 1.5 for the intercept.
# This is weakly informative in the sense that most probability mass is in the lower percentage points (0-25%).
# We know participants never make more than 10% errors in these kind of tasks. 
# For the effect estimates we set all estimates to 0 (could be increase, could be decrease, with a variance of 1.5)


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

#### Fit Inducer Models ####


# we should consider varying non-decision times between the groups
fit_logisticMod_inducer_LW <- brm(formula = m1_LW,
                            family = bernoulli(link = logit),
                            data = Data_y_inducer_LWPC,
                            prior = prior_weakly_informed_LW,
                            warmup = 2000,
                            iter = 22000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 423,
                            control = list(adapt_delta = 0.95),
                            save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_logisticMod_inducer_LW, file = "fit_logisticMod_inducer_LW.rda")
load(file = "fit_logisticMod_inducer_LW.rda")

fit_logisticMod_inducer_IS <- brm(formula = m1_IS,
                            family = bernoulli(link = logit),
                            data = Data_y_inducer_ISPC,
                            prior = prior_weakly_informed_IS,
                            warmup = 2000,
                            iter = 22000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 412,
                            control = list(adapt_delta = 0.95),
                            save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_logisticMod_inducer_IS, file = "fit_logisticMod_inducer_IS.rda")
load(file = "fit_logisticMod_inducer_IS.rda")


### Now that the inducer models are run, we use the posteriors as informed priors for the diagnostic models

# Prior informed List wide
prior_informed_LW <- c(
  prior(normal(-4.4, 0.4), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(-1.67, 0.45), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0.41, 0.46), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(-0.64, 0.46), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(-1.48, 0.38), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0.39, 0.38), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(-0.14, 0.38), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0.72, 0.23), class = sd, coef = Intercept, group = Subject),
  prior(normal(0.46, 0.45), class = sd, coef = Intercept, group = Item)
)

# Prior informed Item Specific
prior_informed_IS <- c(
  prior(normal(-4.51, 0.6), class = Intercept), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(-1.76, 0.53), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(-0.1, 0.54), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(-0.52, 0.5), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(-0.88, 0.25), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0.2, 0.26), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(-0.59, 0.25), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(1.25 , 0.22), class = sd, coef = Intercept, group = Subject),
  prior(normal(0.58 , 0.65), class = sd, coef = Intercept, group = Item)
)


# we should consider varying non-decision times between the groups
fit_logisticMod_diagnostic_LW <- brm(formula = m1_LW,
                                  family = bernoulli(link = logit),
                                  data = Data_y_diagnostic_LWPC,
                                  prior = prior_informed_LW,
                                  warmup = 2000,
                                  iter = 22000,# 20000 is the limit necessary for bridge sampling
                                  cores = 4, seed = 423,
                                  control = list(adapt_delta = 0.95),
                                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                  chains =4
)

save(fit_logisticMod_diagnostic_LW, file = "fit_logisticMod_diagnostic_LW.rda")
#load(file = "fit_logisticMod_diagnostic_LW.rda")

fit_logisticMod_diagnostic_IS <- brm(formula = m1_IS,
                                  family = bernoulli(link = logit),
                                  data = Data_y_diagnostic_ISPC,
                                  prior = prior_informed_IS,
                                  warmup = 2000,
                                  iter = 22000,# 20000 is the limit necessary for bridge sampling
                                  cores = 4, seed = 412,
                                  control = list(adapt_delta = 0.95),
                                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                  chains =4
)

save(fit_logisticMod_diagnostic_IS, file = "fit_logisticMod_diagnostic_IS.rda")
load(file = "fit_logisticMod_diagnostic_IS.rda")









function(proportion){
  
}

samples_logodds <- tibble(alpha =rnorm(100000, -2, 2))
samples_prob <- tibble(p = plogis(rnorm(100000, -3.8, 0.4)))
ggplot(samples_logodds, aes(alpha)) +
  geom_density()
ggplot(samples_prob, aes(p)) +
  geom_density()



a <- get_prior(data = Data_y_inducer_LWPC, m1_LW, family = bernoulli(link = logit))  
m1_LW_1 <- bf(RT ~ 1  + Contrast_LW + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group) 
b <- get_prior(data = Data_y_inducer_LWPC, m1_LW_1, family = shifted_lognormal())  


