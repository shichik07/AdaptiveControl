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
  filter(Subject != "sub-CO_0" &  Subject != "sub-CO_3" & Subject != "sub-CO_7" &  Subject != "sub-CO_9") %>%# exclude because they used both hands
  mutate(Error = case_when( # we define accuracy in terms of errors, easier to interpret
    Correct_A == 1 ~ 0,
    Correct_A == 0 ~1
  )) %>%
  group_by(Subject,Block) %>% # so that we do not exclude trials that follow across blocks or subjects
  mutate(Post_error = lag(Error))%>%
  ungroup() %>%
  filter( RT > 200,
          Post_error ==1)# Remove data with small RT, and Post error trials



# Contrasts only for the list-wide effect only
StroopCon_LW <- hypr(
  CO_Congruency = (CO_LW_incon_C + CO_LW_con_C)/2 ~ (CO_LW_con_I + CO_LW_incon_I)/2, # main effect Congruence healthy controls
  CO_Listwide = (CO_LW_con_C + CO_LW_con_I)/2 ~ (CO_LW_incon_C + CO_LW_incon_I)/2, # main effect Block Listwide control
  CO_LW_Block = (CO_LW_incon_I - CO_LW_incon_C)/2 ~ (CO_LW_con_I - CO_LW_con_C)/2, # interaction listwide effect
  PD_Congruency = (PD_LW_incon_C + PD_LW_con_C)/2 ~ (PD_LW_con_I + PD_LW_incon_I)/2, # main effect Congruence healthy controls
  PD_Listwide = (PD_LW_con_C + PD_LW_con_I)/2 ~ (PD_LW_incon_C + PD_LW_incon_I)/2, # main effect Block Listwide control
  PD_LW_Block = (PD_LW_incon_I - PD_LW_incon_C)/2 ~ (PD_LW_con_I - PD_LW_con_C)/2, # interaction listwide effect
  levels = c("CO_LW_con_C", "CO_LW_con_I", "CO_LW_incon_I", "CO_LW_incon_C", 
             "PD_LW_incon_C", "PD_LW_incon_I", "PD_LW_con_I", "PD_LW_con_C")
  # levels = c("CO_LW_con_C", "CO_LW_con_I", "CO_LW_incon_I", "CO_LW_incon_C", 
  #            "PD_LW_incon_I", "PD_LW_incon_C", "PD_LW_con_C", "PD_LW_con_I")
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

#### Fit Inducer Models ####

# path where we save intermediate steps
#path <- create_folder(folder_name  = "chkpt_folder_fit_inducer_model")

# we should consider varying non-decision times between the groups
fit_inducer_model_LW <- brm(formula = m1_LW,
                            family = shifted_lognormal(),
                            data = Data_y_inducer_LWPC,
                            prior = prior_weakly_informed_LW,
                            warmup = 2000,
                            iter = 12000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 423,
                            control = list(adapt_delta = 0.9),
                            save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)
save(fit_inducer_model_LW, file = "Post_error_fit_inducer_model_LW.rda")


fit_inducer_model_IS <- brm(formula = m1_IS,
                            family = shifted_lognormal(),
                            data = Data_y_inducer_ISPC,
                            prior = prior_weakly_informed_IS,
                            warmup = 2000,
                            iter = 12000,# 20000 is the limit necessary for bridge sampling
                            cores = 4, seed = 412,
                            control = list(adapt_delta = 0.9),
                            save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            chains =4
)

save(fit_inducer_model_IS, file = "Post_error_fit_inducer_model_IS.rda")



#### Fit the Diagnostics ####
load("Post_error_fit_inducer_model_LW.rda")
load("Post_error_fit_inducer_model_IS.rda")
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
  prior(normal(Fixed_Intercept_mean, Fixed_Intercept_sd), class = Intercept, lb = 0), 
  prior(normal(Fixed_Contrast_LWCO_Congruency_mean, Fixed_Contrast_LWCO_Congruency_sd), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_LWCO_Listwide_mean, Fixed_Contrast_LWCO_Listwide_sd), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(Fixed_Contrast_LWCO_LW_Block_mean, Fixed_Contrast_LWCO_LW_Block_sd), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(Fixed_Contrast_LWPD_Congruency_mean, Fixed_Contrast_LWPD_Congruency_sd), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_LWPD_Listwide_mean, Fixed_Contrast_LWPD_Listwide_sd), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(Fixed_Contrast_LWPD_LW_Block_mean, Fixed_Contrast_LWPD_LW_Block_sd), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(Item_mean, Item_sd), class = sd, coef = Intercept, group = Item),
  prior(normal(Subject_mean, Subject_sd), class = sd, coef = Intercept, group = Subject),
  prior(normal(Sigma_mean, Sigma_sd), class = sigma, lb = 0),
  prior(normal(ndt_Exp_groupCO_mean, ndt_Exp_groupCO_sd), class = b, coef = Exp_groupCO, dpar = ndt), 
  prior(normal(ndt_Exp_groupPD_mean, ndt_Exp_groupPD_sd), class = b, coef = Exp_groupPD, dpar = ndt)
)

# now put in our variables into the prior for IS
prior_informed_IS <- c(
  prior(normal(Fixed_Intercept_mean, Fixed_Intercept_sd), class = Intercept, lb = 0), 
  prior(normal(Fixed_Contrast_ISCO_Congruency_mean, Fixed_Contrast_ISCO_Congruency_sd), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_ISCO_Itemspecific_mean, Fixed_Contrast_ISCO_Itemspecific_sd), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(Fixed_Contrast_ISCO_IS_Block_mean, Fixed_Contrast_ISCO_IS_Block_sd), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(Fixed_Contrast_ISPD_Congruency_mean, Fixed_Contrast_ISPD_Congruency_sd), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(Fixed_Contrast_ISPD_Itemspecific_mean, Fixed_Contrast_ISPD_Itemspecific_sd), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(Fixed_Contrast_ISPD_IS_Block_mean, Fixed_Contrast_ISPD_IS_Block_sd), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(Item_mean, Item_sd), class = sd, coef = Intercept, group = Item),
  prior(normal(Subject_mean, Subject_sd), class = sd, coef = Intercept, group = Subject),
  prior(normal(Sigma_mean, Sigma_sd), class = sigma, lb = 0),
  prior(normal(ndt_Exp_groupCO_mean, ndt_Exp_groupCO_sd), class = b, coef = Exp_groupCO, dpar = ndt), 
  prior(normal(ndt_Exp_groupPD_mean, ndt_Exp_groupPD_sd), class = b, coef = Exp_groupPD, dpar = ndt)
)


# we should consider varying non-decision times between the groups
fit_inducer_model_LW_Diagnostics <- brm(formula = m1_LW,
                                        family = shifted_lognormal(),
                                        data = Data_y_diagnostic_LWPC,
                                        prior = prior_informed_LW,
                                        warmup = 2000,
                                        iter = 12000,# 20000 is the limit necessary for bridge sampling
                                        cores = 4, seed = 4233,
                                        control = list(adapt_delta = 0.9),
                                        stanvars = stan_vars_LW,
                                        save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                        chains =4
)

save(fit_inducer_model_LW_Diagnostics, file = "Post_error_fit_inducer_model_LW_Diagnostics.rda")

fit_inducer_model_IS_Diagnostics <- brm(formula = m1_IS,
                                        family = shifted_lognormal(),
                                        data = Data_y_diagnostic_ISPC,
                                        prior = prior_informed_IS,
                                        warmup = 2000,
                                        iter = 12000,# 20000 is the limit necessary for bridge sampling
                                        cores = 4, seed = 4112,
                                        control = list(adapt_delta = 0.9),
                                        stanvars = stan_vars_IS,
                                        save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                        chains =4
)

save(fit_inducer_model_IS_Diagnostics, file = "Post_error_fit_inducer_model_IS_Diagnostics.rda")


#### END Fit Diagnostics

# Now lets load the parameter estimates
load_full_mod <- function(loc, model_t, item_type, effect){
  string <- loc
  #load the RT model - i messed up in the naming a bit here
  if (model_t == "RT"){ 
    string <- file.path(string, "Post_error_fit_inducer_model")
    string <- paste(string, effect, sep = "_")
    if (item_type == "diagnostic"){
      string <- paste(string, "Diagnostics", sep = "_")
    }
    fit <- load(paste(string, ".rda", sep =""))
  } else if (model_t == "Acc"){
    string <- file.path(string, "Post_error_fit_logisticMod")
    string <- paste(string, item_type , effect,sep = "_")
    fit <- load(paste(string, ".rda", sep =""))
  }
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}


# calculate effect estimates in ms
conditional_effect_calc_shift <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)

    model_effects <- tibble(
      MC_C_PD = MC_C_PD,
      MC_I_PD = MC_I_PD,
      MI_C_PD = MI_C_PD,
      MI_I_PD = MI_I_PD,
      MC_C_CO = MC_C_CO,
      MC_I_CO = MC_I_CO,
      MI_C_CO = MI_C_CO,
      MI_I_CO = MI_I_CO,
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      Dif_CO = Control_CO - Control_PD,
      Conflictimprov_PD = MI_I_PD - MC_I_PD,
      Conflictimprov_CO = MI_I_CO - MC_I_CO,
      Improv_dif = Conflictimprov_CO -Conflictimprov_PD)
  } else if (effect == "IS"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    model_effects <- tibble(
      MC_C_PD = MC_C_PD,
      MC_I_PD = MC_I_PD,
      MI_C_PD = MI_C_PD,
      MI_I_PD = MI_I_PD,
      MC_C_CO = MC_C_CO,
      MC_I_CO = MC_I_CO,
      MI_C_CO = MI_C_CO,
      MI_I_CO = MI_I_CO,
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      Dif_CO = Control_CO - Control_PD,
      Conflictimprov_PD = MI_I_PD - MC_I_PD,
      Conflictimprov_CO = MI_I_CO - MC_I_CO,
      Improv_dif = Conflictimprov_CO -Conflictimprov_PD)
  }
  return(model_effects)
}

conditional_effect_calc_acc <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      Dif_CO = Control_CO - Control_PD,
      Conflictimprov_PD = MI_I_PD - MC_I_PD,
      Conflictimprov_CO = MI_I_CO - MC_I_CO,
      Improv_dif = Conflictimprov_CO -Conflictimprov_PD)
  } else if (effect == "IS"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      Dif_CO = Control_CO - Control_PD,
      Conflictimprov_PD = MI_I_PD - MC_I_PD,
      Conflictimprov_CO = MI_I_CO - MC_I_CO,
      Improv_dif = Conflictimprov_CO -Conflictimprov_PD)
  }
  return(model_effects)
}

# function calculate summary statistics from posterior distributions
post_sum_calc <- function(summary_table){
  col_n <- colnames(summary_table)
  model_summary <- tibble(
    parameter = character(),
    mean = numeric(),
    lower95 = numeric(),
    upper95 = numeric()
  )
  for (var_n in col_n){
    var_t <- eval(parse(text = paste('summary_table$', var_n)))
    temp <- tibble(
      parameter = var_n,
      mean = round(mean(var_t), 1),
      lower95 = round(quantile(var_t, probs = 0.025), 1),
      upper95 = round(quantile(var_t, probs = 0.975), 1)
    )
    # join new row with summary table
    model_summary <- bind_rows(model_summary, temp)
  }
  return(model_summary)
}

### First we get all the variables we need to call our readily calculated models


# contrast names seperately

Partial_models_saveloc <- "E:/AdaptiveControl/Data/BehaviorResults"
Full_models_saveloc <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data"

#### Now that we have the table with our BFs, let us load the parameter estimates

#load BF data
BF_Results<- read.table(paste0(Full_models_saveloc, '/BF_results.csv'))


# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(mean = NA , lower95 = NA, upper95 = NA) 

# Load and save the data

mods = c("RT")
item_type = c("inducer","diagnostic")
effect = c("LW", "IS")

for (md in mods) {
  for (eff in effect){
    for (itm in item_type){
      #load model
      model <- load_full_mod(loc = Full_models_saveloc,
                             model_t = md,
                             item_type = itm,
                             effect = eff)
      
      # calculate posterior conditional effect samples
      if (md == "RT"){
        eff_post <- conditional_effect_calc_shift(effect = eff,
                                                  model = model)
      } else if (md == "Acc"){
        eff_post <- conditional_effect_calc_acc(effect = eff,
                                                model = model)
      }
      #translate into summary stats
      sum_t <- post_sum_calc(eff_post)
      
      #integrate summary stats into BF table
      
      temp_t <- sum_t %>% filter(!grepl("Stroop", parameter)) # filter the stroop effects which we are not interested in to show
      for (vars in temp_t$parameter){
        group_in <- substr(vars, start = nchar(vars)-1, stop = nchar(vars))
        if (grepl("Congruency" ,vars)){
          parn = "Congruency"
        } 
        else if (grepl("Block" ,vars)){
          if (eff == "LW"){
            parn = "Listwide"
          } else{
            parn = "Itemspecific"
          }
          
        } 
        else if (grepl("Control" ,vars)) {
          if (eff == "LW"){
            parn = "LW_Block"
          } else{
            parn = "IS_Block"
          }
        }
        else{
          parn = "something"
        }
        #NOTE: by accident when i created the models for the BF calculation, also the models
        # get a single row
        row_var <- sum_t %>% filter(parameter == vars)
        # Write values
        Full_Model_Info$mean[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == parn] <- row_var$mean
        Full_Model_Info$lower95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == parn] <- row_var$lower95
        Full_Model_Info$upper95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == parn] <- row_var$upper95
      }
    }
  }
}


# fit values
conditional_effects(model)
fitted_values <- fitted(model)
head(fitted)

##### Inducer
CO_LW_con_C <-  exp(6 + -0.18*0.5 + 0.5*0.01 + 0.5*-0.05)
CO_LW_con_I <-  exp(6 + -0.18*-0.5 + 0.5*0.01 + -0.5*-0.05)
CO_LW_incon_C <-  exp(6 + -0.18*0.5 - 0.5*0.01 + -0.5*-0.05)
CO_LW_incon_I <-  exp(6 + -0.18*-0.5 - 0.5*0.01 + 0.5*-0.05)

AdapConCO <- (CO_LW_con_C - CO_LW_con_I) - (CO_LW_incon_C - CO_LW_incon_I)

PD_LW_con_C <-  exp(6 + -0.13*0.5 + 0.5*0.01 + 0.5*-0.04)
PD_LW_con_I <-  exp(6 + -0.13*-0.5 + 0.5*0.01 + -0.5*-0.04)
PD_LW_incon_C <-  exp(6 + -0.13*0.5 - 0.5*0.01 + -0.5*-0.04)
PD_LW_incon_I <-  exp(6 + -0.13*-0.5 - 0.5*0.01 + 0.5*-0.04)

AdapConPD <- (PD_LW_con_C - PD_LW_con_I) - (PD_LW_incon_C - PD_LW_incon_I)

##### Diagnostics
CO_LW_con_C_dia <-  exp(6 + -0.18*0.5 + 0.5*0.02 + 0.5*-0.04)
CO_LW_con_I_dia <-  exp(6 + -0.18*-0.5 + 0.5*0.02 + -0.5*-0.04)
CO_LW_incon_C_dia <-  exp(6 + -0.18*0.5 - 0.5*0.02 + -0.5*-0.04)
CO_LW_incon_I_dia <-  exp(6 + -0.18*-0.5 - 0.5*0.02 + 0.5*-0.04)

AdapConCO_dia <- (CO_LW_con_C_dia - CO_LW_con_I_dia) - (CO_LW_incon_C_dia - CO_LW_incon_I_dia)

PD_LW_con_C_dia <-  exp(6 + -0.13*0.5 + 0.5*-0.01 + 0.5*-0.02)
PD_LW_con_I_dia <-  exp(6 + -0.13*-0.5 + 0.5*-0.01 + -0.5*-0.02)
PD_LW_incon_C_dia <-  exp(6 + -0.13*0.5 - 0.5*-0.01 + -0.5*-0.02)
PD_LW_incon_I_dia <-  exp(6 + -0.13*-0.5 - 0.5*-0.01 + 0.5*-0.02)

AdapConPD_dia <- (PD_LW_con_C_dia - PD_LW_con_I_dia) - (PD_LW_incon_C_dia - PD_LW_incon_I_dia)