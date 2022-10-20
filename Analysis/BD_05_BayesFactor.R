#####
# Author: Julius Kricheldorff
# Analysis Behavioral Data Errors - Load and Calculate Models for BF Analysis

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
library(xtable)


#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32936)



Bridge_comp <- function(full_model, deficient_model){
  #calculate marginal logliklihood using bridge sampling
  margLogLik_full <- bridge_sampler(full_model, silent = TRUE)
  margLogLik_null <- bridge_sampler(deficient_model, silent = TRUE)
  # next calculate the BF
  BF <- bayes_factor(margLogLik_full, margLogLik_null)
  return(BF$bf)
}

load_part_mod <- function(loc, model_t, item_type, effect, parameter){
  string <- loc
  if (model_t == "RT"){
    string <- file.path(string, "fit_shifted")
  } else if (model_t == "Acc"){
    string <- file.path(string, "fit")
  }
  string <- paste(string, item_type , effect, "min", parameter,sep = "_")
  fit <- load(paste(string, ".rda", sep =""))
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}

load_full_mod <- function(loc, model_t, item_type, effect){
  string <- loc
  #load the RT model - i messed up in the naming a bit here
  if (model_t == "RT"){ 
    string <- file.path(string, "fit_inducer_model")
    string <- paste(string, effect, sep = "_")
    if (item_type == "diagnostic"){
      string <- paste(string, "Diagnostics", sep = "_")
    }
    fit <- load(paste(string, ".rda", sep =""))
  } else if (model_t == "Acc"){
    string <- file.path(string, "fit_logisticMod")
    string <- paste(string, item_type , effect,sep = "_")
    fit <- load(paste(string, ".rda", sep =""))
  }
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}

BF_calc <- function(fullmod_loc, part_mod_loc){
  # tibble to save the data
  Model_BF <- tibble(
    Model = character(),
    Item_type = character(),
    Effect = character(),
    Group = character(),
    Parameter = character(),
    BF = numeric(),
  )
  mods = c("RT", "Acc")
  item_type = c("inducer","diagnostic")
  effect = c("LW", "IS")
  for (model_type in mods){
    for (itm in item_type){
      for(eff in effect){
        if (eff == "LW"){
          parameter <- c(
            "CO_Congruency",
            "CO_Listwide",
            "CO_LW_Block",
            "PD_Congruency",
            "PD_Listwide",
            "PD_LW_Block"
          )
        } else if (eff == "IS"){
          parameter <- c(
            "CO_Congruency",
            "CO_Itemspecific",
            "CO_IS_Block",
            "PD_Congruency",
            "PD_Itemspecific",
            "PD_IS_Block"
          )
        }
        # load the full model
        full_mod <- load_full_mod(loc = fullmod_loc, 
                                  model_t = model_type,
                                  item_type = itm, 
                                  effect = eff)
        for (par in parameter){
          # load the partial model
          def_mod <- load_part_mod(loc = part_mod_loc, 
                                   model_t = model_type, 
                                   item_type = itm, 
                                   effect = eff, 
                                   parameter = par)
          # claculate the BF
          BF_inc <- Bridge_comp(full_model = full_mod,
                                deficient_model = def_mod)
          #save results as a tibble
          temp <- tibble(
            Model = model_type,
            Item_type = itm,
            Effect = eff,
            Group = substr(par, start = 1, stop = 2),
            Parameter = substr(par, start = 4, stop = nchar(par)),
            BF = BF_inc,
          )
          # join tibble as row 
          Model_BF <- bind_rows(Model_BF, temp)
        }
      }
    }
  }
  return(Model_BF)
}

# calculate effect estimates in ms
conditional_effect_calc_shift <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.25*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.25*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.25*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.25*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.25*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.25*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.25*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.25*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
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
      Control_CO = Stroop_MC_CO - Stroop_MI_CO)
  } else if (effect == "IS"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.25*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.25*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.25*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.25*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.25*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.25*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.25*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.25*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
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
      Control_CO = Stroop_MC_CO - Stroop_MI_CO)
  }
  return(model_effects)
}

conditional_effect_calc_acc <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.25*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.25*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.25*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.25*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.25*m_post$b_Contrast_LWCO_LW_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.25*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.25*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.25*m_post$b_Contrast_LWCO_LW_Block))*100
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
      Control_CO = Stroop_MC_CO - Stroop_MI_CO)
  } else if (effect == "IS"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.25*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.25*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.25*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.25*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.25*m_post$b_Contrast_ISCO_IS_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.25*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.25*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.25*m_post$b_Contrast_ISCO_IS_Block))*100
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
      Control_CO = Stroop_MC_CO - Stroop_MI_CO)
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
      mean = mean(var_t),
      lower95 = quantile(var_t, probs = 0.025),
      upper95 = quantile(var_t, probs = 0.975)
    )
    # join new row with summary table
    model_summary <- bind_rows(model_summary, temp)
  }
  return(model_summary)
}

mean_par <- function(model, parameter, effect) # keep that for later not in this table 
  #get mean pars from model by first getting the summary statistics
  modsum <- summary(model)
  modsum$fixed["Intercept","Estimate"]
  
### First we get all the variables we need to call our readily calculated models


# contrast names seperately

Partial_models_saveloc <- "E:/AdaptiveControl/Data/BehaviorResults"
Full_models_saveloc <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data"

# calculate BFs
BF_Results <- BF_calc( 
        fullmod_loc = Full_models_saveloc, 
        part_mod_loc = Partial_models_saveloc)

#save data
write.table(BF_Results , file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/BF_Results.csv")

# load data again
BF_Results <- read.csv(file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/BF_Results.csv", header = TRUE, sep = "")
#### Now that we have the table with our BFs, let us load the parameter estimates

# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(mean = NA , lower95 = NA, upper95 = NA) 

# Load and save the data

mods = c("RT", "Acc")
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
          par = "Congruency"
        } 
        if (grepl("Block" ,vars)){
            par = "Listwide" 
        } 
        if (grepl("Control" ,vars)) {
            par = "LW_Block"
        }
        #NOTE: by accident when i created the models for the BF calculation, also the models
        # get a single row
        row_var <- sum_t %>% filter(parameter == vars)
        # Write values
        Full_Model_Info$mean[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$mean
        Full_Model_Info$lower95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$lower95
        Full_Model_Info$upper95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$upper95
      }
    }
  }
}

#NOTE: by accident when i created the models for the BF calculation, also the models
#Note: When i created the models for the BF calculation, by accident, I named them the same 
# for the item specific effect as the list wise effect. So the itemspecific model actually is named listwise (variable names are correct though)
# same goes for the interaction. Thus in the above code only LW_Block and Listwise are used. Does not affect the results though.
#For better readability I correct that below
Full_Model_Info_fin <- Full_Model_Info %>%
  mutate(Parameter = case_when(
    Effect == "IS" & Parameter == "Listwide" ~ "Itemspecific",
    Effect == "IS" & Parameter == "LW_Block" ~ "IS_Block",
    TRUE ~ Parameter))

# Export as latex table
print(xtable(Full_Model_Info_fin, type = "latex"), file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/Bheav_Summary.tex")

# Filter only the RT Analysis Results

Full_Model_Info_fin_RT <- Full_Model_Info_fin %>%
  filter(Model == "RT") %>%
  filter(Effect == "LW") %>%
  mutate(mean = round(mean,1)) %>%
  mutate(lower95 = round(lower95,1)) %>%
  mutate(upper95 = round(upper95,1)) %>%
  mutate(BF = round(BF,2))
  
print(xtable(Full_Model_Info_fin_RT, type = "latex"), file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/Bheav_Summary_RT.tex")



# test

load("C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/fit_inducer_model_LW_Diagnostics.rda")
fit_inducer_model_LW_Diagnostics

load("E:/AdaptiveControl/Data/BehaviorResults/fit_diagnostic_LW_min_CO_LW_Block.rda")

val <- Bridge_comp(fit_inducer_model_LW_Diagnostics, fit_model)