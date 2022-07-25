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

# Set a seed for sake of reproducibility
set.seed(32946)

load("fit_inducer_model_LW.rda") # Load data Listwise Inducer
load("fit_inducer_model_LW_Diagnostics.rda") # Load data Listwise Diagnostics
load("fit_inducer_model_IS.rda") # Load data Item-specific Inducer
load("fit_inducer_model_IS_Diagnostics.rda") # Load data Item-specific Diagnostics


# get variable names
variables(fit_inducer_model_LW)
variables(fit_inducer_model_IS)


# posterior samples aus dem Model entnehmen
post_LW_Ind <- posterior_samples(fit_inducer_model_LW)
post_LW_Dia <- posterior_samples(fit_inducer_model_LW_Diagnostics)
post_IS_Ind <- posterior_samples(fit_inducer_model_IS)
post_IS_Dia <- posterior_samples(fit_inducer_model_IS_Diagnostics)

#### Effekte Inducer LW ####

# Haupteffekt Kongruenz ausrechnen (ndt Parameter wird nicht exponentiert und rechnet)
Congruency_ms <- exp((post_LW_Ind$b_Intercept - 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency)+ post_LW_Ind$sigma/2)  - exp((post_LW_Ind$b_Intercept + 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency) + post_LW_Ind$sigma/2)
c(mean = mean(Congruency_ms), quantile(Congruency_ms, probs = c(.025, .975)))

# Marginal Means ausrechnen PD
LW_MC_C_PD <- exp((post_LW_Ind$b_Intercept + 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency + 0.5 * post_LW_Ind$b_Contrast_LWPD_Listwide  + 0.25*post_LW_Ind$b_Contrast_LWPD_LW_Block) + post_LW_Ind$sigma/2)
LW_MC_I_PD <- exp((post_LW_Ind$b_Intercept - 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency + 0.5 * post_LW_Ind$b_Contrast_LWPD_Listwide  - 0.25*post_LW_Ind$b_Contrast_LWPD_LW_Block) + post_LW_Ind$sigma/2)
LW_MI_C_PD <- exp((post_LW_Ind$b_Intercept + 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency - 0.5 * post_LW_Ind$b_Contrast_LWPD_Listwide  - 0.25*post_LW_Ind$b_Contrast_LWPD_LW_Block) + post_LW_Ind$sigma/2)
LW_MI_I_PD <- exp((post_LW_Ind$b_Intercept - 0.5* post_LW_Ind$b_Contrast_LWPD_Congruency - 0.5 * post_LW_Ind$b_Contrast_LWPD_Listwide  + 0.25*post_LW_Ind$b_Contrast_LWPD_LW_Block) + post_LW_Ind$sigma/2)

# Mean und Interval der Marginal Means PD
c(condition = 'LW_MC_C_PD', mean = mean(LW_MC_C_PD), quantile(LW_MC_C_PD, probs = c(.025, .975)))
c(condition = 'LW_MC_I_PD', mean = mean(LW_MC_I_PD), quantile(LW_MC_I_PD, probs = c(.025, .975)))
c(condition = 'LW_MI_C_PD', mean = mean(LW_MI_C_PD), quantile(LW_MI_C_PD, probs = c(.025, .975)))
c(condition = 'LW_MI_I_PD', mean = mean(LW_MI_I_PD), quantile(LW_MI_I_PD, probs = c(.025, .975)))


# Mit den Marginal Means PD
Stroop_MC_PD_LW_Ind <- (LW_MC_I_PD - LW_MC_C_PD)
c(condition = 'Stroop_MC_PD_LW_Ind', mean = mean(Stroop_MC_PD_LW_Ind), quantile(Stroop_MC_PD_LW_Ind, probs = c(.025, .975)))
Stroop_MI_PD_LW_Ind <- (LW_MI_I_PD - LW_MI_C_PD)
c(condition = 'Stroop_MI_PD_LW_Ind', mean = mean(Stroop_MI_PD_LW_Ind), quantile(Stroop_MI_PD_LW_Ind, probs = c(.025, .975)))
LW_Block_PD_Ind <- Stroop_MC_PD_LW_Ind - Stroop_MI_PD_LW_Ind
c(condition = 'Interaction LW_Block', mean = mean(LW_Block_PD_Ind), quantile(LW_Block_PD_Ind, probs = c(.025, .975)))


# Marginal Means ausrechnen CO
LW_MC_C_CO <- exp((post_LW_Ind$b_Intercept + 0.5* post_LW_Ind$b_Contrast_LWCO_Congruency + 0.5 * post_LW_Ind$b_Contrast_LWCO_Listwide  + 0.25*post_LW_Ind$b_Contrast_LWCO_LW_Block) + post_LW_Ind$sigma/2)
LW_MC_I_CO <- exp((post_LW_Ind$b_Intercept - 0.5* post_LW_Ind$b_Contrast_LWCO_Congruency + 0.5 * post_LW_Ind$b_Contrast_LWCO_Listwide  - 0.25*post_LW_Ind$b_Contrast_LWCO_LW_Block) + post_LW_Ind$sigma/2)
LW_MI_C_CO <- exp((post_LW_Ind$b_Intercept + 0.5* post_LW_Ind$b_Contrast_LWCO_Congruency - 0.5 * post_LW_Ind$b_Contrast_LWCO_Listwide  - 0.25*post_LW_Ind$b_Contrast_LWCO_LW_Block) + post_LW_Ind$sigma/2)
LW_MI_I_CO <- exp((post_LW_Ind$b_Intercept - 0.5* post_LW_Ind$b_Contrast_LWCO_Congruency - 0.5 * post_LW_Ind$b_Contrast_LWCO_Listwide  + 0.25*post_LW_Ind$b_Contrast_LWCO_LW_Block) + post_LW_Ind$sigma/2)

# Mean und Interval der Marginal Means CO
c(condition = 'LW_MC_C_CO', mean = mean(LW_MC_C_CO), quantile(LW_MC_C_CO, probs = c(.025, .975)))
c(condition = 'LW_MC_I_CO', mean = mean(LW_MC_I_CO), quantile(LW_MC_I_CO, probs = c(.025, .975)))
c(condition = 'LW_MI_C_CO', mean = mean(LW_MI_C_CO), quantile(LW_MI_C_CO, probs = c(.025, .975)))
c(condition = 'LW_MI_I_CO', mean = mean(LW_MI_I_CO), quantile(LW_MI_I_CO, probs = c(.025, .975)))


# Mit den Marginal Means CO
Stroop_MC_CO_LW_Ind <- (LW_MC_I_CO - LW_MC_C_CO)
c(condition = 'Stroop_MC_CO_LW_Ind', mean = mean(Stroop_MC_CO_LW_Ind), quantile(Stroop_MC_CO_LW_Ind, probs = c(.025, .975)))
Stroop_MI_CO_LW_Ind<- (LW_MI_I_CO - LW_MI_C_CO)
c(condition = 'Stroop_MI_CO_LW_Ind', mean = mean(Stroop_MI_CO_LW_Ind), quantile(Stroop_MI_CO_LW_Ind, probs = c(.025, .975)))
LW_Block_CO_Ind <- Stroop_MC_CO_LW_Ind - Stroop_MI_CO_LW_Ind
c(condition = 'Interaction LW_Block', mean = mean(LW_Block_CO_Ind), quantile(LW_Block_CO_Ind, probs = c(.025, .975)))

remove('LW_MC_C_CO', 'LW_MC_I_CO', 'LW_MI_C_CO', 'LW_MI_I_CO', 'LW_MC_C_PD', 'LW_MC_I_PD', 'LW_MI_C_PD', 'LW_MI_I_PD' )

#### Effekte Diagnostics LW ####

# Haupteffekt Kongruenz ausrechnen (ndt Parameter wird nicht exponentiert und rechnet)
Congruency_ms <- exp((post_LW_Dia$b_Intercept - 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency)+ post_LW_Dia$sigma/2)  - exp((post_LW_Dia$b_Intercept + 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency) + post_LW_Dia$sigma/2)
c(mean = mean(Congruency_ms), quantile(Congruency_ms, probs = c(.025, .975)))

# Marginal Means ausrechnen PD
LW_MC_C_PD <- exp((post_LW_Dia$b_Intercept + 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency + 0.5 * post_LW_Dia$b_Contrast_LWPD_Listwide  + 0.25*post_LW_Dia$b_Contrast_LWPD_LW_Block) + post_LW_Dia$sigma/2)
LW_MC_I_PD <- exp((post_LW_Dia$b_Intercept - 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency + 0.5 * post_LW_Dia$b_Contrast_LWPD_Listwide  - 0.25*post_LW_Dia$b_Contrast_LWPD_LW_Block) + post_LW_Dia$sigma/2)
LW_MI_C_PD <- exp((post_LW_Dia$b_Intercept + 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency - 0.5 * post_LW_Dia$b_Contrast_LWPD_Listwide  - 0.25*post_LW_Dia$b_Contrast_LWPD_LW_Block) + post_LW_Dia$sigma/2)
LW_MI_I_PD <- exp((post_LW_Dia$b_Intercept - 0.5* post_LW_Dia$b_Contrast_LWPD_Congruency - 0.5 * post_LW_Dia$b_Contrast_LWPD_Listwide  + 0.25*post_LW_Dia$b_Contrast_LWPD_LW_Block) + post_LW_Dia$sigma/2)

# Mean und Interval der Marginal Means PD
c(condition = 'LW_MC_C_PD', mean = mean(LW_MC_C_PD), quantile(LW_MC_C_PD, probs = c(.025, .975)))
c(condition = 'LW_MC_I_PD', mean = mean(LW_MC_I_PD), quantile(LW_MC_I_PD, probs = c(.025, .975)))
c(condition = 'LW_MI_C_PD', mean = mean(LW_MI_C_PD), quantile(LW_MI_C_PD, probs = c(.025, .975)))
c(condition = 'LW_MI_I_PD', mean = mean(LW_MI_I_PD), quantile(LW_MI_I_PD, probs = c(.025, .975)))


# Mit den Marginal Means PD
Stroop_MC_PD_LW_Dia <- (LW_MC_I_PD - LW_MC_C_PD)
c(condition = 'Stroop_MC_PD_LW_Dia', mean = mean(Stroop_MC_PD_LW_Dia), quantile(Stroop_MC_PD_LW_Dia, probs = c(.025, .975)))
Stroop_MI_PD_LW_Dia <- (LW_MI_I_PD - LW_MI_C_PD)
c(condition = 'Stroop_MI_PD_LW_Dia', mean = mean(Stroop_MI_PD_LW_Dia), quantile(Stroop_MI_PD_LW_Dia, probs = c(.025, .975)))
LW_Block_PD_Dia <- Stroop_MC_PD_LW_Dia - Stroop_MI_PD_LW_Dia
c(condition = 'Interaction LW_Block', mean = mean(LW_Block_PD_Dia), quantile(LW_Block_PD_Dia, probs = c(.025, .975)))


# Marginal Means ausrechnen CO
LW_MC_C_CO <- exp((post_LW_Dia$b_Intercept + 0.5* post_LW_Dia$b_Contrast_LWCO_Congruency + 0.5 * post_LW_Dia$b_Contrast_LWCO_Listwide  + 0.25*post_LW_Dia$b_Contrast_LWCO_LW_Block) + post_LW_Dia$sigma/2)
LW_MC_I_CO <- exp((post_LW_Dia$b_Intercept - 0.5* post_LW_Dia$b_Contrast_LWCO_Congruency + 0.5 * post_LW_Dia$b_Contrast_LWCO_Listwide  - 0.25*post_LW_Dia$b_Contrast_LWCO_LW_Block) + post_LW_Dia$sigma/2)
LW_MI_C_CO <- exp((post_LW_Dia$b_Intercept + 0.5* post_LW_Dia$b_Contrast_LWCO_Congruency - 0.5 * post_LW_Dia$b_Contrast_LWCO_Listwide  - 0.25*post_LW_Dia$b_Contrast_LWCO_LW_Block) + post_LW_Dia$sigma/2)
LW_MI_I_CO <- exp((post_LW_Dia$b_Intercept - 0.5* post_LW_Dia$b_Contrast_LWCO_Congruency - 0.5 * post_LW_Dia$b_Contrast_LWCO_Listwide  + 0.25*post_LW_Dia$b_Contrast_LWCO_LW_Block) + post_LW_Dia$sigma/2)

# Mean und Interval der Marginal Means CO
c(condition = 'LW_MC_C_CO', mean = mean(LW_MC_C_CO), quantile(LW_MC_C_CO, probs = c(.025, .975)))
c(condition = 'LW_MC_I_CO', mean = mean(LW_MC_I_CO), quantile(LW_MC_I_CO, probs = c(.025, .975)))
c(condition = 'LW_MI_C_CO', mean = mean(LW_MI_C_CO), quantile(LW_MI_C_CO, probs = c(.025, .975)))
c(condition = 'LW_MI_I_CO', mean = mean(LW_MI_I_CO), quantile(LW_MI_I_CO, probs = c(.025, .975)))


# Mit den Marginal Means CO
Stroop_MC_CO_LW_Dia <- (LW_MC_I_CO - LW_MC_C_CO)
c(condition = 'Stroop_MC_CO_LW_Dia', mean = mean(Stroop_MC_CO_LW_Dia), quantile(Stroop_MC_CO_LW_Dia, probs = c(.025, .975)))
Stroop_MI_CO_LW_Dia<- (LW_MI_I_CO - LW_MI_C_CO)
c(condition = 'Stroop_MI_CO_LW_Dia', mean = mean(Stroop_MI_CO_LW_Dia), quantile(Stroop_MI_CO_LW_Dia, probs = c(.025, .975)))
LW_Block_CO_Dia <- Stroop_MC_CO_LW_Dia - Stroop_MI_CO_LW_Dia
c(condition = 'Interaction LW_Block', mean = mean(LW_Block_CO_Dia), quantile(LW_Block_CO_Dia, probs = c(.025, .975)))

remove('LW_MC_C_CO', 'LW_MC_I_CO', 'LW_MI_C_CO', 'LW_MI_I_CO', 'LW_MC_C_PD', 'LW_MC_I_PD', 'LW_MI_C_PD', 'LW_MI_I_PD' )


#### Effekte Inducer IS ####

# Haupteffekt Kongruenz ausrechnen (ndt Parameter wird nicht exponentiert und rechnet)
Congruency_ms <- exp((post_IS_Ind$b_Intercept - 0.5* post_IS_Ind$b_Contrast_LWPD_Congruency)+ post_IS_Ind$sigma/2)  - exp((post_IS_Ind$b_Intercept + 0.5* post_IS_Ind$b_Contrast_LWPD_Congruency) + post_IS_Ind$sigma/2)
c(mean = mean(Congruency_ms), quantile(Congruency_ms, probs = c(.025, .975)))

# Marginal Means ausrechnen PD
IS_MC_C_PD <- exp((post_IS_Ind$b_Intercept + 0.5* post_IS_Ind$b_Contrast_ISPD_Congruency + 0.5 * post_IS_Ind$b_Contrast_ISPD_Itemspecific  + 0.25*post_IS_Ind$b_Contrast_ISPD_IS_Block) + post_IS_Ind$sigma/2)
IS_MC_I_PD <- exp((post_IS_Ind$b_Intercept - 0.5* post_IS_Ind$b_Contrast_ISPD_Congruency + 0.5 * post_IS_Ind$b_Contrast_ISPD_Itemspecific  - 0.25*post_IS_Ind$b_Contrast_ISPD_IS_Block) + post_IS_Ind$sigma/2)
IS_MI_C_PD <- exp((post_IS_Ind$b_Intercept + 0.5* post_IS_Ind$b_Contrast_ISPD_Congruency - 0.5 * post_IS_Ind$b_Contrast_ISPD_Itemspecific  - 0.25*post_IS_Ind$b_Contrast_ISPD_IS_Block) + post_IS_Ind$sigma/2)
IS_MI_I_PD <- exp((post_IS_Ind$b_Intercept - 0.5* post_IS_Ind$b_Contrast_ISPD_Congruency - 0.5 * post_IS_Ind$b_Contrast_ISPD_Itemspecific  + 0.25*post_IS_Ind$b_Contrast_ISPD_IS_Block) + post_IS_Ind$sigma/2)

# Mean und Interval der Marginal Means PD
c(condition = 'IS_MC_C_PD', mean = mean(IS_MC_C_PD), quantile(IS_MC_C_PD, probs = c(.025, .975)))
c(condition = 'IS_MC_I_PD', mean = mean(IS_MC_I_PD), quantile(IS_MC_I_PD, probs = c(.025, .975)))
c(condition = 'IS_MI_C_PD', mean = mean(IS_MI_C_PD), quantile(IS_MI_C_PD, probs = c(.025, .975)))
c(condition = 'IS_MI_I_PD', mean = mean(IS_MI_I_PD), quantile(IS_MI_I_PD, probs = c(.025, .975)))


# Mit den Marginal Means PD
Stroop_MC_PD_IS_Ind <- (IS_MC_I_PD - IS_MC_C_PD)
c(condition = 'Stroop_MC_PD_IS_Ind', mean = mean(Stroop_MC_PD_IS_Ind), quantile(Stroop_MC_PD_IS_Ind, probs = c(.025, .975)))
Stroop_MI_PD_IS_Ind <- (IS_MI_I_PD - IS_MI_C_PD)
c(condition = 'Stroop_MI_PD_IS_Ind', mean = mean(Stroop_MI_PD_IS_Ind), quantile(Stroop_MI_PD_IS_Ind, probs = c(.025, .975)))
IS_Block_PD_Ind <- Stroop_MC_PD_IS_Ind - Stroop_MI_PD_IS_Ind
c(condition = 'Interaction IS_Block', mean = mean(IS_Block_PD_Ind), quantile(IS_Block_PD_Ind, probs = c(.025, .975)))


# Marginal Means ausrechnen CO
IS_MC_C_CO <- exp((post_IS_Ind$b_Intercept + 0.5* post_IS_Ind$b_Contrast_ISCO_Congruency + 0.5 * post_IS_Ind$b_Contrast_ISCO_Itemspecific  + 0.25*post_IS_Ind$b_Contrast_ISCO_IS_Block) + post_IS_Ind$sigma/2)
IS_MC_I_CO <- exp((post_IS_Ind$b_Intercept - 0.5* post_IS_Ind$b_Contrast_ISCO_Congruency + 0.5 * post_IS_Ind$b_Contrast_ISCO_Itemspecific  - 0.25*post_IS_Ind$b_Contrast_ISCO_IS_Block) + post_IS_Ind$sigma/2)
IS_MI_C_CO <- exp((post_IS_Ind$b_Intercept + 0.5* post_IS_Ind$b_Contrast_ISCO_Congruency - 0.5 * post_IS_Ind$b_Contrast_ISCO_Itemspecific  - 0.25*post_IS_Ind$b_Contrast_ISCO_IS_Block) + post_IS_Ind$sigma/2)
IS_MI_I_CO <- exp((post_IS_Ind$b_Intercept - 0.5* post_IS_Ind$b_Contrast_ISCO_Congruency - 0.5 * post_IS_Ind$b_Contrast_ISCO_Itemspecific  + 0.25*post_IS_Ind$b_Contrast_ISCO_IS_Block) + post_IS_Ind$sigma/2)

# Mean und Interval der Marginal Means CO
c(condition = 'IS_MC_C_CO', mean = mean(IS_MC_C_CO), quantile(IS_MC_C_CO, probs = c(.025, .975)))
c(condition = 'IS_MC_I_CO', mean = mean(IS_MC_I_CO), quantile(IS_MC_I_CO, probs = c(.025, .975)))
c(condition = 'IS_MI_C_CO', mean = mean(IS_MI_C_CO), quantile(IS_MI_C_CO, probs = c(.025, .975)))
c(condition = 'IS_MI_I_CO', mean = mean(IS_MI_I_CO), quantile(IS_MI_I_CO, probs = c(.025, .975)))


# Mit den Marginal Means CO
Stroop_MC_CO_IS_Ind <- (IS_MC_I_CO - IS_MC_C_CO)
c(condition = 'Stroop_MC_CO_IS_Ind', mean = mean(Stroop_MC_CO_IS_Ind), quantile(Stroop_MC_CO_IS_Ind, probs = c(.025, .975)))
Stroop_MI_CO_IS_Ind<- (IS_MI_I_CO - IS_MI_C_CO)
c(condition = 'Stroop_MI_CO_IS_Ind', mean = mean(Stroop_MI_CO_IS_Ind), quantile(Stroop_MI_CO_IS_Ind, probs = c(.025, .975)))
IS_Block_CO_Ind <- Stroop_MC_CO_IS_Ind - Stroop_MI_CO_IS_Ind
c(condition = 'Interaction IS_Block', mean = mean(IS_Block_CO_Ind), quantile(IS_Block_CO_Ind, probs = c(.025, .975)))

remove('IS_MC_C_CO', 'IS_MC_I_CO', 'IS_MI_C_CO', 'IS_MI_I_CO', 'IS_MC_C_PD', 'IS_MC_I_PD', 'IS_MI_C_PD', 'IS_MI_I_PD' )


#### Effekte Diagnostics IS ####

# Haupteffekt Kongruenz ausrechnen (ndt Parameter wird nicht exponentiert und rechnet)
Congruency_ms <- exp((post_IS_Dia$b_Intercept - 0.5* post_IS_Dia$b_Contrast_LWPD_Congruency)+ post_IS_Dia$sigma/2)  - exp((post_IS_Dia$b_Intercept + 0.5* post_IS_Dia$b_Contrast_LWPD_Congruency) + post_IS_Dia$sigma/2)
c(mean = mean(Congruency_ms), quantile(Congruency_ms, probs = c(.025, .975)))

# Marginal Means ausrechnen PD
IS_MC_C_PD <- exp((post_IS_Dia$b_Intercept + 0.5* post_IS_Dia$b_Contrast_ISPD_Congruency + 0.5 * post_IS_Dia$b_Contrast_ISPD_Itemspecific  + 0.25*post_IS_Dia$b_Contrast_ISPD_IS_Block) + post_IS_Dia$sigma/2)
IS_MC_I_PD <- exp((post_IS_Dia$b_Intercept - 0.5* post_IS_Dia$b_Contrast_ISPD_Congruency + 0.5 * post_IS_Dia$b_Contrast_ISPD_Itemspecific  - 0.25*post_IS_Dia$b_Contrast_ISPD_IS_Block) + post_IS_Dia$sigma/2)
IS_MI_C_PD <- exp((post_IS_Dia$b_Intercept + 0.5* post_IS_Dia$b_Contrast_ISPD_Congruency - 0.5 * post_IS_Dia$b_Contrast_ISPD_Itemspecific  - 0.25*post_IS_Dia$b_Contrast_ISPD_IS_Block) + post_IS_Dia$sigma/2)
IS_MI_I_PD <- exp((post_IS_Dia$b_Intercept - 0.5* post_IS_Dia$b_Contrast_ISPD_Congruency - 0.5 * post_IS_Dia$b_Contrast_ISPD_Itemspecific  + 0.25*post_IS_Dia$b_Contrast_ISPD_IS_Block) + post_IS_Dia$sigma/2)

# Mean und Interval der Marginal Means PD
c(condition = 'IS_MC_C_PD', mean = mean(IS_MC_C_PD), quantile(IS_MC_C_PD, probs = c(.025, .975)))
c(condition = 'IS_MC_I_PD', mean = mean(IS_MC_I_PD), quantile(IS_MC_I_PD, probs = c(.025, .975)))
c(condition = 'IS_MI_C_PD', mean = mean(IS_MI_C_PD), quantile(IS_MI_C_PD, probs = c(.025, .975)))
c(condition = 'IS_MI_I_PD', mean = mean(IS_MI_I_PD), quantile(IS_MI_I_PD, probs = c(.025, .975)))


# Mit den Marginal Means PD
Stroop_MC_PD_IS_Dia <- (IS_MC_I_PD - IS_MC_C_PD)
c(condition = 'Stroop_MC_PD_IS_Dia', mean = mean(Stroop_MC_PD_IS_Dia), quantile(Stroop_MC_PD_IS_Dia, probs = c(.025, .975)))
Stroop_MI_PD_IS_Dia <- (IS_MI_I_PD - IS_MI_C_PD)
c(condition = 'Stroop_MI_PD_IS_Dia', mean = mean(Stroop_MI_PD_IS_Dia), quantile(Stroop_MI_PD_IS_Dia, probs = c(.025, .975)))
IS_Block_PD_Dia <- Stroop_MC_PD_IS_Dia - Stroop_MI_PD_IS_Dia
c(condition = 'Interaction IS_Block', mean = mean(IS_Block_PD_Dia), quantile(IS_Block_PD_Dia, probs = c(.025, .975)))


# Marginal Means ausrechnen CO
IS_MC_C_CO <- exp((post_IS_Dia$b_Intercept + 0.5* post_IS_Dia$b_Contrast_ISCO_Congruency + 0.5 * post_IS_Dia$b_Contrast_ISCO_Itemspecific  + 0.25*post_IS_Dia$b_Contrast_ISCO_IS_Block) + post_IS_Dia$sigma/2)
IS_MC_I_CO <- exp((post_IS_Dia$b_Intercept - 0.5* post_IS_Dia$b_Contrast_ISCO_Congruency + 0.5 * post_IS_Dia$b_Contrast_ISCO_Itemspecific  - 0.25*post_IS_Dia$b_Contrast_ISCO_IS_Block) + post_IS_Dia$sigma/2)
IS_MI_C_CO <- exp((post_IS_Dia$b_Intercept + 0.5* post_IS_Dia$b_Contrast_ISCO_Congruency - 0.5 * post_IS_Dia$b_Contrast_ISCO_Itemspecific  - 0.25*post_IS_Dia$b_Contrast_ISCO_IS_Block) + post_IS_Dia$sigma/2)
IS_MI_I_CO <- exp((post_IS_Dia$b_Intercept - 0.5* post_IS_Dia$b_Contrast_ISCO_Congruency - 0.5 * post_IS_Dia$b_Contrast_ISCO_Itemspecific  + 0.25*post_IS_Dia$b_Contrast_ISCO_IS_Block) + post_IS_Dia$sigma/2)

# Mean und Interval der Marginal Means CO
c(condition = 'IS_MC_C_CO', mean = mean(IS_MC_C_CO), quantile(IS_MC_C_CO, probs = c(.025, .975)))
c(condition = 'IS_MC_I_CO', mean = mean(IS_MC_I_CO), quantile(IS_MC_I_CO, probs = c(.025, .975)))
c(condition = 'IS_MI_C_CO', mean = mean(IS_MI_C_CO), quantile(IS_MI_C_CO, probs = c(.025, .975)))
c(condition = 'IS_MI_I_CO', mean = mean(IS_MI_I_CO), quantile(IS_MI_I_CO, probs = c(.025, .975)))


# Mit den Marginal Means CO
Stroop_MC_CO_IS_Dia <- (IS_MC_I_CO - IS_MC_C_CO)
c(condition = 'Stroop_MC_CO_IS_Dia', mean = mean(Stroop_MC_CO_IS_Dia), quantile(Stroop_MC_CO_IS_Dia, probs = c(.025, .975)))
Stroop_MI_CO_IS_Dia<- (IS_MI_I_CO - IS_MI_C_CO)
c(condition = 'Stroop_MI_CO_IS_Dia', mean = mean(Stroop_MI_CO_IS_Dia), quantile(Stroop_MI_CO_IS_Dia, probs = c(.025, .975)))
IS_Block_CO_Dia <- Stroop_MC_CO_IS_Dia - Stroop_MI_CO_IS_Dia
c(condition = 'Interaction IS_Block', mean = mean(IS_Block_CO_Dia), quantile(IS_Block_CO_Dia, probs = c(.025, .975)))

remove('IS_MC_C_CO', 'IS_MC_I_CO', 'IS_MI_C_CO', 'IS_MI_I_CO', 'IS_MC_C_PD', 'IS_MC_I_PD', 'IS_MI_C_PD', 'IS_MI_I_PD' )

#### Plot Marginal Mean Effects Proactive Inducer ####

# Get variables into a tibble
MarginalMeans_LW_Ind <- tibble(Stroop_MC_PD_LW_Ind, Stroop_MI_PD_LW_Ind, LW_Block_PD_Ind,Stroop_MC_CO_LW_Ind, Stroop_MI_CO_LW_Ind, LW_Block_CO_Ind)
MarginalMeans_LW_Ind <- MarginalMeans_LW_Ind %>%
  rename(
    Conflict_mCongruent_PD= Stroop_MC_PD_LW_Ind,
    Conflict_mIncongruent_PD= Stroop_MI_PD_LW_Ind ,
    ProactiveC_PD = LW_Block_PD_Ind,
    Conflict_mCongruent_CO = Stroop_MC_CO_LW_Ind,
    Conflict_mIncongruent_CO = Stroop_MI_CO_LW_Ind,
    ProactiveC_CO = LW_Block_CO_Ind
  )

plot_title <- ggtitle("Proactive Control Inducer Items Posterior distributions ",
                      "Marginal Mean Effect Differences and 90% intervals")


vals   <- c(-20, 20, 40, 60, 80, 100)

mcmc_areas(MarginalMeans_LW_Ind,
           pars = c("Conflict_mCongruent_PD", "Conflict_mIncongruent_PD", "ProactiveC_PD","Conflict_mCongruent_CO", "Conflict_mIncongruent_CO", "ProactiveC_CO"),
           prob = 0.9,
           prob_outer = 0.99, # 99%
           point_est = "median"
           ) + vline_0(size = 0.25,color = "darkgray", linetype = 2) + 
  vline_at(vals, linetype = 3, size = 0.25) + plot_title + xlim(c(-20, 120))# +
                                                                  #ggplot2::labs(x = "Average Time in ms "))
                                                       # +
                                                      # ggsave("mm_LW.png"))

#### Plot Marginal Mean Effects Proactive Diagnostics ####

# Get variables into a tibble
MarginalMeans_LW_Dia <- tibble(Stroop_MC_PD_LW_Dia, Stroop_MI_PD_LW_Dia, LW_Block_PD_Dia,Stroop_MC_CO_LW_Dia, Stroop_MI_CO_LW_Dia, LW_Block_CO_Dia)
MarginalMeans_LW_Dia <- MarginalMeans_LW_Dia %>%
  rename(
    Conflict_mCongruent_PD= Stroop_MC_PD_LW_Dia,
    Conflict_mIncongruent_PD= Stroop_MI_PD_LW_Dia ,
    ProactiveC_PD = LW_Block_PD_Dia,
    Conflict_mCongruent_CO = Stroop_MC_CO_LW_Dia,
    Conflict_mIncongruent_CO = Stroop_MI_CO_LW_Dia,
    ProactiveC_CO = LW_Block_CO_Dia
  )

plot_title <- ggtitle("Proactive Control Diagnostic Items Posterior distributions ",
                      "Marginal Mean Effect Differences and 90% intervals")


vals   <- c(-20, 20, 40, 60, 80, 100)

mcmc_areas(MarginalMeans_LW_Dia,
           pars = c("Conflict_mCongruent_PD", "Conflict_mIncongruent_PD", "ProactiveC_PD","Conflict_mCongruent_CO", "Conflict_mIncongruent_CO", "ProactiveC_CO"),
           prob = 0.9,
           prob_outer = 0.99, # 99%
           point_est = "median"
) + vline_0(size = 0.25,color = "darkgray", linetype = 2) + 
  vline_at(vals, linetype = 3, size = 0.25) + plot_title + xlim(c(-20, 120))

#### Plot Marginal Mean Effects Reactive Inducer ####

# Get variables into a tibble
MarginalMeans_IS_Ind <- tibble(Stroop_MC_PD_IS_Ind, Stroop_MI_PD_IS_Ind, IS_Block_PD_Ind,Stroop_MC_CO_IS_Ind, Stroop_MI_CO_IS_Ind, IS_Block_CO_Ind)
MarginalMeans_IS_Ind <- MarginalMeans_IS_Ind %>%
  rename(
    Conflict_mCongruent_PD= Stroop_MC_PD_IS_Ind,
    Conflict_mIncongruent_PD= Stroop_MI_PD_IS_Ind ,
    ReactiveC_PD = IS_Block_PD_Ind,
    Conflict_mCongruent_CO = Stroop_MC_CO_IS_Ind,
    Conflict_mIncongruent_CO = Stroop_MI_CO_IS_Ind,
    ReactiveC_CO = IS_Block_CO_Ind
  )

plot_title <- ggtitle("Reactive Control Inducer Items Posterior distributions ",
                      "Marginal Mean Effect Differences and 90% intervals")


vals   <- c(-20, 20, 40, 60, 80, 100)

mcmc_areas(MarginalMeans_IS_Ind,
           pars = c("Conflict_mCongruent_PD", "Conflict_mIncongruent_PD", "ReactiveC_PD","Conflict_mCongruent_CO", "Conflict_mIncongruent_CO", "ReactiveC_CO"),
           prob = 0.9,
           prob_outer = 0.99, # 99%
           point_est = "median"
) + vline_0(size = 0.25,color = "darkgray", linetype = 2) + 
  vline_at(vals, linetype = 3, size = 0.25) + plot_title + xlim(c(-20, 120))

#### Plot Marginal Mean Effects Reactive Diagnostics ####

# Get variables into a tibble
MarginalMeans_IS_Dia <- tibble(Stroop_MC_PD_IS_Dia, Stroop_MI_PD_IS_Dia, IS_Block_PD_Dia,Stroop_MC_CO_IS_Dia, Stroop_MI_CO_IS_Dia, IS_Block_CO_Dia)
MarginalMeans_IS_Dia <- MarginalMeans_IS_Dia %>%
  rename(
    Conflict_mCongruent_PD= Stroop_MC_PD_IS_Dia,
    Conflict_mIncongruent_PD= Stroop_MI_PD_IS_Dia ,
    ReactiveC_PD = IS_Block_PD_Dia,
    Conflict_mCongruent_CO = Stroop_MC_CO_IS_Dia,
    Conflict_mIncongruent_CO = Stroop_MI_CO_IS_Dia,
    ReactiveC_CO = IS_Block_CO_Dia
  )

plot_title <- ggtitle("Reactive Control Diagnostic Items Posterior distributions ",
                      "Marginal Mean Effect Differences and 90% intervals")


vals   <- c(-20, 20, 40, 60, 80, 100)

mcmc_areas(MarginalMeans_IS_Dia,
           pars = c("Conflict_mCongruent_PD", "Conflict_mIncongruent_PD", "ReactiveC_PD","Conflict_mCongruent_CO", "Conflict_mIncongruent_CO", "ReactiveC_CO"),
           prob = 0.9,
           prob_outer = 0.99, # 99%
           point_est = "median"
) + vline_0(size = 0.25,color = "darkgray", linetype = 2) + 
  vline_at(vals, linetype = 3, size = 0.25) + plot_title + xlim(c(-20, 120))

#### Get a different plotting Layout ####

#### include the data of the IS effect as rows into a tibble ####

var_names <- c("Stroop_MC_PD_IS_Dia", "Stroop_MI_PD_IS_Dia", "IS_Block_PD_Dia", "Stroop_MC_PD_IS_Ind", "Stroop_MI_PD_IS_Ind", "IS_Block_PD_Ind", 
               "Stroop_MC_CO_IS_Dia", "Stroop_MI_CO_IS_Dia", "IS_Block_CO_Dia", "Stroop_MC_CO_IS_Ind", "Stroop_MI_CO_IS_Ind", "IS_Block_CO_Ind",
               "Stroop_MC_PD_LW_Dia", "Stroop_MI_PD_LW_Dia", "LW_Block_PD_Dia", "Stroop_MC_PD_LW_Ind", "Stroop_MI_PD_LW_Ind", "LW_Block_PD_Ind", 
               "Stroop_MC_CO_LW_Dia", "Stroop_MI_CO_LW_Dia", "LW_Block_CO_Dia", "Stroop_MC_CO_LW_Ind", "Stroop_MI_CO_LW_Ind", "LW_Block_CO_Ind")

Marginal_Effects <- tibble()

for (v_n in var_names) {
  temp_data <- eval(parse(text = v_n))
  temp <- tibble(temp_data,
                 Var_name = rep(v_n, length.out = length(Stroop_MC_PD_IS_Ind)),
                 Item_type = rep(str_detect(v_n, "Ind", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind)),
                 Group = rep(str_detect(v_n, "PD", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind)),
                 Effect = rep(str_detect(v_n, "LW", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind)),
                 MC = rep(str_detect(v_n, "MC", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind)),
                 MI = rep(str_detect(v_n, "MI", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind)),
                 Block = rep(str_detect(v_n, "Block", negate = FALSE), length.out = length(Stroop_MC_PD_IS_Ind))
                 ) %>%
                 
    rename(Mean_RT= temp_data) %>%
    mutate(Group = if_else(Group == TRUE, "Parkinson Participants", "Control Participants", missing = NULL)) %>%
    mutate(Item_type = if_else(Item_type == TRUE, "Inducer", "Diagnostic", missing = NULL)) %>%
    mutate(Effect = if_else(Effect == TRUE, "Proactive Control", "Reactive Control", missing = NULL)) %>%
    mutate(Measure = case_when(
      MC == TRUE  ~ "MC",
      MI == TRUE  ~ "MI",
      Block == TRUE ~ "Effect"))
  
  Marginal_Effects <- bind_rows(Marginal_Effects, temp)
  rm(temp)
}

# remove variables which we no longer need
rm(Stroop_MC_PD_IS_Dia, Stroop_MI_PD_IS_Dia, IS_Block_PD_Dia, Stroop_MC_CO_IS_Dia, Stroop_MI_CO_IS_Dia, IS_Block_CO_Dia,
   Stroop_MC_PD_LW_Ind, Stroop_MI_PD_LW_Ind, LW_Block_PD_Ind, Stroop_MC_CO_LW_Ind, Stroop_MI_CO_LW_Ind, LW_Block_CO_Ind)
x_max <- 110 
x_min <- -10

display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 12, name = "Greens")
brewer.pal(n = 12, name = "Greens")

ggplot(NULL, aes(x = Mean_RT, y = Measure, fill = Item_type)) + facet_grid(Effect ~ Group)  + 
  theme_bw(base_size = 14) +
  stat_halfeye(data = Marginal_Effects, side = "right", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 120, by = 20), limits = c(-10, 120)) +
  labs(x = "RT in ms", fill = "Item Type") +
  ggtitle("Marginal Mean Effects Reaction Time") + 
  theme(
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_fill_manual(values = c("#B2182B","#2166AC"))
  # scale_fill_manual(values = c("#B2182B","#2166AC")) +
  # theme(
  #   panel.background = element_rect(colour = "white", fill = "white", size = 0.5, color = "white"),
  #   panel.grid.major = element_line(size = 0.5, colour = "gray"),
  #   panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "gray"),
  #   panel.border = 
  #       )
  
  #stat_halfeye(data = Marginal_Effects %>% filter(Item_type == "Diagnostic"), position = "dodge", side = "left", alpha = 0.5) +
  #stat_halfeye(data = Marginal_Effects %>% filter(Item_type == "Inducer"), position = "dodge", side = "right", alpha = 0.5) 
  
