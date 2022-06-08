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
library(hypr) # package to generate contrasts
library(chkptstanr) #package to interupt and restart sampling with stanr/brms

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
  Congruency = (LW_incon_C + LW_con_C + IS_incon_C + IS_con_C)/4 ~ (LW_con_I + LW_incon_I + IS_con_I + IS_incon_I)/4, # main effect Congruence healthy controls
  Listwide = (LW_con_C + LW_con_I)/2 ~ (LW_incon_C + LW_incon_I)/2, # main effect Block Listwide
  Itemspecific = (IS_con_C + IS_con_I)/2 ~ (IS_incon_C + IS_incon_I)/2, # main effect Itemspecfic
  LW_Block = (LW_incon_I - LW_incon_C)/2 ~ (LW_con_I - LW_con_C)/2, # interaction listwide effect
  IS_Block = (IS_incon_I - IS_incon_C)/2 ~ (IS_con_I - IS_con_C)/2, # interaction itemspecfic effect
  #levels = c("LW_incon_C", "LW_con_C","IS_incon_C", "IS_con_C", "LW_con_I", "LW_incon_I", "IS_con_I", "IS_incon_I")
  levels = c("LW_con_C", "LW_con_I", "LW_incon_I", "LW_incon_C", "IS_incon_I", "IS_con_I", "IS_con_C", "IS_incon_C")
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
  prior(uniform(0,10), class = sigma), # truncated normal
  prior(uniform(0, 10), class = ndt), # pretty hard to guess what an accurate shift parameter would look like
  prior(normal(0,0.5), class = b, coef = Contrast_FCongruency),
  prior(normal(0,0.5), class = b, coef = Contrast_FListwide),
  prior(normal(0,0.5), class = b, coef = Contrast_FItemspecific),
  prior(normal(0,0.5), class = b, coef = Contrast_FLW_Block),
  prior(normal(0,0.5), class = b, coef = Contrast_FIS_Block),
  prior(normal(0, 3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 3), class = sd, coef = Intercept, group = Item)
)

# brmsformula object
m1 <- bf(RT ~ 1 + Contrast_F + (1|Subject) + (1|Item), family = shifted_lognormal()) #lognormal())

# inducer_young <- chkpt_brms(formula = m1,
#                            data = Data_y_inducer,
#                            prior = prior_informed,
#                            warmup = 2000,
#                            iter = 8000, # 20000 is the limit necessary for bridge sampling
#                            iter_per_chkpt = 250 # after how many samples we save
#                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
# )

fit_inducer_young <-brm(formula = m1,
                            data = Data_y_inducer,
                            prior = prior_informed,
                            warmup = 2000,
                            iter = 8000 # 20000 is the limit necessary for bridge sampli
                            #save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                            )
save(inducer_young, file = "inducer_young.rda")
load(file = "inducer_young.rda")
inducer_young

# Check fit of the posterior distribution in general
pp_check(inducer_young, ndraws = 1000) 
# Seems this distribution is not really an ideal fit.. we will have to see about that again

# check predicted max min and median values versus the actual data
pp_check(inducer_young, type = "stat", stat = "min")
pp_check(inducer_young, type = "stat", stat = "max")
pp_check(inducer_young, type = "stat", stat = "median")


# Plot contrast distributions

mcmc_dens(inducer_young, pars = variables(inducer_young)[1:10])

# regardless, let us calculate the effects in milliseconds first
alpha_samples <- as_draws_df(inducer_young)$b_Intercept
Congruency    <- as_draws_df(inducer_young)$b_Contrast_FCongruency  
Listwide      <- as_draws_df(inducer_young)$b_Contrast_FListwide 
itemspecific  <- as_draws_df(inducer_young)$b_Contrast_FItemspecific
LW_Block      <- as_draws_df(inducer_young)$b_Contrast_FLW_Block
IS_Block      <- as_draws_df(inducer_young)$b_Contrast_FIS_Block

# first congruency effects
congruency_ms <- exp(alpha_samples + 0.5* Congruency) - exp(alpha_samples - 0.5* Congruency)
c(mean = mean(congruency_ms), quantile(congruency_ms, probs = c(.025, .975)))

# next we look at the block effects  which is a bit more tricky - be careful in how you constructed your contrasts!
Interaction_LW_MC_congruent <- exp(alpha_samples + 0.5* Congruency + 0.5 * Listwide  + 0.5*LW_Block)
Interaction_LW_MC_incongruent <- exp(alpha_samples - 0.5* Congruency + 0.5 * Listwide - 0.5*LW_Block)
Interaction_LW_MI_congruent <- exp(alpha_samples + 0.5* Congruency - 0.5 * Listwide  + 0.5*LW_Block)
Interaction_LW_MI_incongruent <- exp(alpha_samples - 0.5* Congruency - 0.5 * Listwide - 0.5*LW_Block)

# Now that we have the marginal means, we can calculate the Interaction effect in ms
LW_Block_ms <- (Interaction_LW_MC_incongruent - Interaction_LW_MC_congruent) - (Interaction_LW_MI_incongruent - Interaction_LW_MI_congruent)
c(mean = mean(LW_Block_ms), quantile(LW_Block_ms, probs = c(.025, .975)))


# Now do the Same exact thing for the Item specific effects
Interaction_IS_MC_congruent <- exp(alpha_samples + 0.5* Congruency + 0.5 * itemspecific  + 0.5*IS_Block)
Interaction_IS_MC_incongruent <- exp(alpha_samples - 0.5* Congruency + 0.5 * itemspecific - 0.5*IS_Block)
Interaction_IS_MI_congruent <- exp(alpha_samples + 0.5* Congruency - 0.5 * itemspecific  + 0.5*IS_Block)
Interaction_IS_MI_incongruent <- exp(alpha_samples - 0.5* Congruency - 0.5 * itemspecific - 0.5*IS_Block)
# And again the interaction effect
IS_Block_ms <- (Interaction_IS_MC_incongruent - Interaction_IS_MC_congruent) - (Interaction_IS_MI_incongruent - Interaction_IS_MI_congruent)
c(mean = mean(IS_Block_ms), quantile(IS_Block_ms, probs = c(.025, .975)))


