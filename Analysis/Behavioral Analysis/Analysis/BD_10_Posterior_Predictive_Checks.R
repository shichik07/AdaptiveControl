#####
# Author: Julius Kricheldorff
# Analysis Behavioral Data - Posterior Predictive Checks

####

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
library(plotrix) # to calculate standard error
library(forcats)
library(ggdist)
library(hypr)
library(ggthemes)

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

# dir were we save our models
Full_models_saveloc <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data"

# function to load the models
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

#functions to generate posterior predicitve checks

# generate mean values
post_pred <- function(data2pred, model, number_pred, effect){
 
  # get posterior predictions
  pp_data <- posterior_predict(model, ndraws = number_pred, newdata = data2pred)
  # calculate the mean
  mean_val <- mean(pp_data)
  quant_val <- quantile(pp_data, probs = c(0.1, 0.9))
  return(list(mean_val = mean_val, quant_val = quant_val))
}

# create tibble with predicted, and actual data
pp_df <-  function(data, model, number_pred, effect){
  
  #create a new variable so that we can use this for any effect
   if(effect == "LWPCE"){
     data <- data %>%
       mutate(Contrast = Contrast_LW)
   } else if (effect == "ISPCE"){
     data <- data %>%
       mutate(Contrast = Contrast_IS)
     
   }
  # create a tibble where we save our data
  pp_summary <- tibble(RT = numeric(),
                       Congruency = character(),
                       Group = character(),
                       Item_specific = character(),
                       Analysis_type = character(),
                       Item_type = character(), 
                       data_type = character(),
                       Measure = character())
  
  for (lev in levels(data$Contrast)) {
    print(lev)
    tmp_data <- data %>% dplyr::filter(Contrast == lev)
    pp_res <- post_pred(data2pred = tmp_data, model = model, number_pred = 1000)
    
    # create a temporary tibble for the mean
    temp_tib <- tibble(RT = pp_res$mean_val,
                       Congruency = unique(tmp_data$Congruency),
                       Group = as.character(unique(tmp_data$Exp_group)),
                       Item_specific = unique(tmp_data$Item_specific),
                       Analysis_type = unique(tmp_data$Analysis_type),
                       Item_type =unique(tmp_data$Trl_type),
                       data_type = "predicted",
                       Measure = "mean")
    
    pp_summary <- bind_rows(pp_summary, temp_tib)
    # and for the quantiles
    temp_tib$RT <-  pp_res$quant_val[[2]]
    temp_tib$Measure <- "Q90"
    pp_summary <- bind_rows(pp_summary, temp_tib)          
    temp_tib$RT <-  pp_res$quant_val[[1]]
    temp_tib$Measure <- "Q10"
    pp_summary <- bind_rows(pp_summary, temp_tib) 
    #and the same for the actual data
    temp_tib$RT <- mean(tmp_data$RT)
    temp_tib$Measure <- "mean"
    temp_tib$data_type = "observed"
    pp_summary <- bind_rows(pp_summary, temp_tib)        
    temp_tib$RT <-  quantile(tmp_data$RT, probs = 0.9)[[1]]
    temp_tib$Measure <- "Q90"
    pp_summary <- bind_rows(pp_summary, temp_tib)          
    temp_tib$RT <-  quantile(tmp_data$RT, probs = 0.1)[[1]]
    temp_tib$Measure <- "Q10"
    pp_summary <- bind_rows(pp_summary, temp_tib) 
    
    
                       
                       
    #                    
    # postpred_90 = pp_res$quant_val[[2]],
    # postpred_10 = pp_res$quant_val[[1]],
    # data_mean = mean(tmp_data$RT),
    # data_90 = quantile(tmp_data$RT, probs = 0.9)[[1]],
    # data_10 = quantile(tmp_data$RT, probs = 0.1)[[1]],
    # Congruency = unique(tmp_data$Congruency),
    # Group = as.character(unique(tmp_data$Exp_group)),
    # Item_specific = unique(tmp_data$Item_specific),
    # Analysis_type = unique(tmp_data$Analysis_type),
    # Item_type =unique(tmp_data$Trl_type)
    #)
    pp_summary <- bind_rows(pp_summary, temp_tib)
  }
  return(pp_summary)
}



# check Inducer LW model
LW_ind <- load_full_mod(loc = Full_models_saveloc,
                       model_t = "RT",
                       item_type = "inducer",
                       effect = "LW")

PP_data_LW_ind <-  pp_df(data = Data_y_inducer_LWPC, model = LW_ind, number_pred = 1000,effect = "LWPCE")


# check diagnostic LW model
LW_dia <- load_full_mod(loc = Full_models_saveloc,
                        model_t = "RT",
                        item_type = "diagnostic",
                        effect = "LW")
PP_data_LW_dia <-  pp_df(data = Data_y_diagnostic_LWPC, model = LW_dia, number_pred = 1000,effect = "LWPCE")


# check Inducer IS model
IS_ind <- load_full_mod(loc = Full_models_saveloc,
                        model_t = "RT",
                        item_type = "inducer",
                        effect = "IS")
PP_data_IS_ind <-  pp_df(data = Data_y_inducer_ISPC, model = IS_ind, number_pred = 1000,effect = "ISPCE")


# check diagnostic LW model
IS_dia <- load_full_mod(loc = Full_models_saveloc,
                        model_t = "RT",
                        item_type = "diagnostic",
                        effect = "IS")
PP_data_IS_dia <-  pp_df(data = Data_y_diagnostic_ISPC, model = IS_dia, number_pred = 1000,effect = "ISPCE")

# join the four posterior tibbles
post_sum <- bind_rows(PP_data_LW_ind, PP_data_LW_dia, PP_data_IS_ind, PP_data_IS_dia) %>%
  mutate(dot_size = case_when(
    data_type == "predicted" ~ 4,
    data_type == "observed" ~5
  ))

# determine save path
save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"

# Plot LWPCE effects
LWPCE <- post_sum %>% filter(Item_specific == "List_wide") %>% mutate(Group = as.factor(Group))
levels(LWPCE$Group) <- list(HC = "CO", PD = "PD")
pp_LWPCE <- ggplot(data = LWPCE, aes(x = Analysis_type, y= RT, color = Congruency, group = Congruency, shape = data_type)) +
  facet_grid(Item_type ~ Group) +
  theme_clean(base_family = "Zilla Slab", base_size = 14) +
  geom_point(size = LWPCE$dot_size,  position = position_dodge(width = 0.5), alpha = 0.7) +
  scale_shape_manual(values = c(1,16)) +
  coord_cartesian(ylim = c(0, 1200)) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  labs(title ="Posterior predicitive checks: LWPCE", x = " Conflict Proportion", y = "RT in ms", 
       shape = "Data") +
  theme(
    panel.spacing = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel",
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 6), nrow = 2),
         shape = guide_legend(override.aes = list(size = 6), nrow = 2))

save_n <- "PosteriorPredictive_LWPCE.png"

ggsave(path = save_path,
       plot = pp_LWPCE, 
       filename = save_n,  
       dpi=300, 
       width = 16.5,
       units = c("cm"))


# Now the same for ISPCE 
ISPCE <- post_sum %>% filter(Item_specific != "List_wide") %>%
  mutate(Analysis_type = case_when(
    Analysis_type == "main_con" ~ "MC",
    Analysis_type == "main_incon" ~ "MI"))%>% 
  mutate(Group = as.factor(Group))
levels(ISPCE$Group) <- list(HC = "CO", PD = "PD")
  
pp_ISPCE <-ggplot(data = ISPCE, aes(x = Analysis_type, y= RT, color = Congruency, group = Congruency, shape = data_type)) +
  facet_grid(Item_type ~ Group) +
  theme_clean(base_family = "Zilla Slab", base_size = 14) +
  geom_point(size = LWPCE$dot_size,  position = position_dodge(width = 0.5), alpha = 0.7) +
  scale_shape_manual(values = c(1,16)) +
  coord_cartesian(ylim = c(0, 1200)) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  labs(title ="Posterior predicitive checks: ISPCE", x = " Conflict Proportion", y = "RT in ms", 
       shape = "Data") +
  theme(
    panel.spacing = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel",
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 6), nrow = 2),
         shape = guide_legend(override.aes = list(size = 6), nrow = 2))

save_n <- "PosteriorPredictive_ISPCE.png"

ggsave(path = save_path,
       plot = pp_ISPCE, 
       filename = save_n,  
       dpi=300, 
       width = 16.5,
       units = c("cm"))


#### For the random slop models

a <- load("fit_ShiftedLog_COOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic # Load data Listwise Inducer

LW_dia_CO <- fit_ShiftedLog_COOnly_diagnostic_LW

pp_check(LW_dia_CO, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(LW_dia_CO, ndraws = 100, type = "dens_overlay")
# Looks a little off
pp_check(LW_dia_CO, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(LW_dia_CO, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(LW_dia_CO, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_LW")
# looks very good

a <- load("fit_ShiftedLog_PDOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic # Load data Listwise Inducer

LW_dia_PD <- fit_ShiftedLog_PDOnly_diagnostic_LW

pp_check(LW_dia_PD, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(LW_dia_PD, ndraws = 100, type = "dens_overlay")
# Looks very good
pp_check(LW_dia_PD, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(LW_dia_PD, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(LW_dia_PD, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_LW")
# looks very good


# To sum up for me, the random intercept model could be improved.
# But if I look at the estimates of the random slope model for my effects of interest (the interaction),
# they vary by one to two milliseconds and would not change the conclusions drawn. I will not get a substantial
# improvement by using separate models for each group or including random slopes - thus I will go with the simple model.

