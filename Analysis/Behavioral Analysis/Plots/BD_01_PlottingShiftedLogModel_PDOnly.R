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
library(plotrix) # to calculate standard error
library(haven)
library(ggpubr)
library(forcats) # so we can simply reorder the variables with fct_inorder

# Set a seed for sake of reproducibility
set.seed(32944)

load("fit_ShiftedLog_PDOnly_inducer_LW.rda") # Load data Listwise Inducer
load("fit_ShiftedLog_PDOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic

# posterior samples aus dem Model entnehmen
post_LW_Ind <- posterior_samples(fit_ShiftedLog_PDOnly_inducer_LW)
post_LW_Dia <- posterior_samples(fit_ShiftedLog_PDOnly_diagnostic_LW)

#load the clinical data first
Clinic <- read.table(file = "participants.tsv", header = TRUE)
Clinic <- as_tibble(Clinic)
 
path1 = "Z:/JuliusKricheldorff/AdaptiveControl/BehavioralData"
files = list.files(path = path1, pattern = "AdapCon_PD", full.names = TRUE, recursive = FALSE)
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




## Clean up the data a little bit
data <- df %>% 
  mutate(number_1 = as.character(number_1))%>%
  mutate(number_2 = as.character(number_2))%>%
  unite('Item', number_1:number_2, remove = FALSE) %>%
  mutate(Item = as_factor(Item)) %>%
  mutate(Item_specific = ifelse(Analysis_type == "main_con" | Analysis_type == "main_incon", 'Item_spec', 'List_wide' )) %>%
  mutate( Congruency_main = ifelse(Analysis_type == "main_con"| Analysis_type == "MC", "M_Congruent", "M_Incongruent")) %>% 
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main, UPDRS) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response) %>%
  filter(Subject != "sub-PD_c(\"11\", \"2\")" &  Subject != "sub-PD_16") #exclude participant 16 and the partial dataset by participant 11



#### Effekte Inducer LW ####
# find unique numbers of particpants
part_nr <- unique(data$Subject)
a <- unique(data$UPDRS)

posterior_cor <- tibble(
  Subject = character(),
  UPDRS = numeric(),
  MC_Stroop = numeric(),
  MI_Stroop = numeric(),
  CC_cost = numeric(),
  CI_cost = numeric(),
  ProactiveC = numeric(),
  Item = character(),
  paired = numeric()
)
p_num <- 0

# Calculate median of the marginal mean effect distributions
for (part in part_nr) {
  p_num <- p_num +1
  # First get the parameternames
  var_intercept <- paste("post_LW_Ind$`r_Subject[",part, ",Intercept]`", sep = "")
  var_Congruency <- paste("post_LW_Ind$`r_Subject[",part, ",Contrast_LWPD_Congruency]`", sep = "")
  var_Listwide <- paste("post_LW_Ind$`r_Subject[",part, ",Contrast_LWPD_Listwide]`", sep = "")
  var_Interaction <- paste("post_LW_Ind$`r_Subject[",part, ",Contrast_LWPD_LW_Block]`", sep = "")
  
  # calculate marginal in the inducer items  
  LW_MC_C_PD_ind <- exp((post_LW_Ind$b_Intercept + eval(parse(text = var_intercept))) +
                      0.5* (post_LW_Ind$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) +
                      0.5 * (post_LW_Ind$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  +
                      0.5* (post_LW_Ind$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                      post_LW_Ind$sigma/2)
  LW_MC_I_PD_ind <- exp((post_LW_Ind$b_Intercept + eval(parse(text = var_intercept))) - 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) + 
                      0.5 * (post_LW_Ind$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                      post_LW_Ind$sigma/2)
  LW_MI_C_PD_ind <- exp((post_LW_Ind$b_Intercept + eval(parse(text = var_intercept))) + 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                      0.5 * (post_LW_Ind$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                      post_LW_Ind$sigma/2)
  LW_MI_I_PD_ind <- exp((post_LW_Ind$b_Intercept + eval(parse(text = var_intercept))) - 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                      0.5 * (post_LW_Ind$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  + 
                      0.5* (post_LW_Ind$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                      post_LW_Ind$sigma/2)
  
  
  # Calculate difference Interaction
  Stroop_MC_PD_LW_Ind <- (LW_MC_I_PD_ind - LW_MC_C_PD_ind)
  Stroop_MI_PD_LW_Ind <- (LW_MI_I_PD_ind - LW_MI_C_PD_ind)
  LW_Block_PD_Ind <- Stroop_MC_PD_LW_Ind - Stroop_MI_PD_LW_Ind
  
  # calculate reduction congruency cost
  CC_cost_Ind <- (LW_MC_C_PD_ind - LW_MI_C_PD_ind)
  CI_cost_Ind <- (LW_MI_I_PD_ind - LW_MC_I_PD_ind)
  
  # calculate marginal in the diagnostic items  
  LW_MC_C_PD_Dia <- exp((post_LW_Dia$b_Intercept + eval(parse(text = var_intercept))) +
                          0.5* (post_LW_Dia$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) +
                          0.5 * (post_LW_Dia$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  +
                          0.5* (post_LW_Dia$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                          post_LW_Dia$sigma/2)
  LW_MC_I_PD_Dia <- exp((post_LW_Dia$b_Intercept + eval(parse(text = var_intercept))) - 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) + 
                          0.5 * (post_LW_Dia$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                          post_LW_Dia$sigma/2)
  LW_MI_C_PD_Dia <- exp((post_LW_Dia$b_Intercept + eval(parse(text = var_intercept))) + 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                          0.5 * (post_LW_Dia$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  - 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                          post_LW_Dia$sigma/2)
  LW_MI_I_PD_Dia <- exp((post_LW_Dia$b_Intercept + eval(parse(text = var_intercept))) - 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_Congruency + eval(parse(text = var_Congruency))) - 
                          0.5 * (post_LW_Dia$b_Contrast_LWPD_Listwide + eval(parse(text = var_Listwide)))  + 
                          0.5* (post_LW_Dia$b_Contrast_LWPD_LW_Block + eval(parse(text = var_Interaction))) + 
                          post_LW_Dia$sigma/2)
  
  
  # Calculate difference Interaction
  Stroop_MC_PD_LW_Dia <- (LW_MC_I_PD_Dia - LW_MC_C_PD_Dia)
  Stroop_MI_PD_LW_Dia <- (LW_MI_I_PD_Dia - LW_MI_C_PD_Dia)
  LW_Block_PD_Dia <- Stroop_MC_PD_LW_Dia - Stroop_MI_PD_LW_Dia
  
  # calculate reduction congruency cost
  CC_cost_Dia <- (LW_MC_C_PD_Dia - LW_MI_C_PD_Dia)
  CI_cost_Dia <- (LW_MI_I_PD_Dia - LW_MC_I_PD_Dia)
  
  #updrs score
  m_score <- pull(data %>% filter(Subject == part) %>%  distinct(UPDRS))
  
  # save inducer variables
  posterior_cor <- posterior_cor %>% add_row(Subject = part, 
                            UPDRS = m_score, 
                            MC_Stroop = mean(Stroop_MC_PD_LW_Ind),
                            MI_Stroop = mean(Stroop_MI_PD_LW_Ind),
                            CC_cost = mean(CC_cost_Ind),
                            CI_cost = mean(CI_cost_Ind),
                            ProactiveC = mean(LW_Block_PD_Ind),
                            Item = "Inducer",
                            paired = p_num
                            )
  
  # save diagnostic variables
  posterior_cor <- posterior_cor %>% add_row(Subject = part, 
                                             UPDRS = m_score, 
                                             MC_Stroop = mean(Stroop_MC_PD_LW_Dia),
                                             MI_Stroop = mean(Stroop_MI_PD_LW_Dia),
                                             CC_cost = mean(CC_cost_Dia),
                                             CI_cost = mean(CI_cost_Dia),
                                             ProactiveC = mean(LW_Block_PD_Dia),
                                             Item = "Diagnostic",
                                             paired = p_num
  )
  
}


# now let us plot as a test
sp <- ggscatter(posterior_cor, x = "ProactiveC", y = "UPDRS",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
                ) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(Item))) +
  labs(y = "UPDRS-III", x= "Proactive Control Effect in ms") +
  ggtitle("Association UPDRS and Individual Marginal Mean Effects - PD Group")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 20, label.y = 25)

# New graph diagnostic to inducer by lines
ggplot(data = posterior_cor, aes(x = fct_inorder(Item), y = ProactiveC, color = Item))  +
  geom_line(aes(group = paired),
            color="black") +
  geom_point(aes(fill=Item),size=3) +
  theme_bw(base_size = 14) +
  labs(x = "Item Type", y = "Proactive Control Effect in ms") +
  ggtitle("Individual Marginal Mean Effects - PD Group") + 
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_fill_manual(values = c("#B2182B","#2166AC"))

posterior_cor %>%
  filter(Item == "Diagnostic") %>%
  summarize(mean_eff = mean(ProactiveC), mean_CC = mean(CC_cost), mean_CI= mean(CI_cost))