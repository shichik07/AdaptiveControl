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
library(cowplot)

# Set a seed for sake of reproducibility
set.seed(32944)

load("fit_ShiftedLog_PDOnly_inducer_LW.rda") # Load data Listwise Inducer
load("fit_ShiftedLog_PDOnly_diagnostic_LW.rda") # Load data Listwise Diagnostic

# posterior samples aus dem Model entnehmen
post_LW_Ind <- posterior_samples(fit_ShiftedLog_PDOnly_inducer_LW)
post_LW_Dia <- posterior_samples(fit_ShiftedLog_PDOnly_diagnostic_LW)

#load the clinical data first
Clinic <- read.table(file = "participants_full1.tsv", header = TRUE)
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
      if (Clinic$educationyears[updrs_ind] == "n/a"){
        t$educationyears = NA
      } else {
        t$educationyears <- as.numeric(Clinic$educationyears[updrs_ind])
      }
      if (Clinic$disease_duration[updrs_ind] == "n/a"){
        t$disease_duration = NA
      } else {
        t$disease_duration <- as.numeric(Clinic$disease_duration[updrs_ind])
      }
      if (Clinic$H_and_Y[updrs_ind] == "n/a"){
        t$H_and_Y = NA
      } else {
        t$H_and_Y <- as.numeric(Clinic$H_and_Y[updrs_ind])
      }
      if (Clinic$Aff_Side[updrs_ind] == "n/a"){
        t$Aff_Side = NA
      } else {
        t$Aff_Side <- as.character(Clinic$Aff_Side[updrs_ind])
      }
      if (Clinic$Subtype[updrs_ind] == "n/a"){
        t$Subtype = NA
      } else {
        t$Subtype <- as.character(Clinic$Subtype[updrs_ind])
      }
      if (Clinic$LEDD[updrs_ind] == "n/a"){
        t$LEDD = NA
      } else {
        t$LEDD <- as.character(Clinic$LEDD[updrs_ind])
      }
      t$MMST <- as.numeric(Clinic$MMST[updrs_ind])
    }
  } else {
    t$UPDRS = NA
    #t$disease_duration = NA
    #t$H_and_Y = NA
    #t$Aff_Side = NA
    #t$Subtype = NA
  }
  if  (substr(subject_id, 5,6) == "CO"){
    updrs_ind <- which(Clinic$participant_id == subject_id)
    if (Clinic$educationyears[updrs_ind] == "n/a"){
      t$educationyears = NA
    } else {
      t$educationyears <- as.numeric(Clinic$educationyears[updrs_ind])
    }
    t$MMST <- as.numeric(Clinic$MMST[updrs_ind])
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
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, 
         Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, 
         Item_specific, Congruency_main, UPDRS, H_and_Y, disease_duration, Aff_Side, Subtype, LEDD) %>%
  mutate(LEDD = as.numeric(LEDD)) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response) %>%
  filter(Subject != "sub-PD_c(\"11\", \"2\")" &  Subject != "sub-PD_16") #exclude participant 16 and the partial dataset by participant 11



#### Effekte Inducer LW ####
# find unique numbers of particpants
part_nr <- unique(data$Subject)
a <- unique(data$UPDRS)

posterior_cor <- tibble(
  Subject = character(),
  UPDRS = numeric(),
  HnY = numeric(),
  Dis_Dur = numeric(),
  LEDD = numeric(),
  Aff_Side = character(),
  Subtype = character(),
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
  
  #Pull clinical scores
  m_score <- pull(data %>% filter(Subject == part) %>%  distinct(UPDRS))
  HnY     <- pull(data %>% filter(Subject == part) %>%  distinct(H_and_Y))
  dis_dur <- pull(data %>% filter(Subject == part) %>%  distinct(disease_duration))
  Aff_Side <- pull(data %>% filter(Subject == part) %>%  distinct(Aff_Side))
  Subtype <- pull(data %>% filter(Subject == part) %>%  distinct(Subtype))
  LEDD <- pull(data %>% filter(Subject == part) %>%  distinct(LEDD))
  
  
  # save inducer variables
  posterior_cor <- posterior_cor %>% add_row(Subject = part, 
                                             UPDRS = m_score, 
                                             HnY = HnY,
                                             Dis_Dur = dis_dur,
                                             Aff_Side =Aff_Side,
                                             Subtype = Subtype,
                                             LEDD = LEDD,
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
                                             HnY = HnY,
                                             Dis_Dur = dis_dur,
                                             Aff_Side =Aff_Side,
                                             Subtype = Subtype,
                                             LEDD = LEDD,
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
  theme_bw(base_size = 8) +
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


# Plot correlation disease duration and proactive control
sp_disease_dur <- ggscatter(posterior_cor, y = "ProactiveC", x = "Dis_Dur",
                            add = "reg.line", # Add regressin line
                            size = 1,
                            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                            conf.int = TRUE # Add confidence interval,
) + 
  scale_x_continuous(breaks = c(0, 2,4,6,8,10,12)) +
  facet_grid(cols = vars(fct_inorder(Item))) +
  labs(x = "Disease Duration in years", y= "Proactive Control in ms") +
  ggtitle("Disease Duration by Effect Estimates") +
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )  + stat_cor(method = "pearson", label.x = 0, label.y = 90, size = 3)

# Plot correlation Levodopa eqivalence daily dosis and proactive control
sp_LEDD <- ggscatter(posterior_cor, y = "ProactiveC", x = "LEDD",
                     add = "reg.line",  # Add regressin line
                     size = 1,
                     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval,
) + 
  facet_grid(cols = vars(fct_inorder(Item))) +
  labs(x = "LEDD in mg", y= "Proactive Control in ms") +
  ggtitle("LEDD by Effect Estimates") +
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )
# Add correlation coefficient
sp_LEDD <- sp_LEDD + stat_cor(method = "pearson", label.x = 200, label.y = 90, size = 3)

# Next we plot proactive control by affected side
posterior_cor <- posterior_cor %>% 
  mutate(Aff_Side = as.factor(Aff_Side)) %>%
  mutate(new_Aff_Side = recode_factor(Aff_Side, 'l' = "left", 'r'="right", 'b'="both"))

# First select the colors we want to use
SD_cols <- RColorBrewer::brewer.pal(8, "Dark2")[5:8]

sides <- ggplot(data = posterior_cor %>% drop_na(new_Aff_Side), 
                aes(x = fct_inorder(new_Aff_Side), y = ProactiveC, color = fct_inorder(new_Aff_Side)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(Item))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  labs(x = "Affected Side", y = "Proactive Control in ms") +
  ggtitle("Affected Side by Effect Estimates") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Affected Side",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )+
  scale_color_manual(values = SD_cols) +
  scale_fill_manual(values = SD_cols)

# Perform t_test right versus left
left_side_I <- posterior_cor %>%
  filter(Aff_Side == "l",
         Item == "Inducer")
right_side_I <- posterior_cor %>%
  filter(Aff_Side == "r",
         Item == "Inducer")
t.test(left_side_I$ProactiveC, right_side_I$ProactiveC, alternative = "two.sided", var.equal = FALSE)

# and for the diagnostic items
left_side_D <- posterior_cor %>%
  filter(Aff_Side == "l",
         Item != "Inducer")
right_side_D <- posterior_cor %>%
  filter(Aff_Side == "r",
         Item != "Inducer")
t.test(left_side_D$ProactiveC, right_side_D$ProactiveC, alternative = "two.sided", var.equal = FALSE)

# Next we plot proactive control by disease Subtype
posterior_cor <- posterior_cor %>% 
  mutate(Subtype = as.factor(Subtype)) %>%
  mutate(new_Subtype = recode_factor(Subtype, 'trem_d' = "TD", 'akin_reg'="AR", 'eqiv_type'="ET"))
# First select the colors we want to use
DS_cols <- RColorBrewer::brewer.pal(3, "Dark2")[1:3]

subtype <- ggplot(data = posterior_cor %>% drop_na(new_Subtype), 
                  aes(x = fct_inorder(new_Subtype), y = ProactiveC, color = fct_inorder(new_Subtype)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(Item))) +
  geom_boxplot(alpha = 0.5, show.legend =FALSE) +
  labs(x = "Disease Subtype", y = "Proactive Control in ms") +
  ggtitle("Disease Subtype by Effect Estimates") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Disease Subtype",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )+
  scale_color_manual(values = DS_cols) +
  scale_fill_manual(values = DS_cols)

# Perform t_test right versus left
AD_I <- posterior_cor %>%
  filter(new_Subtype == "AD",
         Item == "Inducer")
ET_I <- posterior_cor %>%
  filter(new_Subtype == "ET",
         Item == "Inducer")
t.test(AD_I$ProactiveC, ET_I$ProactiveC, alternative = "two.sided", var.equal = FALSE)

# and for the diagnostic items
AD_D <- posterior_cor %>%
  filter(new_Subtype == "AD",
         Item != "Inducer")
ET_D <- posterior_cor %>%
  filter(new_Subtype == "ET",
         Item != "Inducer")
t.test(AD_D$ProactiveC, ET_D$ProactiveC, alternative = "two.sided", var.equal = FALSE)


# Lastly we plot by H&Y
# First select the colors we want to use
my_cols <- RColorBrewer::brewer.pal(9, "PuBu")[6:9]
#Plot
HNY <- ggplot(data = posterior_cor %>% drop_na(HnY) %>%
                mutate(HnY = as.factor(HnY)), 
              aes(x = HnY, y = ProactiveC, color = HnY))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(Item))) +
  geom_boxplot(alpha = 0.5, show.legend =FALSE) +
  labs(x = "H&Y Stage", y = "Proactive Control in ms") +
  ggtitle("H&Y Stage by Effect Estimates") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "H&Y Stage",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )+
  scale_color_manual(values = my_cols) +
  scale_fill_manual(values = my_cols)



ggplot(data = posterior_cor %>% drop_na(HnY) %>%
         mutate(HnY = as.factor(HnY)), 
       aes(x = HnY,  fill = HnY)) + 
  geom_bar() +
  facet_grid(cols = vars(fct_inorder(Item)))  +
  
  theme_bw(base_size = 14) +
  labs(x = "H&Y", y = "Proactive Control in ms") +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) 
save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"

Clinical_var_cors <- ggarrange(sp_disease_dur, sp_LEDD,
                               labels = c("A", "B"),
                               label.y = 0.97,
                               ncol = 1, nrow = 2, align = "h", widths = c(1,1),
                               #heights = c(4.5, -0.5, 4.5),
                               font.label = list(size = 12, color = "black", face = "bold", family = NULL),
                               common.legend = TRUE, legend = "right")

save_n <- "CorrealtionPlots_ClinicalVars.png"
Clinical_var_cors
ggsave(path = save_path, filename = save_n,  dpi=600,
       width = 8,
       height = 10.5,
       units = c("cm"))

Clinical_var_types <- ggarrange(sides, subtype, HNY ,
                                labels = c("A", "B", "C" ),
                                label.y = 0.97,
                                font.label = list(size = 12, color = "black", face = "bold", family = NULL),
                                ncol = 1, nrow = 3, align = "h", widths = c(1,1))
#heights = c(4.5, -0.5, 4.5),
#common.legend = TRUE, legend = "right")

save_n <- "Association_Plots_ClinicalVars.png"
Clinical_var_types
ggsave(path = save_path, filename = save_n,  dpi=600,
       width = 8,
       height = 12.5,
       units = c("cm"))

# redo with cowplot
left_row <- plot_grid(sp_disease_dur, sp_LEDD,
                      labels = c("A", "B"),
                      ncol = 1)

right_row <- plot_grid(sides, subtype, HNY ,
                       labels = c("C", "D", "E" ),
                       ncol = 1)
# join both
full_plot <- plot_grid(left_row, right_row, labels = c("", ""))

save_n <- "ClinicalVars.png"
full_plot
ggsave(path = save_path, filename = save_n,  dpi=600,
       width = 16,
       height = 16,
       units = c("cm"))

# Get statistics Subtype
posterior_cor %>% group_by(Subtype, Item) %>%
  summarize(part_nr = n())

# Get statistics Side
posterior_cor %>% group_by(Aff_Side, Item) %>%
  summarize(part_nr = n())

# Get statistics Disease Duration
dat <- posterior_cor %>% group_by(Item) %>%
  summarize(mean_dur = mean(Dis_Dur, na.rm = TRUE), SD_Dur = sd(Dis_Dur, na.rm = TRUE),
            mean_LEDD = mean(LEDD), sd_LEDD = sd(LEDD),
            mean_UPDRS = mean(UPDRS,na.rm = TRUE), sd_UPDRS = sd(UPDRS, na.rm = TRUE),
            mean_HnY = mean(HnY), sd_HnY = sd(HnY)) 