#####
# Author: Julius Kricheldorff
# Plotting average Theta band activity over FCz FC1 and FC2 electrodes

####
setwd('E:/AdaptiveControl/Data/BehaviorResults')

# Load packages
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstudioapi)
library(xtable)
library(ggplot2)
library(ggpubr)
library(ggdist)
library(ggpattern)
library(plotrix) # to calculate standard error

# Set a seed for sake of reproducibility
set.seed(32946)

# save location
save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/EEG/"

# load the data

path1<- "E:/AdaptiveControl/Data/TF_for_R"
files = list.files(path = path1, pattern = "Theta_data", full.names = TRUE, recursive = FALSE)
df = tibble()

for (x in files) {
  t <- read.csv(x, header = TRUE)
  t <- as_tibble(t)
  df <- bind_rows(df, t)
}

# We need to add congruency as a factor
Data1 <- df %>%
  mutate(Congruency = case_when(
    substr(Condition, start = 4, stop = 4) == "I" ~ "Incongruent",
    substr(Condition, start = 4, stop = 4) == "C" ~ "Congruent")) %>%
  mutate(Prop_type = substr(Condition, start = 1, stop = 2)) %>%
  mutate(Group = case_when(
    Group == "CO" ~ "HC",
    TRUE ~ as.character(Group)
  ))#rename CO to HC for consistency


# We also need to include the period of significant difference for the cluster
# For that I need to create a dataset by hand with the respective intervals
Data <- Data1 %>%
  #mutate(sig = NA) %>%
  mutate(sig = case_when(
    Group == "HC" & Lock == "SL" & Ana_type == "ISPC" & Time >= 0.392 & Time <= 0.752 ~ -1, 
    Group == "HC" & Lock == "RL" & Ana_type == "ISPC" & Time >= -0.44 & Time <= -0.104 ~ -1,
    Group == "HC" & Lock == "SL" & Ana_type == "LWPC" & Time >= 0.352 & Time <= 0.768 ~ -1,
    Group == "HC" & Lock == "RL" & Ana_type == "LWPC" & Time >= -0.6 & Time <= -0.048 ~ -1,
    Group == "PD" & Lock == "SL" & Ana_type == "LWPC" & Time >= 0.328 & Time <= 0.712 ~ -1,
    Group == "PD" & Lock == "RL" & Ana_type == "LWPC" & Time >= -0.6 & Time <= -0.144 ~ -1
  ))


# Next Lets Plot the data Stimulus_locked
SL_ISPC <- Data %>%
  filter(Ana_type == "ISPC",
         Lock == "SL")
hjust_var = 1


Plt_sl_ISPCE <- ggplot(data = SL_ISPC, aes(x= Time , y = Pw, color = Prop_type, linetype = Congruency)) +
  facet_grid(. ~ Group ) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin= Pw-SE, ymax= Pw + SE, fill = Condition), alpha = 0.1, linetype = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#E16A86","#00A6CA")) +
  scale_fill_manual(values = c("#E16A86","#E16A86", "#00A6CA","#00A6CA")) +
  labs(y= "Power in dB", x = "Time in s", color = "Congruency Proportion", linetype = "Item Type") +
  #ggtitle("4-8 Hz Power Frontal Electrodes ISPCE") + 
  scale_x_continuous(breaks = c(-0.25, 0, 0.5, 1 ),
                     labels = c(-0.25, 0,  0.5, 1)) +
  coord_cartesian(ylim = c( -3.6, -2)) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold", angle = 45, hjust = hjust_var),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  geom_line(aes(x = Time, y = sig-2.5),na.rm = TRUE,size = 1.5, color = "black", show.legend =FALSE) +
  guides(fill = "none") +
  guides(color = guide_legend(override.aes = list(
                                                  size = 6,
                                                  linetype = c(0, 0),
                                                  fill = c("#E16A86","#00A6CA"),
                                                  alpha = 1)))

save_n <- "Average_theta_ISPC_SL_new.png"
ggsave(path = save_path, filename = save_n,  dpi=600)




RL_ISPC <- Data %>%
  filter(Ana_type == "ISPC",
         Lock != "SL")


Plt_rl_ISPCE <- ggplot(data = RL_ISPC, aes(x= Time , y = Pw, color = Prop_type, linetype = Congruency)) +
  facet_grid(. ~ Group ) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin= Pw-SE, ymax= Pw + SE, fill = Condition), alpha = 0.1, linetype = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#E16A86","#00A6CA")) +
  scale_fill_manual(values = c("#E16A86","#E16A86", "#00A6CA","#00A6CA")) +
  labs(y= "Power in dB", x = "Time in s", color = "Congruency Proportion", linetype = "Item Type") +
  #ggtitle("4-8 Hz Power Frontal Electrodes ISPCE") + 
  theme_bw(base_size = 9) +
  scale_x_continuous(breaks = c(-1, -0.5,  0, 0.25),
                     labels = c(-1, -0.5,  0, 0.25)) +
  coord_cartesian(ylim = c( -3.6, -2)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold", angle = 45, hjust = hjust_var),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  guides(fill = "none")  +
  geom_line(aes(x = Time, y = sig-2.5),na.rm = TRUE,size = 1.5, color = "black", show.legend =FALSE) +
  guides(color = guide_legend(override.aes = list(
    size = 6,
    linetype = c(0, 0),
    fill = c("#E16A86","#00A6CA"),
    alpha = 1)))


save_n <- "Average_theta_ISPC_RL_new.png"
ggsave(path = save_path, filename = save_n,  dpi=600)


SL_LWPC <- Data %>%
  filter(Ana_type == "LWPC",
         Lock == "SL")

Plt_sl_LWPCE <- ggplot(data = SL_LWPC, aes(x= Time , y = Pw, color = Prop_type, linetype = Congruency)) +
  facet_grid(. ~ Group ) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin= Pw-SE, ymax= Pw + SE, fill = Condition), alpha = 0.1, linetype = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#E16A86","#00A6CA")) +
  scale_fill_manual(values = c("#E16A86","#E16A86", "#00A6CA","#00A6CA")) +
  labs(y= "Power in dB", x = "Time in s", color = "Congruency Proportion", linetype = "Item Type") +
  #ggtitle("4-8 Hz Power Frontal Electrodes LWPCE") + 
  scale_x_continuous(breaks = c(-0.25, 0, 0.5, 1 ),
    labels = c(-0.25, 0,  0.5, 1)) + 
  coord_cartesian(ylim = c( -4.2, -1.8)) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold", angle = 45, hjust = hjust_var),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  guides(fill = "none")  +
  geom_line(aes(x = Time, y = sig-3.1),na.rm = TRUE,size = 1.5, color = "black", show.legend =FALSE) +
  guides(color = guide_legend(override.aes = list(
    size = 6,
    linetype = c(0, 0),
    fill = c("#E16A86","#00A6CA"),
    alpha = 1)))

save_n <- "Average_theta_LWPC_SL_new.png"
ggsave(path = save_path, filename = save_n,  dpi=600)



RL_LWPC <- Data %>%
  filter(Ana_type == "LWPC",
         Lock != "SL")


Plt_rl_LWPCE <- ggplot(data = RL_LWPC, aes(x= Time , y = Pw, color = Prop_type, linetype = Congruency)) +
  facet_grid(. ~ Group ) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin= Pw-SE, ymax= Pw + SE, fill = Condition), alpha = 0.1, linetype = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#E16A86","#00A6CA")) +
  scale_fill_manual(values = c("#E16A86","#E16A86", "#00A6CA","#00A6CA")) +
  labs(y= "Power in dB", x = "Time in s", color = "Congruency Proportion", linetype = "Item Type") +
  #ggtitle("4-8 Hz Power Frontal Electrodes LWPCE") + 
  scale_x_continuous(breaks = c(-1, -0.5,  0, 0.25),
                     labels = c(-1, -0.5,  0, 0.25))  +
  coord_cartesian(ylim = c( -4.2, -1.8)) +
  theme_bw(base_size = 9) +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold", angle = 45, hjust = hjust_var),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  guides(fill = "none")  +
  geom_line(aes(x = Time, y = sig-3.1),na.rm = TRUE,size = 1.5, color = "black", show.legend =FALSE) +
  guides(color = guide_legend(override.aes = list(
    size = 6,
    linetype = c(0, 0),
    fill = c("#E16A86","#00A6CA"),
    alpha = 1)))


save_n <- "Average_theta_LWPC_RL_new.png"
ggsave(path = save_path, filename = save_n,  dpi=600)


# Now let us combine the plots for the ISPC data
LWPCE_comb <- ggarrange(Plt_sl_LWPCE, Plt_rl_LWPCE,
                        ncol = 2, nrow = 1, align = "h", common.legend = TRUE, legend = "bottom", widths = c(1,1))

save_n <- "Average_theta_LWPC.png"
ggsave(
  path = save_path,
  LWPCE_comb,
  filename = save_n,
  #device = agg_png,
  width = 16,
  height = 5.5,
  units = "cm", dpi = 600
)

ISPCE_comb <- ggarrange(Plt_sl_ISPCE, Plt_rl_ISPCE,
                        ncol = 2, nrow = 1, align = "h", common.legend = TRUE, legend = "bottom",widths = c(1,1))

save_n <- "Average_theta_ISPC.png"
ggsave(
  path = save_path,
  ISPCE_comb,
  filename = save_n,
  #device = agg_png,
  width = 16, 
  height = 5.5,
  units = "cm", dpi = 600
)


#### For the Next part we want to compare model shifted log-normal model outputs with our Theta data
# For this analysis we collapsed theta data across the cluster electrodes, across the time
# where the cluster differed significantly and for each subject seperately. We then calculated 
# the simple conditional means for each subjects and subtracted MC - congruent from MI -congruent 
# and did the same for the incongruent items. The idea was to investigate whether particular
# context effects drive particular electrophysiological context adaptions. The previouse
# plots may indicate such a relationship.

# First we need to extract the subject specific effects of the diagnostic items for the LWPCE effect
# which we created in the BD_06_RandomSlope_LWPC_Dia_Results.R file

path_behav_sum <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data"
behav_n <-  "Participant_specific_effect_summary.csv"

Behav_sum <- read.csv(file.path(path_behav_sum,behav_n), header = TRUE)
Behav_sum <- as_tibble(Behav_sum)
# Next let us load the theta data
path1<- "E:/AdaptiveControl/Data/TF_for_R"
files = list.files(path = path1, pattern = "Theta_cluster_averaged", full.names = TRUE, recursive = FALSE)
df2 = tibble()

for (x in files) {
  tca <- read.csv(x, header = TRUE)
  tca <- as_tibble(tca)
  df2 <- bind_rows(df2, tca)
}

# next merge the Adaptive control and theta Info
Theta_X_Behav <- df2 %>%
  filter(Ana_type == "LWPC") %>% # only look at the "LWPC effects here"
  mutate(Control = NA) # create an empty variable where we will save the effects of proactive control 

# Now merge
for (ind in seq(1,nrow(Behav_sum))){
  Theta_X_Behav$Control[(Theta_X_Behav$Part == Behav_sum$Subject[ind])] = Behav_sum$ProactiveC[ind]
  print(Behav_sum$Subject[ind])
  
}

# Let us first look at the stimulus locked incongruent trials
LWPC_SL_INC <- Theta_X_Behav %>%
  filter(Lock == "SL",
         Condition == "Incongruent") 

ggplot(data = LWPC_SL_INC, aes(x = Pw, y = Control, color = Group)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  labs(x= "Power in dB", x = "Effect in ms", Color = "Group") +
  ggtitle("Incongruent difference SL") 

# Let us next look at the stimulus locked congruent trials
LWPC_SL_C <- Theta_X_Behav %>%
  filter(Lock == "SL",
         Condition == "Congruent") 

ggplot(data = LWPC_SL_C, aes(x = Pw, y = Control, color = Group)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  labs(x= "Power in dB", x = "Effect in ms", Color = "Group") +
  ggtitle("Congruent difference SL") 
  

# Let us next look at the respones locked incongruent trials
LWPC_RL_INC <- Theta_X_Behav %>%
  filter(Lock == "RL",
         Condition == "Incongruent") 

ggplot(data = LWPC_RL_INC, aes(x = Pw, y = Control, color = Group)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  labs(x= "Power in dB", x = "Effect in ms", Color = "Group") +
  ggtitle("Incongruent difference RL") 
  

# Let us next look at the response locked congruent trials
LWPC_RL_C <- Theta_X_Behav %>%
  filter(Lock == "RL",
         Condition == "Congruent") 

ggplot(data = LWPC_RL_C, aes(x = Pw, y = Control, color = Group)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) +
  labs(x= "Power in dB", x = "Effect in ms", Color = "Group") +
  ggtitle("Congruent difference RL") 
  