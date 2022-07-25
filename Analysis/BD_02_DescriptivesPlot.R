# Plot Descriptives Reaction Time

# Analysis on the behavioral data of healthy control participants and participants with Parkinson's disease
setwd('C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data')

# Load packages
library(brms)
library(ggplot2)
library(dplyr)
library(haven)
library(tidyr)
library(ggdist)
library(ggpattern)

# Set a seed for sake of reproducibility
set.seed(32923)

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
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>%  # only correct responses
  mutate(Exp_group = if_else(Exp_group == "PD", "Parkinson Participants", "Control Participants", missing = NULL)) %>%
  mutate(Effect = case_when(
    Analysis_type == "main_con" | Analysis_type == "main_incon"  ~ "Reactive Control",
    Analysis_type == "MC" | Analysis_type == "MI"  ~ "Proactive Control")) %>%
  mutate(Measure = case_when(
    Analysis_type == "main_con" & Congruency == "congruent"  ~ "MC",
    Analysis_type == "main_con" & Congruency == "incongruent"  ~ "MC",
    Analysis_type == "main_incon" & Congruency == "congruent"  ~ "MI",
    Analysis_type == "main_incon" & Congruency == "incongruent"  ~ "MI",
    Analysis_type == "MC" & Congruency == "congruent"  ~ "MC",
    Analysis_type == "MC" & Congruency == "incongruent"  ~ "MC",
    Analysis_type == "MI" & Congruency == "congruent"  ~ "MI",
    Analysis_type == "MI" & Congruency == "incongruent"  ~ "MI"
  )) #%>%
  #dplyr::filter(Trl_type == "Diagnostic") # only correct responses
    
      

# Plot Reaction Times

rbp <- ggplot(data = Data, aes(x = Measure, y = RT, fill = Congruency)) +
  facet_grid(Exp_group ~ Effect)  + 
  theme_bw(base_size = 14) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 1), side = "left", adjust = 1, width = 0.35, height = 1, justification  = 2, .width = 0, point_colour = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(width = 1)) +
  geom_boxplot(alpha = 0.7, width = 0.3, position = position_dodge(width = 1), outlier.color = NA) +
  coord_cartesian(ylim = c(150, 1250)) +
  labs(y= "RT in ms", x = "Measure", fill = "Item Type") +
  ggtitle("Reaction Time Distributions") + 
  theme(
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel",
    panel.background = element_blank()
  ) +
  scale_fill_manual(values = c("#74C476","#00441B")) 

## we want to annotate the mean difference between conditions - to know how - we have to deconstruct the plot:
my_plot <- ggplot_build(rbp)

## create dataframe with parameter of interest
annotation_df <- data.frame(
  x = my_plot$data[[2]][["x"]],
  PANEL = my_plot$data[[2]][["PANEL"]],
  group =  my_plot$data[[2]][["group"]]
)

# create better frame again
x_pos1 = c(annotation_df$x[1], annotation_df$x[3])
x_pos2 = c(annotation_df$x[2], annotation_df$x[4])
Con_RT <- Data %>% filter(Congruency == "congruent") %>% 
  group_by(Exp_group, Effect, Measure) %>%
  summarise(medianRT = median(RT))

Incon_RT <- Data %>% filter(Congruency == "incongruent") %>% 
  group_by(Exp_group, Effect, Measure) %>%
  summarise(medianRT = median(RT))

Dif_RT <- Con_RT 
Dif_RT$RT_Dif <- round(Incon_RT$medianRT - Con_RT$medianRT)

anno_df <- data.frame(
  PANEL = rep(seq(1,4,by = 1), times = 1, each = 2),
  start = rep(x_pos1, times = 4, each = 1), # set the starting point of the comparison
  end   = rep(x_pos2, times = 4, each = 1), # set the starting point of the comparison
  y_pos = 200, # on the y-axis
  label = as.factor(paste(Dif_RT$RT_Dif, ' ms')),
  Congruency = c("congruent", "incongruent"), 
  Measure = c("MC", "MI"),
  Exp_group = Dif_RT$Exp_group,
  Effect = Dif_RT$Effect
)

## Now let us plot again

rbp + ggsignif::geom_signif(
  data = anno_df,
  aes(xmin = start, xmax = end, annotations = label, y_position = y_pos),
  manual = TRUE,
  textsize = 3.5, vjust = 1.8, tip_length = - 0.02 
)