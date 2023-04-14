#####
# Author: Julius Kricheldorff
# Exploring time series effects of adaptive control in the LWPCE effect for PD and CO seperately

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
  select(response_time_Keyboard_response, Age, Congruency, Gender, Analysis_type, Exp_group, Subject, Trl_type, correct_Keyboard_response, Block, Item, Item_specific, Congruency_main, UPDRS, total_trial) %>%
  dplyr::rename(RT = response_time_Keyboard_response, Correct_A = correct_Keyboard_response) %>%
  filter(Subject != "sub-PD_c(\"11\", \"2\")" &  Subject != "sub-PD_16") %>%#exclude participant 16 and the partial dataset by participant
  filter(Subject != "sub-CO_0" &  Subject != "sub-CO_3" & Subject != "sub-CO_7" &  Subject != "sub-CO_9") %>% # exclude because they used both hands
  mutate(Error = case_when( # we define accuracy in terms of errors, easier to interpret
    Correct_A == 1 ~ 0,
    Correct_A == 0 ~1
  )) 

Data_one_P <- Data %>%
  filter(Subject == "sub-CO_1",
         Block != "practice",
         Item_specific == "List_wide") %>%
  group_by(Analysis_type, Block, Trl_type, Congruency) %>%
  summarize(N = n())
  

Data_correct <- Data %>%
  filter( RT > 200) %>% # Remove data with small RT, otherwise messes with our non-decision parameter 
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group)) %>% # rank items by group
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency)%>%
  mutate(rank = dense_rank(total_trial)) %>% # get rank of individual trials
  mutate(quant_rank = ntile(rank,4)) %>% # get quantile rank of items
  ungroup() %>%
  mutate(Analysis_type = case_when(
    Analysis_type == "main_con" | Analysis_type == "MC" ~ "MC",
    Analysis_type == "main_incon" | Analysis_type == "MI" ~ "MI"
  )) %>%
  mutate(Congruency = case_when(
    Congruency == "incongruent" ~"Incongruent",
    Congruency == "congruent" ~ "Congruent"
  ))  %>%
  mutate(Exp_group = case_when(
    Exp_group == "PD" ~ "Parkinson",
    Exp_group == "CO" ~ "Control"
  )) %>%
  mutate(Item_specific = case_when(
    Item_specific == "Item_spec" ~"ISPCE",
    Item_specific == "List_wide" ~"LWPCE"
  ))


#### Next we want to calculate the average reaction time by percentile of items:
# Here we go into 

# SUM_RT  <- Data_correct %>% 
#   group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
#   mutate(Av_RT = mean(RT)) %>%
#   ungroup() %>% 
#   group_by(Exp_group, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
#   mutate(grand_RT = mean(RT)) %>% # calculate grand average 
#   ungroup() %>%
#   group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
#   mutate(corrected_RT = RT - Av_RT + grand_RT) %>%# correct RT for within subject variance
#   summarize(mean_RT = mean(RT), mean_cRT = mean(corrected_RT)) %>% 
#   ungroup() %>%
#   # get group variable again
#   mutate(Exp_group = substr(Subject, start = 5, stop = 6)) %>%
#   group_by(Exp_group, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
#   summarize(Group_RT = mean(mean_RT), SE = sd(mean_RT)/sqrt(30))
#   

SUM_RT  <- Data_correct %>% 
  group_by(Subject, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  summarize(Av_RT = mean(RT)) %>%
  ungroup() %>%
  # get group variable again
  mutate(Exp_group = substr(Subject, start = 5, stop = 6)) %>%
  mutate(Exp_group = case_when(
    Exp_group == "PD" ~ "Parkinson",
    Exp_group == "CO" ~ "Control"
  )) %>%
  group_by(Exp_group, Item_specific, Trl_type,  Analysis_type, Congruency, quant_rank) %>%
  summarize(Group_RT = mean(Av_RT), SE = sd(Av_RT)/sqrt(30)) %>% 
  ungroup()
    

## Next let us plot the results for the LWPCE effect first
save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"
LWPCE_SUM <- SUM_RT %>%
  filter(Item_specific == "LWPCE")


ggt_LWPCE <- ggplot(data = LWPCE_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
  facet_grid(Exp_group ~ Analysis_type) +
  scale_x_continuous(labels=c("1" = "Q1", "2" = "Q2", "3" = "Q3", "4" = "Q4"), 
                     minor_breaks = c(1, 2, 3 , 4)) +
  geom_line(aes(linetype = Trl_type), size = 0.5) + 
  geom_point() +
  coord_cartesian(ylim = c(600, 820)) +
  labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item") +
  ggtitle(NULL,subtitle = "RT over Time") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  scale_y_continuous(
    breaks = seq(650, 800, by = 50),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  scale_linetype_manual(values = c(4,1))  +
  guides(color = guide_legend(override.aes = list(shape = 15, 
                                                  size = 6,
                                                  linetype = 0)))
  #geom_errorbar(aes(ymin = Group_RT - SE, ymax = Group_RT + SE), width =.08)


save_n <- "RTbyTime_LWPCE.png"
  
ggsave(path = save_path, filename = save_n,  dpi=300)
  

  
ISPCE_SUM <- SUM_RT %>%
  filter(Item_specific != "LWPCE")
  
ggt_ISPCE <- ggplot(data = ISPCE_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
  facet_grid(Exp_group ~ Analysis_type) +
  scale_x_continuous(labels=c("1" = "Q1", "2" = "Q2", "3" = "Q3", "4" = "Q4"), 
                     minor_breaks = c(1, 2, 3 , 4)) +
  geom_line(aes(linetype = Trl_type), size = 0.5) + 
  geom_point() +
  coord_cartesian(ylim = c(600, 820)) +
  labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item") +
  #ggtitle("RT over Time") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  scale_y_continuous(
    breaks = seq(650, 800, by = 50),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    #strip.text.x = element_text(face = "bold"),
    strip.text.x = element_blank(), # remove facet to combine plots
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  scale_linetype_manual(values = c(4,1)) +
  guides(color = guide_legend(override.aes = list(shape = 15, 
                                                  size = 6,
                                                  linetype = 0)))
#geom_errorbar(aes(ymin = Group_RT - SE, ymax = Group_RT + SE), width =.08)


save_n <- "RTbyTime_ISPCE.png"
  
ggsave(path = save_path, filename = save_n,  dpi=300)
  
# Next we plot the descriptives plot

# first all panels
Data_RT_ISPCE <- Data_correct 

rbp <- ggplot(data = Data_RT_ISPCE, aes(x = Analysis_type, y = RT, fill = Congruency)) +
  facet_grid(Exp_group ~ Item_specific)  + 
  theme_bw(base_size = 12) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 1), side = "left", adjust = 1, width = 0.35, height = 1, justification  = 2, .width = 0, point_colour = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(width = 1)) +
  geom_boxplot(alpha = 0.7, width = 0.22, position = position_dodge(width = 1), outlier.color = NA) +
  labs(y= "RT in ms", x = "Measure", fill = "Item") +
  #ggtitle("RT Distributions") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  scale_y_continuous(
    breaks = seq(500, 1500, by = 300),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
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
  ) +
  scale_fill_manual(values = c("#74C476","#00441B")) +
  coord_cartesian(ylim = c(150, 1550))



# then the effects seperately

Data_RT_ISPCE <- Data_correct %>%
  filter(Item_specific == "ISPCE")

rbp_ISPCE <- ggplot(data = Data_RT_ISPCE, aes(x = Analysis_type, y = RT, fill = Congruency)) +
  facet_grid(Exp_group ~ Item_specific)  + 
  theme_bw(base_size = 12) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 1), side = "left", 
               adjust = 1, width = 0.35, height = 1, justification  = 2, .width = 0, point_colour = NA, show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(width = 1)) +
  geom_boxplot(alpha = 0.7, width = 0.3, position = position_dodge(width = 1), outlier.color = NA, show.legend = FALSE) +
  labs(y= "RT in ms", x = "Measure", fill = "Item") +
  ggtitle(NULL, subtitle = "ISPCE") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  scale_y_continuous(
    breaks = seq(500, 1500, by = 300),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(face = "bold",
                                 size = 13),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none",
    #legend.text = element_text(face = "bold"),
    #legend.title = element_text(face = "bold"),
    #strip.text.x = element_text(hjust = 0,
      #face="bold"),
    strip.text.x = element_blank(), # remove facet to combine
    strip.text.y = element_blank(),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_fill_manual(values = c("#74C476","#00441B")) +
  coord_cartesian(ylim = c(150, 1550))

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
Con_RT <- Data_correct %>% filter(Congruency == "Congruent") %>% 
  group_by(Exp_group, Item_specific, Analysis_type) %>%
  summarise(medianRT = median(RT))

Incon_RT <- Data_correct %>% filter(Congruency == "Incongruent") %>% 
  group_by(Exp_group, Item_specific, Analysis_type) %>%
  summarise(medianRT = median(RT))

Dif_RT <- Con_RT 
Dif_RT$RT_Dif <- round(Incon_RT$medianRT - Con_RT$medianRT)

anno_df <- data.frame(
  PANEL = rep(seq(1,4,by = 1), times = 1, each = 2),
  start = rep(x_pos1, times = 4, each = 1), # set the starting point of the comparison
  end   = rep(x_pos2, times = 4, each = 1), # set the starting point of the comparison
  y_pos = 1300, # on the y-axis
  label = as.factor(paste(Dif_RT$RT_Dif, ' ms')),
  Congruency = c("Congruent", "Incongruent"), 
  Analysis_type = c("MC", "MI"),
  Exp_group = Dif_RT$Exp_group,
  Item_specific = Dif_RT$Item_specific
)

anno_df_ISPCE <- as_tibble(anno_df) %>%
  filter(Item_specific == "ISPCE") %>%
  mutate(PANEL = case_when(
    PANEL == 3 ~ 2, 
    TRUE ~ PANEL
  ))

## Now let us plot again

rbp_ISPCE <- rbp_ISPCE + ggsignif::geom_signif(
  data = anno_df_ISPCE,
  aes(xmin = start, xmax = end, annotations = label, y_position = y_pos),
  manual = TRUE,
  textsize = 3.5, vjust = -0.1, tip_length = 0.02 
)
save_n <- "Descriptive_ISPCE_RT.png"

ggsave(path = save_path, filename = save_n,  dpi=300)

# Now let us plot the LWPCE data

Data_RT_LWPCE <- Data_correct %>%
  filter(Item_specific != "ISPCE")

rbp_LWPCE <- ggplot(data = Data_RT_LWPCE, aes(x = Analysis_type, y = RT, fill = Congruency)) +
  facet_grid(Exp_group ~ Item_specific)  + 
  theme_bw(base_size = 12) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 1), side = "left", adjust = 1, width = 0.35, 
               height = 1, justification  = 2, .width = 0, point_colour = NA, show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(width = 1)) +
  geom_boxplot(alpha = 0.7, width = 0.3, position = position_dodge(width = 1), outlier.color = NA, show.legend = FALSE) +
  labs(y= "RT in ms", x = "Measure", fill = "Item") +
  ggtitle("LWPCE", subtitle = "RT") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  #scale_x_discrete(position = "top") +
  scale_y_continuous(
    breaks = seq(500, 1500, by = 300),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    #legend.text = element_text(face = "bold"),
    #legend.title = element_text(face = "bold"),
    #strip.text.x = element_text(face="bold"),
    strip.text.x = element_blank(), # remove facet to combine
    strip.text.y = element_blank(),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.title.position  =  "panel"
  ) +
  scale_fill_manual(values = c("#74C476","#00441B")) +
  coord_cartesian(ylim = c(150, 1550))

anno_df_LWPCE <- as_tibble(anno_df) %>%
  filter(Item_specific != "ISPCE") %>%
  mutate(PANEL = case_when(
    PANEL == 3 ~ 2, 
    TRUE ~ PANEL
  ))

## Now let us plot again

rbp_LWPCE <- rbp_LWPCE + ggsignif::geom_signif(
  data = anno_df_LWPCE,
  aes(xmin = start, xmax = end, annotations = label, y_position = y_pos),
  manual = TRUE,
  textsize = 3.5, vjust = -0.1, tip_length = + 0.02 
)
save_n <- "Descriptive_LWPCE_RT.png"

ggsave(path = save_path, filename = save_n,  dpi=300)

# Now let us combine the plots

LWPCE_comb <- ggarrange(rbp_LWPCE + rremove("legend"), ggt_LWPCE,
          labels = c("A", "B"),
          ncol = 2, nrow = 1, align = "h", widths = c(1,2))
annotate_figure(LWPCE_comb,
                top =text_grob( "List-wide Proportion Congruency Effect", face = "bold", size = 20))



save_n <- "Descriptive_LWPCE_combined.png"

ggsave(path = save_path, filename = save_n,  dpi=300)



# Now let us combine the plots

ISPCE_comb <- ggarrange(rbp_ISPCE + rremove("legend"), ggt_ISPCE,
          labels = c("A", "B"),
          ncol = 2, nrow = 1, align = "h", widths = c(1,2))

annotate_figure(ISPCE_comb,
                top =text_grob( "Item-specific Proportion Congruency Effect", face = "bold", size = 20))



save_n <- "Descriptive_ISPCE_combined.png"

ggsave(path = save_path, filename = save_n,  dpi=300)

# Do it for all the data
Descriptives <- ggarrange(rbp_LWPCE + rremove("legend") + rremove("x.text") + rremove("xlab"), 
                          ggt_LWPCE + rremove("x.text") + rremove("xlab"), 
                          NULL, NULL,
                          rbp_ISPCE + rremove("legend"), ggt_ISPCE,
          labels = c("A", "B", "C", "D", NULL, NULL),
          label.y = 0.97,
          vjust = 0.2,
          ncol = 2, nrow = 3, align = "h", widths = c(3,4),
          heights = c(4.5, -0.5, 4.5),
          common.legend = TRUE, legend = "right")

annotate_figure(Descriptives,
                top = text_grob( "Reaction Time Distributions", face = "bold", size = 15))


save_n <- "Descriptives.png"

ggsave(path = save_path, filename = save_n,  dpi=300,
       width = 16.5,
       height = 19,
       units = c("cm"))

# To be done:
#   Move MC MI from the RT bottom to the top
# Show what is LWPCE and what is ISPCE 

# Poster Plot


ggt_LWPCE2 <- ggplot(data = LWPCE_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
  facet_grid(Exp_group ~ Analysis_type) +
  scale_x_continuous(labels=c("1" = "Q1", "2" = "Q2", "3" = "Q3", "4" = "Q4"), 
                     minor_breaks = c(1, 2, 3 , 4)) +
  geom_line(aes(linetype = Trl_type), size = 0.5) + 
  geom_point(size = 0.5) +
  coord_cartesian(ylim = c(600, 820)) +
  labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item") +
  ggtitle(NULL,subtitle = "RT over Time") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8) +
  scale_y_continuous(
    breaks = seq(650, 800, by = 50),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  scale_linetype_manual(values = c(4,1))  +
  guides(color = guide_legend(override.aes = list(shape = 15, 
                                                  size = 4,
                                                  linetype = 0)))


ggt_ISPCE2 <- ggplot(data = ISPCE_SUM, aes(x= quant_rank , y = Group_RT, color = Congruency, linetype = Trl_type)) + 
  facet_grid(Exp_group ~ Analysis_type) +
  scale_x_continuous(labels=c("1" = "Q1", "2" = "Q2", "3" = "Q3", "4" = "Q4"), 
                     minor_breaks = c(1, 2, 3 , 4)) +
  geom_line(aes(linetype = Trl_type), size = 0.5) + 
  geom_point(size = 0.5) +
  coord_cartesian(ylim = c(600, 820)) +
  labs(y= "RT in ms", x = "Quartile", Color = "Congruency", linetype = "Item") +
  #ggtitle("RT over Time") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8) +
  scale_y_continuous(
    breaks = seq(650, 800, by = 50),
    expand = c(.001, .001)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    #legend.key = element_rect(fill = "White", colour = "black"),
    #strip.text.x = element_text(face = "bold"),
    strip.text.x = element_blank(), # remove facet to combine plots
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_color_manual(values = c("#74C476","#00441B")) +
  scale_linetype_manual(values = c(4,1)) +
  guides(color = guide_legend(override.aes = list(shape = 15, 
                                                  size = 4,
                                                  linetype = 0)))
#geom_errorbar(aes(ymin = Group_RT - SE, ymax = Group_RT + SE), width =.08)


Descriptives2 <- ggarrange(ggt_LWPCE2 + labs(subtitle = "LWPCE") + theme(plot.subtitle = element_text(face = "bold")) + rremove("x.text") + rremove("xlab"), 
                           ggt_ISPCE2 + labs(subtitle = "ISPCE") + theme(plot.subtitle = element_text(face = "bold")),
                          labels = c("A", "B"),
                          label.y = 0.97,
                          vjust = 0.2,
                          ncol = 1, nrow = 2, align = "h", widths = c(1,1),
                          #heights = c(4.5, -0.5, 4.5),
                          common.legend = TRUE, legend = "right")

annotate_figure(Descriptives2,
                top = text_grob( "RT over time", face = "bold", size = 10))


save_n <- "DescriptivesPoster.png"

ggsave(path = save_path, filename = save_n,  dpi=300,
       width = 8,
       height = 10.5,
       units = c("cm"))
