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
library(forcats)
library(ggdist)
library(ggtext)
library(colorspace)
library(ragg)
library(ggpp) # so we can use the function dodgennudge
library(ggpubr)

# Set a seed for sake of reproducibility
set.seed(32946)


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

# calculate effect estimates in ms
conditional_effect_calc_shift <- function(effect, model, Item){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # Initialize matrix were we save our data
  model_effects <- tibble(Mean_RT = numeric(),
                          Effect_name = character(),
                          Item_type = character(),
                          Group = character(),
                          Effect = character(),
                          Plot_type = character(),
  )
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
    
    Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2
    Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2
    Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2
    Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2
    Stroop_MC_PD = MC_I_PD - MC_C_PD
    Stroop_MI_PD = MI_I_PD - MI_C_PD
    Stroop_MC_CO = MC_I_CO - MC_C_CO
    Stroop_MI_CO = MI_I_CO - MI_C_CO
    Control_PD = Stroop_MC_PD - Stroop_MI_PD
    Control_CO = Stroop_MC_CO - Stroop_MI_CO
    #Conflict_imp_CO = MC_I_CO - MI_I_CO
    #Facilitation_imp_CO = MI_C_CO - MC_C_CO
    #Conflict_imp_PD = MC_I_PD - MI_I_PD
    #Facilitation_imp_PD = MI_C_PD - MC_C_PD
    #CC_PD = MC_C_PD - MI_C_PD,
    #CI_CO = MI_I_CO - MC_I_CO,
    #CI_PD = MI_I_PD - MC_I_PD)
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MC_PD, eff_n = "Stroop_MC_PD", IT = Item,
                                 Grp = "Parkinson", ef = "LWPCE", Plt = "MC")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MI_PD, eff_n = "Stroop_MI_PD", IT = Item,
                                 Grp = "Parkinson", ef = "LWPCE", Plt = "MI")
    model_effects <- join_postib(main = model_effects, eff_RT = Control_PD, eff_n = "Control_PD", IT = Item,
                                 Grp = "Parkinson", ef = "LWPCE", Plt = "Difference")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Conflict_imp_PD", IT = Item,
    #                              Grp = "Parkinson", ef = "LWPCE", Plt = "Conflict_imp")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Facilitation_imp_PD", IT = Item,
    #                              Grp = "Parkinson", ef = "LWPCE", Plt = "Facilitation_imp")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MC_CO, eff_n = "Stroop_MC_CO", IT = Item,
                                 Grp = "Control", ef = "LWPCE", Plt = "MC")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MI_CO, eff_n = "Stroop_MC_CO", IT = Item,
                                 Grp = "Control", ef = "LWPCE", Plt = "MI")
    model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Control_CO", IT = Item,
                                 Grp = "Control", ef = "LWPCE", Plt = "Difference")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Conflict_imp_CO", IT = Item,
    #                              Grp = "Control", ef = "LWPCE", Plt = "Conflict_imp")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Facilitation_imp_CO", IT = Item,
    #                              Grp = "Control", ef = "LWPCE", Plt = "Facilitation_imp")
    #CC_CO = MC_C_CO - MI_C_CO,
    #CC_PD = MC_C_PD - MI_C_PD,
    #CI_CO = MI_I_CO - MC_I_CO,
    #CI_PD = MI_I_PD - MC_I_PD)
  } else if (effect == "IS"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2
    Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2
    Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2
    Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2
    Stroop_MC_PD = MC_I_PD - MC_C_PD
    Stroop_MI_PD = MI_I_PD - MI_C_PD
    Stroop_MC_CO = MC_I_CO - MC_C_CO
    Stroop_MI_CO = MI_I_CO - MI_C_CO
    Control_PD = Stroop_MC_PD - Stroop_MI_PD
    Control_CO = Stroop_MC_CO - Stroop_MI_CO
    #Conflict_imp_CO = MC_I_CO - MI_I_CO
    #Facilitation_imp_CO = MI_C_CO - MC_C_CO
    #Conflict_imp_PD = MC_I_PD - MI_I_PD
    #Facilitation_imp_PD = MI_C_PD - MC_C_PD
    
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MC_PD, eff_n = "Stroop_MC_PD", IT = Item,
                                 Grp = "Parkinson", ef = "ISPCE", Plt = "MC")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MI_PD, eff_n = "Stroop_MI_PD", IT = Item,
                                 Grp = "Parkinson", ef = "ISPCE", Plt = "MI")
    model_effects <- join_postib(main = model_effects, eff_RT = Control_PD, eff_n = "Control_PD", IT = Item,
                                 Grp = "Parkinson", ef = "ISPCE", Plt = "Difference")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Conflict_imp_PD", IT = Item,
    #                              Grp = "Parkinson", ef = "ISCE", Plt = "Conflict_imp")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Facilitation_imp_PD", IT = Item,
    #                              Grp = "Parkinson", ef = "ISPCE", Plt = "Facilitation_imp")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MC_CO, eff_n = "Stroop_MC_CO", IT = Item,
                                 Grp = "Control", ef = "ISPCE", Plt = "MC")
    model_effects <- join_postib(main = model_effects, eff_RT = Stroop_MI_CO, eff_n = "Stroop_MC_CO", IT = Item,
                                 Grp = "Control", ef = "ISPCE", Plt = "MI")
    model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Control_CO", IT = Item,
                                 Grp = "Control", ef = "ISPCE", Plt = "Difference")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Conflict_imp_CO", IT = Item,
    #                              Grp = "Control", ef = "ISPCE", Plt = "Conflict_imp")
    # model_effects <- join_postib(main = model_effects, eff_RT = Control_CO, eff_n = "Facilitation_imp_CO", IT = Item,
    #                              Grp = "Control", ef = "ISPCE", Plt = "Facilitation_imp")
    #CC_CO = MC_C_CO - MI_C_CO,
    #CC_CO = MC_C_CO - MI_C_CO,
    #CC_PD = MC_C_PD - MI_C_PD,
    #CI_CO = MI_I_CO - MC_I_CO,
    #CI_PD = MI_I_PD - MC_I_PD)
  }
  return(model_effects)
}

join_postib <- function(main, eff_RT, eff_n, IT, Grp, ef, Plt, Blck){
  temp <- tibble(Mean_RT = eff_RT,
                 Effect_name = eff_n,
                 Item_type = IT,
                 Group = Grp,
                 Effect = ef,
                 Plot_type = Plt
  )
  main <- bind_rows(main, temp)
  return(main)
  
}

mod = c("RT")
item_type = c("inducer","diagnostic")
effect = c("LW", "IS")

Conditional_posterior <- tibble(Mean_RT = numeric(),
                                Effect_name = character(),
                                Item_type = character(),
                                Group = character(),
                                Effect = character(),
                                Plot_type = character(),
)

for (eff in effect){
  for (itm in item_type){
    #load model
    model <- load_full_mod(loc = Full_models_saveloc,
                           model_t = mod,
                           item_type = itm,
                           effect = eff)
    
    # calculate posterior conditional effect samples
    eff_post <- conditional_effect_calc_shift(effect = eff,
                                                model = model,
                                              Item = itm)
    # join datasets 
    Conditional_posterior <- bind_rows( Conditional_posterior, eff_post)
    
    #
  }
}

# Now to plot
Conditional_posterior$Plot_type <- factor(Conditional_posterior$Plot_type,levels = c("Difference", "MI", "MC"))
Conditional_posterior$Effect <- factor(Conditional_posterior$Effect,levels = c("LWPCE", "ISPCE"))


Control_eff <- ggplot(data = Conditional_posterior, aes(y = Plot_type, x = Mean_RT, fill = Item_type)) +
  facet_grid(Group ~ Effect)  + 
  theme_bw(base_size = 12) +
  stat_halfeye(aes(fill = Item_type, fill_ramp = stat(level)), .width = c(.66, .95), 
               position = position_dodge(width = 0.8), side = "right", size = 2) +
  scale_fill_ramp_discrete(na.translate = FALSE, range = c(.5, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = round(..x..)),
        #color = Plot_type,
        #color = after_scale(darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    color = "white",
    size =3,
    vjust = -0.6,
    position = position_dodge(width = 0.8)
  ) +
  coord_cartesian(xlim = c(-10, 120)) +
  labs(x= "RT in ms", y = "Congruency Effect", fill = "Item Type", fill_ramp = "Interval") +
  ggtitle("Estimated Marginal Mean Effects") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  theme(
    panel.spacing = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face="bold"),
    strip.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  ) +
  scale_fill_manual(values = c("#B2182B","#2166AC"), labels = c("Diagnostic", "Inducer")) +
  scale_color_manual(values = c("#B2182B","#2166AC"), labels = c("Diagnostic", "Inducer")) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                  size = 6,
                                                  linetype = 0)),
         fill_ramp= guide_legend(override.aes = list(shape = NA, 
                                                    linetype = 0)))

save_path <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Figures/Behavior/"

save_n <- "Posterior_Summary_new1.png"

ggsave(path = save_path, filename = save_n,  dpi=600,
       width = 16.5,
       units = c("cm"))

#Poster plot
new_dat <- Conditional_posterior %>%
  mutate(Plot_type = recode_factor(Plot_type, "Difference" = "Diff",
                                   "MC" = "MC",
                                   "MI" = "MI")) 

Control_eff2 <- ggplot(data =new_dat, aes(y = Plot_type, x = Mean_RT, fill = Item_type)) +
  facet_grid(Group ~ Effect)  + 
  stat_halfeye(aes(fill = Item_type, fill_ramp = stat(level)), .width = c(.66, .95), 
               position = position_dodge(width = 0.78), side = "right", scale = 1, size = 1, height = 0.8) +
  scale_fill_ramp_discrete(na.translate = FALSE, range = c(.5, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = round(..x..)),
    #color = Plot_type,
    #color = after_scale(darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    color = "white",
    size =3,
    vjust = -0.6,
    position = position_dodge(width = 0.78)
  ) +
  coord_cartesian(xlim = c(-10, 120)) +
  labs(x= "RT in ms", y = "Congruency Effect", fill = "Item Type", fill_ramp = "Interval") +
  ggtitle("Estimated Marginal Mean Effects") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  theme(
    panel.spacing = unit(2, "lines"),
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face="bold"),
    strip.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel",
    legend.position="bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.direction = "vertical"
  ) +
  scale_fill_manual(values = c("#B2182B","#2166AC"), labels = c("Diagnostic", "Inducer")) +
  scale_color_manual(values = c("#B2182B","#2166AC"), labels = c("Diagnostic", "Inducer")) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 linetype = 0)),
         fill_ramp= guide_legend(override.aes = list(shape = NA, 
                                                     linetype = 0)))

Descriptives2 <- ggarrange(Control_eff2,
                           #labels = c("A"),
                           ncol = 1, nrow = 1, align = "h", widths = c(2)) +  bgcolor("white")



save_n <- "Posterior_Summary_poster.tiff"

ggsave(path = save_path, filename = save_n,  dpi=600,
       width = 16,
       height = 17,
       units = c("cm"))
