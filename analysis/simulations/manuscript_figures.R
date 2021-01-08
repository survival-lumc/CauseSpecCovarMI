devtools::load_all()

library(tidyverse)
library(scales)
library(RColorBrewer)

# Set ggplot theme
theme_set(
  theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
)

palo <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

all_estims <- readRDS("analysis/simulations//all_estims.rds")

all_preds <- readRDS("analysis/simulations/all_preds.rds") %>% 
  .[n != 500]


# Manuscript figures ------------------------------------------------------



dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    var == "X.1" &
    beta1 != "0" &
    miss_mech == "MAR" & 
    X_level == "continuous" &
    !(analy %in% c("ref", "ch1")) 
]

p <- ggplot_nlp(dat = dat_nlp,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
           point_dodge = 0.7,
           text_size = 4, 
           pointsize = 1.5, 
           height_steps = 0.01,
           height_betw_steps = 0.015,
           step_labels = c("beta1" = "Beta1 = {0.5, 1}", 
                           "prop_miss" = "Missing % = {10, 50}", 
                           "haz_shape" = "Hazard shapes = {similar, different}", 
                           "eta1" = "Mech. strength = {weak, strong}"),
           top_step = -0.15) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = c(0, -0.05, -0.1, - 0.15)) +
  theme(legend.position = "bottom")

p

ggsave(plot = p, filename = "MAR_B1_NLP.eps")

#library(latex2exp)
#test <- TeX("$\beta_1$")


dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    var == "X.1" &
    beta1 != "0" &
    miss_mech == "MAR_GEN" & 
    X_level == "continuous" &
    !(analy %in% c("ref", "ch1")) 
]

p <- ggplot_nlp(dat = dat_nlp,
                estim = "bias", 
                method_var = "analy", 
                true = 0, 
                step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
                point_dodge = 0.7,
                text_size = 4, 
                pointsize = 1.5, 
                height_steps = 0.02,
                height_betw_steps = 0.04,
                step_labels = c("beta1" = "Beta1 = {0.5, 1}", 
                                "prop_miss" = "Missing % = {10, 50}", 
                                "haz_shape" = "Hazard shapes = {similar, different}", 
                                "eta1" = "Mech. strength = {weak, strong}"),
                top_step = -0.225) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = c(0.1, 0, -0.1, - 0.2)) +
  theme(legend.position = "bottom")

p

ggsave(plot = p, filename = "MAR-T_B1_NLP.eps")


# Prediction
p <- all_preds[m %in% c(0, 50) & 
            state != "EFS" &
            beta1 != "0" &
            #eta1 %in% c("Weak") &
            miss_mech == "MAR" & 
            times == "5 years" &  
            `combo-X_Z` %in% c("-1SD_X-Z_-1SD") & 
            X_level == "continuous"& 
            prop_miss == "50%" &
            !(analy %in% c("ref", "ch1")) 
          
] %>% 
  #.[, bias := bias * 100] %>% 
  ggplot_nlp(dat = .,
             estim = "rmse", 
             method_var = "analy", 
             true = 0, 
             step_factors = c("state", "beta1", 
                              "haz_shape", "eta1"),# "eta1"), 
             text_size = 4, 
             pointsize = 1.5,
             point_dodge = 0.7,
             height_steps = 0.001,
             height_betw_steps = 0.003,
             step_labels = c("beta1" = "Beta1 = {0.5, 1}", 
                             "state" = "State = {REL, NRM}", 
                             "haz_shape" = "Hazard shapes = {similar, different}", 
                             "eta1" = "Mech. strength = {weak, strong}"),
             top_step = -0.004) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Root mean square error (RMSE)") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02)) +
  theme(legend.position = "bottom")


p

ggsave(plot = p, filename = "MAR_pred_NLP.eps")




# SE ----------------------------------------------------------------------

dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    #var == "Z.2" &
    beta1 != "0" &
    miss_mech == "MAR" & 
    X_level == "continuous" &
    !(analy %in% c("ref", "ch1")) 
]

p <- ggplot_nlp(dat = dat_nlp,
                estim = "cover", 
                method_var = "analy", 
                true = 0, 
                step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
                point_dodge = 0.7,
                text_size = 4, 
                pointsize = 1.5, 
                height_steps = 0.01,
                height_betw_steps = 0.015,
                step_labels = c("beta1" = "Beta1 = {0.5, 1}", 
                                "prop_miss" = "Missing % = {10, 50}", 
                                "haz_shape" = "Hazard shapes = {similar, different}", 
                                "eta1" = "Mech. strength = {weak, strong}"),
                top_step = -0.15) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(breaks = c(0, -0.05, -0.1, - 0.15)) +
  theme(legend.position = "bottom") +
  facet_wrap(~ var)

p
# Jeetje ------------------------------------------------------------------






dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    var == "Z.2" &
    beta1 != "0" &
    miss_mech != "MCAR" & 
    #X_level == "binary" &
    !(analy %in% c("ref", "ch1")) 
]

ggplot_nlp(dat = dat_nlp,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
           point_dodge = 0.75,
           text_size = 3, 
           pointsize = 2, 
           step_labels = NULL,
           top_step = NULL) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  facet_grid(miss_mech~ X_level)






dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    var == "Z" &
    beta1 != "0" &
    miss_mech == "MAR_GEN" & 
    X_level == "continuous" &
    !(analy %in% c("ref", "ch1")) 
]

ggplot_nlp(dat = dat_nlp,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
           point_dodge = 0.75,
           text_size = 3, 
           pointsize = 2, 
           step_labels = NULL,
           top_step = NULL) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias")  


# Try with facets here?
dat_nlp <- all_estims[
  m %in% c(0, 50) & 
    n == 2000 &
    #var == "X.1" &
    beta1 != "0" &
    miss_mech == "MAR" & 
    X_level == "continuous" &
    !(analy %in% c("ref", "ch1")) 
]

ggplot_nlp(dat = dat_nlp,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
           point_dodge = 0.75,
           text_size = 3, 
           pointsize = 2, 
           step_labels = NULL,
           top_step = NULL) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  facet_wrap(~ var)


#  Prediciton -------------------------------------------------------------






all_preds[m %in% c(0, 50) & 
            state != "EFS" &
            beta1 != "0" &
            #eta1 %in% c("Weak") &
            miss_mech == "MAR" & 
            times == "5 years" &  
            #`combo-X_Z` %in% c(#"1_X-Z_-1SD",
            #                   "+1SD_X-Z_+1SD") & 
            X_level == "continuous"& 
            prop_miss == "50%" &
            !(analy %in% c("ref", "ch1")) 
          
] %>% 
  ggplot_nlp(dat = .,
             estim = "rmse", 
             method_var = "analy", 
             true = 0, 
             step_factors = c("beta1", "state",
                              "haz_shape", "eta1"),# "eta1"), 
             text_size = 3, 
             pointsize = 2,
             top_step = -0.025,
             height_betw_steps = 0.025,
             height_steps = 0.015) +
  facet_wrap(~ `combo-X_Z`)






all_preds[m %in% c(0, 50) & 
            state != "EFS" &
            beta1 != "0" &
            #eta1 %in% c("Weak") &
            miss_mech == "MAR" & 
            times == "5 years" &  
            `combo-X_Z` %in% c("-1SD_X-Z_-1SD") & 
            X_level == "continuous"& 
            prop_miss == "50%" &
            !(analy %in% c("ref", "ch1")) 
          
] %>% 
  #.[, bias := bias * 100] %>% 
  ggplot_nlp(dat = .,
             estim = "rmse", 
             method_var = "analy", 
             true = 0, 
             step_factors = c("state", "beta1", 
                              "haz_shape", "eta1"),# "eta1"), 
             text_size = 3, 
             pointsize = 2) + 
             #top_step = -0.025 * 100,
             #height_betw_steps = 0.025 * 100,
             #height_steps = 0.015 * 100) +
  facet_wrap(~ `combo-X_Z`, nrow = 2)



