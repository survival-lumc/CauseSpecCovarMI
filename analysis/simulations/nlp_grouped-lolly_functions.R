##*************************************##
## Functions to make nested loop plots ##
##     and grouped lolly plots         ##
##          using ggplot2              ##
##*************************************##



# Read in estims file as test data
all_estims <- readRDS("analysis/simulations/all_estims.rds")

# Load nessary packages
devtools::load_all()
library(ggplot2)
library(RColorBrewer)
library(scales)

# For colors later
palo <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
theme_set(theme_bw(base_size = 14))


# NLP ---------------------------------------------------------------------

# Keep only relevant subset
dat_nlp_test <- all_estims[m %in% c(0, 50) & 
                           var == "X.1" &
                           beta1 != "0" &
                           miss_mech != "MCAR" #&
                           #X_level == "continuous"
                           ]


ggplot_nlp(dat = dat_nlp_test,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
           point_dodge = 1,
           text_size = 3, 
           pointsize = 1.5, 
           step_labels = NULL,
           top_step = NULL) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  facet_grid(X_level ~ miss_mech)
  #scale_y_continuous(breaks = seq(-0.15, 0.05, by = 0.05))
  
#facet_grid(X_level ~ miss_mech)




# To do:
# - break up into functions
# - add faceting
# - add warning if duplicated, co you add facets
# - Test, test



# (Grouped) lolly ---------------------------------------------------------



# Test grouped
all_estims[m %in% c(0, 50) & 
             haz_shape == "similar" &
             beta1 == "1" & 
             eta1 == "Strong"] %>% 

ggplot_lolly(dat = .,
             estim = "cove", 
             method_var = "analy", 
             true = 0.95, 
             group = "prop_miss", 
             facets = c("miss_mech * X_level", "var"),
             mcarlo_se = "mcarlo_se_cover",
             dodge = 1,
             scales = "free") +
  scale_color_brewer(palette="Dark2") +
  xlab("Methods") +
  ylab("Coverage (Monte-Carlo SE)") + 
  ggtitle("Coverage with haz_shape = similar, eta1 = strong") 



all_estims[m %in% c(0, 50) &
             haz_shape == "similar" &
             beta1 == "1" & 
             miss_mech == "MCAR" &
             X_level == "continuous" &
             var == "X.1" &
             prop_miss == "50%"] %>% 
  ggplot_lolly(dat = .,
               estim = "cover", 
               method_var = "analy", 
               true = 0.95, 
               mcarlo_se = "mcarlo_se_cover") +
  scale_color_brewer(palette="Dark2") +
  xlab("Methods") +
  ylab("Coverage (Monte-Carlo SE)")




# Candystick plot ---------------------------------------------------------






all_estims[m %in% c(0, 50) &
             haz_shape == "similar" &
             beta1 == "1" & 
             eta1 == "Strong" & 
             miss_mech != "MCAR" &
             prop_miss == "50%"] %>% 
  ggplot_estimates(dat = ., 
                   estim = "cover", 
                   method_var = "analy", 
                   se = "mcarlo_se_cover", 
                   true = 0.95, 
                   facets = c("var"), 
                   conf.int = 0.95,
                   lims_estim = c(0.75, 1)) +
  scale_color_brewer(palette="Dark2")




