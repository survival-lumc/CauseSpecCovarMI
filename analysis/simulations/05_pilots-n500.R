##************************************##
## Comparing pilots (n = 500 vs 2000) ##
##************************************##

devtools::load_all()
library(tidyverse)
library(scales)
library(RColorBrewer)

# Set ggplot theme
theme_set(
  theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
)

scenarios <- readRDS("inst/testdata/scenarios.rds")
palo <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))



# Read-in data ------------------------------------------------------------


all_estims <- readRDS("analysis/simulations/all_estims.rds")
all_preds <- readRDS("analysis/simulations/all_preds.rds")



# Variability of ests and preds -------------------------------------------

as.character(c(1:14, 120:133))

all_estims %>% 
  filter(m <= 100 & m != 0,
         scen_num %in% as.character(120:133)
) %>% 
  ggplot(aes(m, emp_se, 
             col = interaction(scen_num, analy)),
         group = interaction(scen_num, analy)) +
  geom_point(size = 2) +
  geom_line(size = 1, alpha = 0.5) +
  xlab("Number of imputed datasets (m)") +
  ylab("Mean model SE") +
  facet_grid(var ~ X_level, scales = "free") + #+
  theme(legend.position = "none")


all_preds[m <= 100 & 
            m != 0 &
            scen_num %in% as.character(120:131) &
            `combo-X_Z` %in% c("-1SD_X-Z_-1SD") &
            times == "5 years" &
            state == "NRM"] %>% 
  ggplot(aes(m, emp_se, 
             col = analy),
         group = analy) +
  geom_point(size = 2) +
  scale_x_continuous("Number of imputed datasets (m)", 
                     guide = guide_axis(n.dodge = 2)) +
  geom_line(size = 1, alpha = 0.5) +
  facet_wrap(. ~ scen_num, scales = "free") + 
  #xlab("Number of imputed datasets (m)") +
  ylab("Empirical SE predicted probabilities") #+
  #theme(legend.position = "none")


# The rest

all_estims %>% 
  .[m %in% c(0, 100) &
      scen_num %in% as.character(c(1:14, 120:133)) &
      miss_mech != "MCAR"
] %>% 
  .[, eta1 := droplevels(eta1)] %>% 
ggplot_nlp(dat = .,
           estim = "bias", 
           method_var = "analy", 
           true = 0, 
           step_factors = c("beta1", "n", "eta1"),
           point_dodge = 1,
           text_size = 3, 
           pointsize = 1.5, 
           step_labels = NULL,
           top_step = NULL) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  ylab("Bias") +
  facet_grid(miss_mech ~ var)
#scale_y_continuous(breaks = seq

