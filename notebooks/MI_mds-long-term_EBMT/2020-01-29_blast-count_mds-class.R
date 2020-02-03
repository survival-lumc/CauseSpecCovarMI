#********************************#
# Blast percentage and MDS class #
#********************************#


if (suppressWarnings(!require("pacman", character.only = T))) {
  install.packages("pacman", dep = T)
}; library(pacman)

pacman::p_load(foreign, tidyverse, ggpubr)


# Read-in EBMT data -------------------------------------------------------


# dat <- ...


# Check correspondence of MDS classes -------------------------------------


# Extract name of blast variables in data
blast_vars <- str_subset(names(dat), "bm_|pb_")

# Numerical summary
num_summary <- lapply(blast_vars, function(blast_var) {
  by(dat, dat$mdsclass, function(x) summary(x[, blast_var]))
})

names(num_summary) <- blast_vars

# Call for example: 
num_summary$bm_allo1


# Graphical summary - increase 'adjust' parameter for smoother density
blast_plots <- lapply(blast_vars, function(blast) {
  
  dat %>% 
    ggplot(aes(x = .[, blast], fill = mdsclass)) +
    geom_density(alpha= .5, adjust = 1, na.rm = T) +
    xlim(c(0, 75)) +
    ggtitle(blast) +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) 
})

# Make a grid
grid_blast_plots <- ggarrange(plotlist = blast_plots, 
                              nrow = 2, ncol = 2, 
                              common.legend = T, legend = "top")

annotate_figure(grid_blast_plots, 
                left = "Density", bottom = "Blast percentage")



# Is pmax(bm_allo1, bm_diag1) used to classify? ------------------------


# Take maximum and plot
dat %>% 
  mutate(bm_blast_max = pmax(bm_allo1, bm_diag1, na.rm = T)) %>% 
  ggplot(aes(x = bm_blast_max, fill = mdsclass)) +
  geom_density(alpha= .5, adjust = 1, na.rm = T) +
  xlim(c(0, 75)) +
  theme(legend.position = "top")


