##*************************************##
## Functions to make nested loop plots ##
##     and grouped lolly plots         ##
##          using ggplot2              ##
##*************************************##


# Read in estims file as test data
all_estims <- readRDS("analysis/simulations/all_estims.rds")

# Load nessary packages
devtools::load_all()
library(tidyverse)
library(RColorBrewer)
library(scales)

# For colors later
palo <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))


# NLP ---------------------------------------------------------------------

# Keep only relevant subset
dat_nlp_test <- all_estims %>% 
  filter(
    m %in% c(0, 50),
    var == "X.1",
    beta1 != "0",
    #miss_mech != "MCAR"
    miss_mech == "MAR",
    X_level == "continuous"
  )




# Parameters for a function
plot_nlp <- function(dat,
                     estim,
                     method_var,
                     true,
                     step_factors,
                     stretch = 1,
                     pointsize = 2,
                     point_dodge = 0.75,
                     text_size = 4,
                     step_labels = NULL) {
  
  # Make checks here
  
  # Treat step factors
  num_steps <- length(step_factors)
  step_list <- lapply(step_factors, function(step_var) dat[, step_var])
  
  # Create labels, if step_labels is NULL
  lab_steps <- sapply(step_factors, function(step_var) {
    steps <- levels(droplevels(dat[, step_var]))
    paste0(step_var, ": ", paste(steps, collapse = ", "))
  }, USE.NAMES = T)
  
  # Set up positions and bounds of steps
  range_data <- range(dat[, estim], na.rm = T)
  height_steps <- 0.1 * diff(range_data)
  height_betw_steps <- (2 / 3) * height_steps
  top_step <- range_data[1] - (1 / 4) * diff(range_data)
  bounds <- sapply(
    top_step - 0:(num_steps - 1) * (height_steps + height_betw_steps), 
    function(top) c(top - height_steps, top)
  )
  colnames(bounds) <- step_factors
  
  # Base data for the plot
  dat_nlp <- data.table::copy(dat) %>% 
    data.table::setDT() %>% 
    # Lex order true ensures steps biggest from large to small
    .[, scenario := as.numeric(
      interaction(step_list, lex.order = T, drop = T)
    ) * stretch]
  
  # Data-frame for steps and their labels
  steps_scaled <- sapply(step_factors, function(col) {
    scales::rescale(as.numeric(dat[, col]), to = bounds[, col])
  }) %>% 
    data.table::data.table() %>% 
    .[, scenario := dat_nlp$scenario] %>% 
    .[order(scenario)] 
  
  # Complete df, first expand last step
  dat_steps <- steps_scaled[.N, ] %>% 
    .[, scenario := scenario + stretch] %>% 
    rbind(., steps_scaled) %>% 
    
    # Put in long format
    data.table::melt.data.table(
      id.vars = "scenario", 
      measure.vars = step_factors,
      variable.name = "step_facts", 
      value.name = "scaled_vals" 
    ) %>% 
    
    # Remove duplicates
    .[, .SD[.N], by = .(scenario, step_facts)] %>% 
    
    # Put label in between steps
    .[, labpos := ifelse(
      scenario == min(scenario), 
      max(scaled_vals) + height_betw_steps / 2,
      NA
    ), by = step_facts] %>% 
    
    # Match the label
    .[scenario == min(scenario), ':=' (
      lab = lab_steps[match(step_facts, names(lab_steps))],
      lab_xpos = scenario + stretch / 5
    ), by = step_facts]
  
  # Prep plot breaks
  lowest_step_levels <- droplevels(dat[, step_factors[length(step_factors)]])
  len_seq <- length(levels(lowest_step_levels)) * stretch
  gridline_seq <- seq(len_seq, max(dat_nlp$scenario), 
                      by = len_seq)
  
  # Prep certain vars for plotting
  estim <- rlang::sym(estim)
  method_var <- rlang::sym(method_var)
  if (is.numeric(true)) {
    dat_nlp$true <- true
  } else true <- rlang::sym(true)

  
  p <- dat_nlp %>% 
    ggplot(aes(scenario + stretch / 2, !!estim)) +
    geom_text(
      data = dat_steps,
      aes(x = lab_xpos, y = labpos, label = lab),
      size = text_size,
      na.rm = T, hjust = 0
    ) +
    geom_point(
      aes(shape = !!method_var,  col = as.factor(scenario)),
      position = position_dodge(point_dodge * stretch), 
      size = pointsize, alpha = 0.8
    ) +
    geom_linerange(
      aes(ymin = !!true, ymax = !!estim, 
          xmin = scenario, xmax = scenario,
          col = as.factor(scenario),
          linetype = !!method_var, group = !!method_var),
      position = position_dodge(point_dodge * stretch)
    ) +
    geom_hline(aes(yintercept = !!true), linetype = "dashed") + 
    geom_step(
      data = dat_steps,
      aes(x = scenario, y = scaled_vals, group = step_facts)
    ) +
    scale_color_manual(values = palo(
      length(unique(dat_steps$scenario))
    )) +
    guides(colour = F) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", 
          axis.ticks.x = element_blank()) +
    coord_cartesian(expand = 0, 
                    ylim = c(min(bounds) - height_betw_steps, 
                             range_data[2] + height_betw_steps)) +
    scale_x_continuous(NULL, breaks = gridline_seq, labels = NULL) 
  
  return(p)
}

all_estims %>% 
  filter(
    m %in% c(0, 50),
    var == "X.1",
    beta1 != "0",
    #miss_mech != "MCAR"
    miss_mech == "MAR",
    X_level == "continuous"
  ) %>% 
  plot_nlp(
    dat = .,
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("prop_miss", "haz_shape", "eta1", "beta1"),
    stretch = 4,
    point_dodge = 0.8,
    text_size = 4, 
    pointsize = 2.5, 
    step_labels = NULL
  ) + theme_minimal()
  #xlab("Bias") +
  #facet_grid(miss_mech ~ X_level)


# To do:
# - break up into functions
# - add option of legend
# - add faceting
# - add warning if duplicated, co you add facets
# - Test, test
# - Make function for grouped lolly
# - Function for hybrid of lolly and candy with  bias


# (Grouped) lolly ---------------------------------------------------------


