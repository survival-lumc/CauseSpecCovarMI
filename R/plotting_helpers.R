##*************************************##
##      Plotting helper functions:     ##
##   new NLP and grouped lolly plots + ##
##     other classic one               ##
##*************************************##


# To-do: add function documentation


# Simples/grouped lolly plots ---------------------------------------------


ggplot_lolly <- function(dat,
                         estim,
                         method_var,
                         true,
                         mcarlo_se = NULL, # This would be a MCSE
                         group = NULL,
                         facets = NULL,
                         conf.int = 0.95,
                         dodge = 0.75,
                         scales = "fixed",
                         lims_estim = NULL) {
  
  # Check if true is numeric or a variable
  if (is.numeric(true)) dat$true <- true
  else true <- rlang::sym(true)
  
  # Convert characters
  estim <- rlang::sym(estim)
  method_var <- rlang::sym(method_var)
  
  # Build basic lolly plot - base on groups
  if (!is.null(group)) {
    
    group <- rlang::sym(group)
    
    p <- dat %>% 
      ggplot2::ggplot(ggplot2::aes(x = !!method_var, 
                                   y = !!estim, 
                                   col = !!method_var,
                                   group = !!group, 
                                   shape = !!group, 
                                   linetype = !!group)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = !!true), 
                          linetype = "dashed") 
    
  } else {
    p <- dat %>%
      ggplot2::ggplot(ggplot2::aes(x = !!method_var, 
                                   y = !!estim, 
                                   col = !!method_var)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = !!true),
                          linetype = "dashed")
  }
  
  # Add points and line range
  p <- p +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(dodge)) + 
    ggplot2::geom_linerange(ggplot2::aes(ymin = !!true,
                                         ymax = !!estim, 
                                         xmin = !!method_var, 
                                         xmax = !!method_var),
                            position = ggplot2::position_dodge(dodge))
  
  
  # Add standard errors
  if (!is.null(mcarlo_se)) {
    
    mcarlo_se <- rlang::sym(mcarlo_se)
    crit <- qnorm((1 - conf.int) / 2, lower.tail = FALSE)
    
    p <- p + 
      ggplot2::geom_point(ggplot2::aes(x = !!method_var,
                                       y = !!estim + crit * !!mcarlo_se),
                          position = ggplot2::position_dodge(dodge), 
                          shape = 41,
                          size = 1.5) + 
      ggplot2::geom_point(ggplot2::aes(x = !!method_var, 
                                       y = !!estim - crit * !!mcarlo_se),
                          position = ggplot2::position_dodge(dodge), 
                          shape = 40,
                          size = 1.5)
  }
  
  # Check facets
  if (!is.null(facets)) {
    
    if (length(facets) > 2)
      stop("Facets should be of length 1 or 2. If you want to include another var,
           use an interaction as one of the facets e.g. X * Z")
    if (length(facets) == 2) {
      facets <- as.formula(paste(facets, collapse = " ~ ")) 
    } else facets <- as.formula(paste0(". ~ ", facets))
    
    p <- p + ggplot2::facet_grid(facets, scales = scales) 
  } 
  
  # Final touches
  p <- p + 
    ggplot2::coord_flip(ylim = lims_estim) +
    ggplot2::guides(colour = FALSE)
  
  return(p)
}


# New NLP plots -----------------------------------------------------------


get_nlp_steps <- function(dat,
                          estim,
                          step_factors,
                          step_labels,
                          top_step,
                          height_betw_steps,
                          height_steps) {
  
  # Extract number of step factors
  num_steps <- length(step_factors)
  
  # Create labels for the steps
  if (is.null(step_labels)) {
    step_labels <- sapply(step_factors, function(step_var) {
      
      steps <- levels(droplevels(dat[[step_var]]))
      label <- paste0(step_var, ": ", paste(steps, collapse = ", "))
      return(label)
    }, USE.NAMES = T)
  }
  
  # Set up height of steps and space between them, based on range
  range_data <- range(dat[[estim]], na.rm = T)
  
  if (is.null(height_steps)) height_steps <- 0.1 * diff(range_data)
  if (is.null(height_betw_steps)) height_betw_steps <- (2 / 3) * height_steps
  
  # Location of top of the first step
  if (is.null(top_step)) top_step <- range_data[1] - (1 / 4) * diff(range_data)
  
  # Compute bounds, 'top - heigh_steps' is lower bound of steps
  bounds <- sapply(
    top_step - 0:(num_steps - 1) * (height_steps + height_betw_steps), 
    function(top) c(top - height_steps, top)
  )
  colnames(bounds) <- step_factors
  
  # Make list
  obj <- list("bounds" = bounds,
              "betw_steps" = height_betw_steps,
              "labs" = step_labels,
              "num_steps" = num_steps,
              "range" = range_data)
  
  return(obj)
}


prep_nlp_data <- function(dat,
                          step_factors,
                          bounds_obj) {
  
  #' @import data.table
  
  # Lex order true ensures steps biggest from large to small
  dat_nlp <- dat[, scenario := as.numeric(
    interaction(.SD, lex.order = T, drop = T)
  ), .SDcols = step_factors] 
  
  # Scale the steps based on bounds
  steps_scaled <- sapply(step_factors, function(col) {
    scales::rescale(as.numeric(dat[[col]]), to = bounds_obj$bounds[, col])
  }) %>% 
    data.table::data.table() %>%
    
    # Add scenario and true, for plotting
    .[, ':=' (
      scenario = dat_nlp$scenario,
      true = dat_nlp$true
    )] %>% 
    .[order(scenario)] 
  
  # Complete df, first expand last step so steps don't end abruptly
  dat_steps <- steps_scaled[.N, ] %>% 
    .[, scenario := scenario + 1] %>% 
    rbind(., steps_scaled) %>% 
    
    # Put in long format
    data.table::melt.data.table(
      id.vars = "scenario", 
      measure.vars = step_factors,
      variable.name = "step_ID", 
      value.name = "scaled_vals" 
    ) %>% 
    
    # Remove duplicates
    .[, .SD[.N], by = .(scenario, step_ID)] %>% 
    
    # Put label in middle of first step
    .[, labpos := ifelse(
      scenario == min(scenario), 
      max(scaled_vals) + bounds_obj$betw_steps / 2,
      NA
    ), by = step_ID] %>% 
    
    # Match label to its ID
    .[scenario == min(scenario), ':=' (
      lab = bounds_obj$labs[match(step_ID, names(bounds_obj$labs))],
      lab_xpos = scenario + (1 / 5)
    ), by = step_ID]
  
  # Make list
  obj <- list("df_main" = dat_nlp,
              "df_steps" = dat_steps)
  
  return(obj)
}


# Parameters for a function
ggplot_nlp <- function(dat,
                       estim,
                       method_var,
                       true,
                       step_factors,
                       pointsize = 2,
                       point_dodge = 0.75,
                       text_size = 4,
                       top_step = NULL,
                       height_betw_steps = NULL,
                       height_steps = NULL,
                       step_labels = NULL) {
  
  
  # Coerce input to data.table (creates a copy)
  dat <- data.table::as.data.table(dat)
  
  # Read-in key parameters
  estim <- rlang::sym(estim)
  method_var <- rlang::sym(method_var)
  if (is.numeric(true)) {
    dat$true <- true
  } else true <- rlang::sym(true)
  
  # Create labels, if step_labels is NULL
  bounds_obj <- get_nlp_steps(dat = dat, 
                              estim = estim, 
                              step_factors, 
                              step_labels = step_labels, 
                              top_step = top_step,
                              height_betw_steps = height_betw_steps,
                              height_steps = height_steps)
  
  # Prepare data for plotting
  dat_obj <- prep_nlp_data(dat = dat,
                           step_factors = step_factors,
                           bounds_obj = bounds_obj)
  
  # Prepare gridline breaks based on steps and scenario
  levels_lowest_step <- droplevels(dat[[step_factors[bounds_obj$num_steps]]])
  gridline_by <- length(levels(levels_lowest_step))
  gridline_seq <- seq(from = gridline_by, 
                      to = max(dat_obj$df_main$scenario), 
                      by = gridline_by)
  
  
  # Begin plot
  p <- dat_obj$df_main %>% 
    ggplot2::ggplot(aes(x = scenario + 0.5, y = !!estim)) +
    
    # Add step labels
    ggplot2::geom_text(
      data = dat_obj$df_steps,
      ggplot2::aes(x = lab_xpos, y = labpos, label = lab),
      size = text_size,
      na.rm = T, hjust = 0
    ) +
    
    # Add points
    geom_point(
      ggplot2::aes(shape = !!method_var,  col = as.factor(scenario)),
      position = ggplot2::position_dodge(point_dodge), 
      size = pointsize#, alpha = 0.8
    ) +
    
    # Add dashed line for true
    ggplot2::geom_step(data = dat_obj$df_steps, 
                       ggplot2::aes(x = scenario, y = !!true), 
                       linetype = "dashed") +
    
    # Add the vertical lines
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = !!true, 
                   ymax = !!estim,
                   xmin = scenario, 
                   xmax = scenario,
                   col = as.factor(scenario),
                   linetype = !!method_var, 
                   group = !!method_var),
      position =  ggplot2::position_dodge(point_dodge),
      size = pointsize * .25
    ) +
    
    # There are the step functions
    ggplot2::geom_step(
      data = dat_obj$df_steps,
      ggplot2::aes(x = scenario, y = scaled_vals, group = step_ID)
    ) +
    ggplot2::guides(colour = F) +
    ggplot2::theme(legend.position = "bottom", axis.ticks.x = element_blank()) +
    
    # Zoom in on plot
    ggplot2::coord_cartesian(
      expand = 0, 
      ylim = c(min(bounds_obj$bounds) - bounds_obj$betw_steps, 
               bounds_obj$range[2] + bounds_obj$betw_steps)
    ) +
    
    # Add proper gridlines
    ggplot2::scale_x_continuous(name = "Scenarios", breaks = gridline_seq, labels = NULL) 
  
  return(p)
}


# 'Classic' estimates plots -----------------------------------------------


ggplot_estimates <- function(dat,
                             estim,
                             method_var,
                             se,
                             true,
                             facets = NULL,
                             scales = "fixed",
                             conf.int = 0.95,
                             lims_estim = NULL) {
  
  # Check if true is numeric or a variable
  if (is.numeric(true)) dat$true <- true
  else true <- rlang::sym(true)
  
  # Convert characters
  estim <- rlang::sym(estim)
  method_var <- rlang::sym(method_var)
  se <- rlang::sym(se)
  
  # Calculate critical val
  crit <- qnorm((1 - conf.int) / 2, lower.tail = FALSE)
  
  # Build basic plots
  p <- dat %>% 
    ggplot2::ggplot(ggplot2::aes(x = !!method_var,
                                 y = !!estim, 
                                 col = !!method_var)) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::geom_segment(aes(y = !!estim - crit * !!se,
                              yend = !!estim + crit * !!se,
                              x = !!method_var,
                              xend = !!method_var), 
                          size = 2,
                          alpha = .5)+ 
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(aes(yintercept = !!true), 
                        linetype = "dashed")  
  
  
  # Check facets
  if (!is.null(facets)) {
    
    if (length(facets) > 2)
      stop("Facets should be of length 1 or 2. If you want to include another var,
           use an interaction as one of the facets e.g. X * Z")
    if (length(facets) == 2) {
      facets <- as.formula(paste(facets, collapse = " ~ ")) 
    } else facets <- as.formula(paste0(". ~ ", facets))
    
    p <- p + ggplot2::facet_grid(facets, scales = scales) 
  } 
  
  # Final touches
  p <- p + 
    ggplot2::coord_cartesian(ylim = lims_estim) +
    ggplot2::guides(colour = FALSE)
  
  return(p)
}



# smcfcs: function to plot convergence ------------------------------------

ggplot_smcfcs_converg <- function(imps_obj,
                                  smformula,
                                  dat) {
  
  # Extract meta data
  m <- dim(imps_obj$smCoefIter)[1]
  iters <- dim(imps_obj$smCoefIter)[3]
  
  # Number of comp risk models
  K <- length(smformula)
  
  # Get column names for smcoefiter
  coef_names_list <- lapply(X = 1:K, FUN = function(k) {
    rhs <- gsub(x = smformula[k], pattern = ".*~", replacement = "")
    
    model_mat <- stats::model.matrix(
      object <- as.formula(paste0("~ 1 +", rhs)), 
      data <- dat
    )
    
    model_mat <- model_mat[, !(colnames(model_mat) %in% "(Intercept)")]
    
    coef_names_modk <- paste0(colnames(model_mat), ".", as.character(k))
  })
  
  # Unlist for names of smcoefiter
  coef_names <- unlist(coef_names_list)
  
  # Diagnostics of convergence\
  ests_list <- lapply(X = 1:m, function(i) {
    
    coef_dat <- as.data.frame(t(imps_obj$smCoefIter[i, ,]))
    coef_dat$iter <- 1:iters
    coef_dat$imp <- i 
    
    return(coef_dat)
  })
  
  ests_df <- data.table::setDT(do.call(rbind.data.frame, ests_list))
  
  # Set names
  data.table::setnames(ests_df, new = c(coef_names, "iter", "imp"))
  
  # Make plot
  p <- data.table::melt.data.table(
    data = ests_df, 
    variable.name = "covar",
    value.name = "value", 
    id.vars = c("imp","iter")
  ) %>% 
    ggplot2::ggplot(ggplot2::aes(iter, value, col = factor(imp))) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_wrap(~ covar) 
  
  return(p)
}

