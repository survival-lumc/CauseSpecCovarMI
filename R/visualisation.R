##********************************##
## Functions for general plotting ##
##********************************##


# Nested-loop plots -------------------------------------------------------


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
    step_labels <- vapply(
      X = step_factors, 
      FUN = function(step_var) {
        steps <- levels(droplevels(dat[[step_var]]))
        paste0(step_var, ": ", paste(steps, collapse = ", "))
      }, 
      FUN.VALUE = character(1), 
      USE.NAMES = TRUE
    )
  }
  
  # Set up height of steps and space between them, based on range
  range_data <- range(dat[[estim]], na.rm = TRUE)
  if (is.null(height_steps)) height_steps <- 0.1 * diff(range_data)
  if (is.null(height_betw_steps)) height_betw_steps <- (2 / 3) * height_steps
  if (is.null(top_step)) top_step <- range_data[1] - (1 / 4) * diff(range_data)
  
  # Compute bounds, 'top - heigh_steps' is lower bound of steps
  bounds <- vapply(
    X = top_step - 0:(num_steps - 1) * (height_steps + height_betw_steps), 
    FUN = function(top) c(top - height_steps, top),
    FUN.VALUE = numeric(2)
  )
  colnames(bounds) <- step_factors
  
  # Make list
  obj <- list(
    "bounds" = bounds,
    "betw_steps" = height_betw_steps,
    "labs" = step_labels,
    "num_steps" = num_steps,
    "range" = range_data
  )
  
  return(obj)
}


prep_nlp_data <- function(dat,
                          step_factors,
                          bounds_obj) {
  
  # For checks
  . <- scenario <- step_ID <- labpos <- scaled_vals <- NULL
  
  # Lex order true ensures steps biggest from large to small
  dat_nlp <- dat[, scenario := as.numeric(
    interaction(.SD, lex.order = TRUE, drop = TRUE)
  ), .SDcols = step_factors] 
  
  # Scale the steps based on bounds
  steps_scaled <- vapply(
    X = step_factors, 
    FUN = function(col) scales::rescale(as.numeric(dat[[col]]), to = bounds_obj$bounds[, col]),
    FUN.VALUE = numeric(nrow(dat))
  ) %>% 
    data.table::data.table()

  steps_scaled[, ':=' (
    scenario = dat_nlp$scenario,
    true = dat_nlp$true
  )]
  
  data.table::setorder(steps_scaled, scenario)

  # Complete df, first expand last step so steps don't end abruptly
  dat_steps <- steps_scaled[.N, ] 
  dat_steps[, scenario := scenario + 1]

  # Put in long format
  dat_steps_long <- data.table::melt.data.table(
    data = rbind(dat_steps, steps_scaled),
    id.vars = "scenario", 
    measure.vars = step_factors,
    variable.name = "step_ID", 
    value.name = "scaled_vals" 
  )  
  
  # Remove duplicates
  dat_steps_long <- dat_steps_long[, .SD[.N], by = .(scenario, step_ID)]
    
  # Put label in middle of first step
  dat_steps_long[, labpos := ifelse(
    scenario == min(scenario), 
    max(scaled_vals) + bounds_obj$betw_steps / 2,
    NA
  ), by = step_ID]  
    
  # Match label to its ID
  dat_steps_long[scenario == min(scenario), ':=' (
    lab = bounds_obj$labs[match(step_ID, names(bounds_obj$labs))],
    lab_xpos = scenario + (1 / 5)
  ), by = step_ID]
  
  # Make list
  obj <- list("df_main" = dat_nlp, "df_steps" = dat_steps_long)
  
  return(obj)
}



#' Custom ggplot-based nested loop plots
#'
#' @param dat Data frame with simulation results
#' @param estim Character for column containing estimates
#' @param method_var Character for column containing method
#' @param true Character for column containing true data-generating 
#' estimand value, or a scalar e.g. 0
#' @param step_factors Vector of step factors for bottom of the NL
#' @param pointsize Size of points
#' @param point_dodge Measure of horizontal dodge
#' @param text_size Size of the step labels
#' @param top_step Position of top step (in scale of data)
#' @param height_betw_steps Height between NLP step (on scale of data)
#' @param height_steps Height of the NLP steps themselves (on scale of data)
#' @param step_labels Optional vector of labels for the step functions.
#' Defaults to variable name and the factor labels
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' # Later here
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
  dat <- data.table::data.table(dat)
  
  # Read-in key parameters
  estim <- rlang::sym(estim)
  method_var <- rlang::sym(method_var)
  if (is.numeric(true)) dat$true <- true else true <- rlang::sym(true)
  
  # Create labels, if step_labels is NULL
  bounds_obj <- get_nlp_steps(
    dat = dat, 
    estim = estim, 
    step_factors, 
    step_labels = step_labels, 
    top_step = top_step,
    height_betw_steps = height_betw_steps,
    height_steps = height_steps
  )
  
  # Prepare data for plotting
  dat_obj <- prep_nlp_data(
    dat = dat,
    step_factors = step_factors,
    bounds_obj = bounds_obj
  )
  
  # Prepare gridline breaks based on steps and scenario
  levels_lowest_step <- droplevels(dat[[step_factors[bounds_obj$num_steps]]])
  gridline_by <- length(levels(levels_lowest_step))
  gridline_seq <- seq(
    from = gridline_by, 
    to = max(dat_obj$df_main$scenario), 
    by = gridline_by
  )
  
  # Begin plot
  p <- dat_obj$df_main %>% 
    ggplot2::ggplot(ggplot2::aes(x = .data$scenario + 0.5, y = !!estim)) +
    
    # Add step labels
    ggplot2::geom_text(
      data = dat_obj$df_steps,
      ggplot2::aes(x = .data$lab_xpos, y = .data$labpos, label = .data$lab),
      size = text_size,
      na.rm = TRUE, hjust = 0
    ) +
    
    # Add points
    ggplot2::geom_point(
      ggplot2::aes(shape = !!method_var,  col = as.factor(.data$scenario)),
      position = ggplot2::position_dodge(point_dodge), 
      size = pointsize#, alpha = 0.8
    ) +
    
    # Add dashed line for true
    ggplot2::geom_step(
      data = dat_obj$df_steps, 
      ggplot2::aes(x = .data$scenario, y = !!true), 
      linetype = "dashed"
    ) +
    
    # Add the vertical lines
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = !!true, 
        ymax = !!estim,
        xmin = .data$scenario, 
        xmax = .data$scenario,
        col = as.factor(.data$scenario),
        linetype = !!method_var, 
        group = !!method_var
      ),
      position = ggplot2::position_dodge(point_dodge),
      size = pointsize * .25
    ) +
    
    # There are the step functions
    ggplot2::geom_step(
      data = dat_obj$df_steps,
      ggplot2::aes(
        x = .data$scenario, 
        y = .data$scaled_vals, 
        group = .data$step_ID
      )
    ) +
    ggplot2::guides(colour = F) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.ticks.x = ggplot2::element_blank()
    ) +
    
    # Zoom in on plot
    ggplot2::coord_cartesian(
      expand = 0, 
      ylim = c(
        min(bounds_obj$bounds) - bounds_obj$betw_steps, 
        bounds_obj$range[2] + bounds_obj$betw_steps
      )
    ) +
    
    # Add proper gridlines
    ggplot2::scale_x_continuous(name = "Scenarios", breaks = gridline_seq, labels = NULL) 
  
  return(p)
}



# Forest plots ------------------------------------------------------------


#...
