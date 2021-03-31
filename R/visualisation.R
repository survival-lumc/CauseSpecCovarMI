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
#' @param step_factors Vector of variable names to use as step factors at 
#' bottom of the NLP.
#' @param pointsize Size of points
#' @param point_dodge Measure of horizontal dodge
#' @param text_size Size of the step labels
#' @param top_step Position of top step (on scale of data)
#' @param height_betw_steps Height between NLP step (on scale of data)
#' @param height_steps Height of the NLP steps themselves (on scale of data)
#' @param step_labels Optional named vector of labels for the step functions.
#' Defaults to variable name and the factor labels
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' 
#' sim_results <- CauseSpecCovarMI::regr_results
#' 
#' # Plot example code here...
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


# (Will be kept internal - otherwise needs function for data dictionary too)
# maybe later :)

#' Ggplot2-based grouped forest plot 
#' 
#' Used to visualise results from regression analyses. Function is not yet for 
#' general purpose use since it beforehand require a detailed 'data dictionary' data 
#' frame with variables labels/levels.
#' 
#' @param dat Dataframe on which the models were run
#' @param dictionary Data dictionary containing variable names/labels,
#' levels names/labels for categorical variables, general counts and event numbers
#' @param results Named list of data.frame each containing summarised regression results
#' 
#' @return Ggplot2 object (as returned by patchwork)
#' 
#' @noRd
ggplot_grouped_forest <- function(dat,
                                  dictionary,
                                  results,
                                  event = "Relapse", 
                                  form,
                                  breaks_x = c(0.5, 1, 1.5, 2, 3),
                                  lims_x = c(0.65, 4.5)) {
  
  # Prepare plot df
  df_plot <- prepare_forest_df(
    dat = dat,
    dictionary = dictionary,
    form = form,
    results = results
  )
  
  # Prepare labels and counts
  if (event == "Relapse") {
    event_counts <- "count_REL"
  } else event_counts <- "count_NRM"
  
  # Prepare df for table part of forest plot
  vars_table <- c("levels_lab", "count", event_counts, "graph_x", "colour_row")
  expr_table <- parse(text = paste0(".(", paste(vars_table, collapse = ","), ")"))
  table_df <- unique(df_plot[, eval(expr_table)]) 
  max_x <- max(table_df$graph_x) + 2
  
  # Make table
  table_plot <- table_df %>% 
    ggplot2::ggplot(ggplot2::aes(y = .data$levels_lab)) + 
    ggplot2::geom_rect(
      xmin = -Inf, 
      xmax = Inf,
      ggplot2::aes(
        fill = .data$colour_row,
        ymin = .data$graph_x - 0.5,
        ymax = .data$graph_x + 0.5
      )
    ) +
    ggplot2::geom_text(ggplot2::aes(label = .data$levels_lab), x= 0, hjust = 0, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$count), x= 1, hjust = 1, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = !!rlang::sym(event_counts)), x = 1.3, hjust = 1, na.rm = TRUE) +
    ggplot2::annotate("text", x = 0, y = max_x, label = "Variable", hjust = 0, fontface = "bold") +
    ggplot2::annotate("text", x = 1, y = max_x, label = "n", hjust = 1, fontface = "bold") +
    ggplot2::annotate("text", x = 1.3, y = max_x, label = "# Events", hjust = 1, fontface = "bold") +
    ggplot2::coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 40)) +
    ggplot2::scale_fill_manual(values = c("gray90", "white"), guide = "none") +
    ggplot2::scale_y_discrete(limits = rev(levels(table_df$levels_lab))) +
    ggplot2::theme_void(base_size = 14)
  
  # Make the right hand side
  estimates_plot <- df_plot %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$levels_lab, 
        y = .data$estimate, 
        group = .data$method
      )
    ) +
    ggplot2::scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      trans = "log", 
      breaks = breaks_x
    ) + 
    ggplot2::geom_rect(
      ymin = -Inf, 
      ymax = Inf,
      ggplot2::aes(
        fill = .data$colour_row,
        xmin = .data$graph_x - 0.5,
        xmax = .data$graph_x + 0.5
      )
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = `2.5 %`, 
        ymax = `97.5 %`, 
        xmin = .data$levels_lab, 
        xmax = .data$levels_lab, 
        col = .data$method
      ),
      position = position_dodge(width = 0.75), 
      size = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(col = .data$method, shape = .data$method), 
      position = position_dodge(width = 0.75), 
      size = 1.25,
      na.rm = TRUE
    ) +
    ggplot2::coord_flip(ylim = lims_x, xlim = c(0, 40)) +
    ggplot2::geom_segment(x = 0, y = log(1), xend = max_x - 1, yend = log(1), linetype = "dashed") +
    ggplot2::geom_segment(x = 0, y = log(2), xend = max_x - 1, yend = log(2), linetype = "dotted") +
    ggplot2::annotate("text", x = max_x, y = 1, label = event, hjust = 0.5, fontface = "bold") +
    ggplot2::scale_x_discrete(limits = rev(levels(table_df$levels_lab))) +
    ggplot2::scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
    ggplot2::scale_shape_discrete(na.translate = FALSE) +
    ggplot2::scale_fill_manual(values = c("gray90", "white"), guide = "none") +
    ggplot2::guides(
      colour = ggplot2::guide_legend("Method", byrow = TRUE, reverse = TRUE), 
      shape = ggplot2::guide_legend("Method", byrow = TRUE, reverse = TRUE)
    ) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(colour = "black", size = 12),
      axis.ticks.x = ggplot2::element_line(colour = "grey10"),
      axis.ticks.length.x = ggplot2::unit(.15, "cm"),
      axis.line.x = ggplot2::element_line(colour = "grey10"),
      axis.text.x = ggplot2::element_text(colour = "black", size = 10)
    )
  
  # Bring together with patchwork
  forest_plot <- table_plot + 
    estimates_plot + 
    patchwork::plot_layout(guides = "collect", widths = c(1, 1)) & 
    ggplot2::theme(legend.position = "top", plot.margin = ggplot2::margin(0, 0, 0, 2))
  
  return(forest_plot)
}

prepare_forest_df <- function(dat,
                              dictionary, 
                              form,
                              results) {
  
  # Combine list of regression results - add also NA row with reference categories
  res <- data.table::rbindlist(
    lapply(X = results, FUN = function(mod_summary) {
      reflevels_add_summary(summ = mod_summary, dat = dat, form = form)
    }), 
    fill = TRUE,
    idcol = "method"
  )
  
  res[, method := factor(method, levels = names(results))]
  
  # Get predictors from RHS
  patt <- paste0("(", paste(unique(dictionary$var_name), collapse = "|"), ")")
  
  # Set var_name and levels according to data dict?
  # Possibly along with eventual function create_data_dictionary()
  res[, ':=' (
    var_name = regmatches(x = term, m = regexpr(pattern = patt, text = term)),
    levels = gsub(pattern = patt, replacement = "", x = term)
  )]
  
  res[levels == "", levels := NA_character_]
  df_merged <- merge(res, dictionary, by = c("var_name", "levels"))
  
  # Add rows for variable names of factors (once for each method)
  new_rows <- df_merged[!is.na(levels_lab), list(
    var_name = unique(var_name)
  ), by = list(var_label, method)]
  
  df_merged_expanded <- rbind(df_merged, new_rows, fill = TRUE)
  df_merged_expanded[is.na(level_num), level_num := 0]
  
  # Prepare the left most column of forest table, use tabs for levels
  df_merged_expanded[, levels_lab := ifelse(
    is.na(levels_lab), 
    var_label, 
    paste0("    ", levels_lab)
  )]
  
  # Sort a) alphabetic var labels, then by levels order, then method
  data.table::setorder(df_merged_expanded, var_label, level_num, method)
  
  # Unique() will take order they were sorted in 
  df_merged_expanded[, levels_lab := factor(levels_lab, levels = unique(levels_lab))]
  
  # Prepare the alternating grey/white rows
  df_merged_expanded[, graph_x := as.numeric(levels_lab)]
  df_merged_expanded[order(levels_lab), colour_row := rep(
    x = c("white", "gray"), length.out = .N
  ), by = .(method)]
  
  return(df_merged_expanded)
}
