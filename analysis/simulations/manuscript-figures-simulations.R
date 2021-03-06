##*********************************************##
## Make simulation result plots for manuscript ##
##*********************************************##

library(CauseSpecCovarMI)

# Manuscript plots --------------------------------------------------------


# Read-in summarised to save time
regr_results <- data.table::data.table(CauseSpecCovarMI::regr_results)
preds_results <- data.table::data.table(CauseSpecCovarMI::preds_results)

# Global theme
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12))

# Larger pointsize after review
ptsize <- 1.75


# Figure 1 ----------------------------------------------------------------


fig1 <- regr_results[
  m %in% c(0, 50) &
    n == 2000 & 
    beta1 != "0" & 
    miss_mech == "MAR" & 
    X_level == "continuous" & 
    var == "X.1" &
    !(analy %in% c("ch1", "ref"))
] %>% 
  ggplot_nlp(
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
    point_dodge = 0.7,
    text_size = 4, 
    pointsize = ptsize, 
    height_steps = 0.015,
    height_betw_steps = 0.03,
    top_step = -0.175,
    step_labels = c(
      "beta1" = "Beta1 = {0.5, 1}",
      "prop_miss" = "Missing % = {10, 50}", 
      "haz_shape" = "Hazard shapes = {similar, different}",
      "eta1" = "Mech. strength = {weak, strong}"
    )
  ) +
  ggplot2::scale_shape_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::scale_linetype_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend("Method"),
    linetype = ggplot2::guide_legend("Method")
  ) +
  ggplot2::scale_y_continuous(breaks = c(0, -0.05, -0.1, - 0.15)) +
  ggplot2::ylab("Bias") +
  ggplot2::theme(legend.position = "bottom")

fig1
ggplot2::ggsave(
  plot = fig1, 
  filename = "analysis/figures/MAR_B1_NLP.pdf", 
  units = "in",
  width = 6.24,
  height = 4.65
)


# Figure 2 ----------------------------------------------------------------


fig2 <- regr_results[
  m %in% c(0, 50) &
    n == 2000 & 
    beta1 != "0" & 
    miss_mech == "MAR_GEN" & 
    X_level == "continuous" & 
    var == "X.1" &
    !(analy %in% c("ch1", "ref"))
] %>% 
  ggplot_nlp(
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
    point_dodge = 0.7,
    text_size = 4, 
    pointsize = ptsize,
    top_step = -0.225,
    height_steps = 0.02,
    height_betw_steps = 0.04,
    step_labels = c(
      "beta1" = "Beta1 = {0.5, 1}",
      "prop_miss" = "Missing % = {10, 50}", 
      "haz_shape" = "Hazard shapes = {similar, different}",
      "eta1" = "Mech. strength = {weak, strong}"
    )
  ) +
  ggplot2::scale_shape_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::scale_linetype_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend("Method"),
    linetype = ggplot2::guide_legend("Method")
  ) +
  ggplot2::scale_y_continuous(breaks = c(0.1, 0, -0.1, -0.2)) +
  ggplot2::ylab("Bias") +
  ggplot2::theme(legend.position = "bottom")

fig2
ggplot2::ggsave(plot = fig2, filename = "analysis/figures/MAR-T_B1_NLP.pdf", 
                units = "in",
                width = 6.24,
                height = 4.65)


# Figure 3 ----------------------------------------------------------------


fig3 <- preds_results[
  m %in% c(0, 50) & 
    state != "EFS" &
    n == 2000 & 
    beta1 != "0" &
    miss_mech == "MAR" & 
    times == "5 years" &  
    `combo-X_Z` %in% c("-1SD_X-Z_-1SD") & 
    X_level == "continuous"& 
    prop_miss == "50%" &
    !(analy %in% c("ref", "ch1"))
] %>% 
  ggplot_nlp(
    estim = "rmse", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("state", "beta1", "haz_shape", "eta1"),
    text_size = 4, 
    pointsize = ptsize,
    point_dodge = 0.7,
    height_steps = 0.001,
    height_betw_steps = 0.003,
    step_labels = c(
      "beta1" = "Beta1 = {0.5, 1}", 
      "state" = "State = {REL, NRM}", 
      "haz_shape" = "Hazard shapes = {similar, different}", 
      "eta1" = "Mech. strength = {weak, strong}"
    ),
    top_step = -0.004
  ) +
  ggplot2::scale_shape_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::scale_linetype_discrete(
    breaks = c("CCA", "ch12", "ch12_int", "smcfcs"),
    labels = c("CCA", expression(CH[12]), expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend("Method"),
    linetype = ggplot2::guide_legend("Method")
  ) +
  ggplot2::ylab("Root mean square error (RMSE)") +
  ggplot2::scale_y_continuous(breaks = c(0, 0.01, 0.02)) +
  ggplot2::theme(legend.position = "bottom")


fig3
ggplot2::ggsave(plot = fig3, filename = "analysis/figures/MAR_pred_NLP.pdf", 
                units = "in",
                width = 6.24,
                height = 4.65)

