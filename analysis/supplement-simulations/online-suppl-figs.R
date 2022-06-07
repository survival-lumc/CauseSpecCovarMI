# 2.1.2
p1 <- ggplot_nlp(
  dat = regr_results[m %in% c(0, 50) & n == 2000 & beta1 != "0" & 
                       miss_mech == "MCAR" & X_level == "continuous"],
  estim = "rmse", 
  method_var = "analy", 
  true = 0, 
  step_factors = c("beta1", "prop_miss", "haz_shape", "eta1"),
  point_dodge = 0.7,
  text_size = 4, 
  pointsize = 1.5, 
  height_steps = 0.03,
  height_betw_steps = 0.075,
  top_step = -0.1,
  step_labels = c(
    "beta1" = "Beta1 = {0.5, 1}",
    "prop_miss" = "Missing % = {10, 50}", 
    "haz_shape" = "Hazard shapes = {similar, different}",
    "eta1" = "Mech. strength = {weak, strong}"
  )
) +
  ggplot2::guides(shape = guide_legend("Method"),
                  linetype = guide_legend("Method")) +
  theme_bw(base_size = 12) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  facet_wrap(. ~ var_label, nrow = 4, ncol = 1, labeller = label_parsed) +
  ylab("RMSE") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)#,
    #strip.text = element_text(size=15)
  ) +
  ggplot2::scale_shape_discrete(
    breaks = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs"),
    labels = c("Ref", "CCA", expression(CH[1]), expression(CH[12]), 
               expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::scale_linetype_discrete(
    breaks = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs"),
    labels = c("Ref", "CCA", expression(CH[1]), expression(CH[12]), 
               expression(CH["12,int"]), "smcfcs")
  ) 

ggplot2::ggsave(
  plot = p1, 
  filename = "analysis/figures/github_212_regr.pdf", 
  units = "in",
  width = 7, 
  height = 10, 
  dpi = 300
) 

# 1.2.1
p2 <- preds_results[
  beta1 != "0" &
    miss_mech == "MAR" & 
    times == "5 years" &  
    `combo-X_Z` %in% c("mean_X-Z_mean", "-1SD_X-Z_-1SD", "+1SD_X-Z_+1SD") & 
    X_level == "continuous" & 
    prop_miss == "50%"
] %>% 
  ggplot_nlp(
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("state", "beta1", "haz_shape", "eta1"),
    text_size = 4, 
    pointsize = 1.5,
    point_dodge = 0.7,
    height_steps = 0.005,
    height_betw_steps = 0.01,
    step_labels = c(
      "beta1" = "Beta1 = {0.5, 1}", 
      "state" = "State = {REL, NRM}", 
      "haz_shape" = "Hazard shapes = {similar, different}", 
      "eta1" = "Mech. strength = {weak, strong}"
    ),
    top_step = -0.025
  ) +
  facet_wrap(
    . ~ `combo-X_Z`, nrow = 3, ncol = 1,
    labeller = labeller(`combo-X_Z` = patient_labs)
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend("Method"),
    linetype = ggplot2::guide_legend("Method")
  ) +
  ggplot2::ylab("Bias") +
  ggplot2::scale_y_continuous(breaks = c(-0.025, 0, 0.025)) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_shape_discrete(
    breaks = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs"),
    labels = c("Ref", "CCA", expression(CH[1]), expression(CH[12]), 
               expression(CH["12,int"]), "smcfcs")
  ) +
  ggplot2::scale_linetype_discrete(
    breaks = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs"),
    labels = c("Ref", "CCA", expression(CH[1]), expression(CH[12]), 
               expression(CH["12,int"]), "smcfcs")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


ggplot2::ggsave(
  plot = p2, 
  filename = "analysis/figures/github_121_pred.pdf", 
  units = "in",
  width = 7, 
  height = 10, 
  dpi = 300
) 
