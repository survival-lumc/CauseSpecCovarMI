devtools::load_all()
library(data.table)
library(ggplot2)
library(rsimsum)
library(tidyverse)
theme_set(theme_bw(base_size = 14))
source("analysis/supplement-simulations/scenarios-supplemental.R")
scenarios_raw[seq_len(2), ]

breslow_files <- list.files(
  path = "analysis/supplement-simulations/",
  pattern = "^suppl-sims_breslow*",
  full.names = TRUE
) 

breslow_scen_list <- lapply(breslow_files, function(file) {
  rbindlist(readRDS(file), idcol = "rep")
}) 

breslow_all <- rbindlist(breslow_scen_list, idcol = "scenario")
colnames(breslow_all)

breslow_summary <- breslow_all[, .(
  n = .N,
  est = mean(estimate),
  se = mean(std.error),
  se_mcse = sqrt(var(std.error^2) / (4 * .N * mean(std.error^2))),
  emp_se = sd(estimate),
  cover = mean(`2.5 %` < true & true < `97.5 %`),
  bias = mean(estimate - true),
  rmse = sqrt(mean((estimate - true)^2)),
  rmse_mcse = CauseSpecCovarMI::rmse_mcse(estimate, true, .N),
  warns = mean(warns)
), by = .(var, scenario, analy, true, beta1, X_level, miss_mech, eta1)] %>% # others later...
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n),
    cover_mcse = sqrt((cover * (1 - cover)) / n),
    beta1 = factor(beta1),
    analy = factor(
      analy, levels = rev(c("ref", "CCA", "imp_true", "imp_true_int",
                            "imp_nels", "imp_nels_int", "imp_breslow",
                            "imp_breslow_int", "smcfcs"))
    )
  )]


ggplot_lolly(
  dat = breslow_summary,
  estim = "bias",
  method_var = "analy",
  true = 0,
  mcarlo_se = "bias_mcse", 
  facet = "var ~ beta1"
) 


# Try multisimsum
simsum_way <- rsimsum::multisimsum(
  data = breslow_all,
  par = "var",
  estvarname = "estimate", 
  se = "std.error",
  by = "beta1",
  ref = "ref",
  methodvar = "analy",
  true = "true"
)


full_results <- map_dfr(
  .x = 5:160,
  .f = function(i) {
    s <- rsimsum::multisimsum(
      data = breslow_all[rep <= i],
      par = "var",
      estvarname = "estimate", 
      se = "std.error",
      by = "beta1",
      ref = "ref",
      methodvar = "analy",
      true = "true"
    )

    s <- summary(s)
    results <- tidy(s)
    results <- filter(results, stat == "bias" & var == "X.1") %>%
      mutate(i = i)
    return(results)
  }
)

library(hrbrthemes) # Provide some nicer themes and colour scales
ggplot(full_results, aes(x = i, y = est, color = analy)) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, fill = analy), alpha = 1 / 5) + 
  geom_line() +
  #scale_color_ipsum() +
  #theme_ipsum_rc(base_size = 12) +
  #theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
  labs(x = "Repetition #", y = "Estimated Bias", color = "Method",
       fill = "Method",
       title = "Bias over (cumulative) repetition number") +
  facet_wrap("beta1")


ggplot(full_results, aes(x = i, y = mcse, color = analy)) +
  geom_line() +
  geom_hline(yintercept = 0.01, color = "red", linetype = "dotted") +
  facet_wrap("beta1")


# Three comp --------------------------------------------------------------



threecomp_files <- list.files(
  path = "analysis/supplement-simulations/",
  pattern = "^suppl-sims_three*",
  full.names = TRUE
) 

threecomp_scen_list <- lapply(threecomp_files, function(file) {
  rbindlist(readRDS(file), idcol = "rep")
}) 

threecomp_all <- rbindlist(threecomp_scen_list, idcol = "scenario")
colnames(threecomp_all)

threecomp_summary <- threecomp_all[, .(
  n = .N,
  est = mean(estimate),
  se = mean(std.error),
  se_mcse = sqrt(var(std.error^2) / (4 * .N * mean(std.error^2))),
  emp_se = sd(estimate),
  cover = mean(`2.5 %` < true & true < `97.5 %`),
  bias = mean(estimate - true),
  rmse = sqrt(mean((estimate - true)^2)),
  rmse_mcse = CauseSpecCovarMI::rmse_mcse(estimate, true, .N),
  warns = mean(warns)
), by = .(var, scenario, analy, true, beta1, X_level, miss_mech, eta1)] %>% # others later...
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n),
    cover_mcse = sqrt((cover * (1 - cover)) / n),
    beta1 = factor(beta1),
    analy = factor(
      analy, levels = rev(c("ref", "CCA", "imp_ch123", "imp_ch123_int",
                            "smcfcs"))
    ),
    eta1 = factor(eta1)
  )]


# NLP plot
threecomp_summary |>#[var == "X.1"] |> 
  ggplot_nlp(
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("beta1", "X_level", "eta1"),
    point_dodge = 0.7,
    text_size = 4, 
    #pointsize = ptsize,
    #top_step = -0.225,
    #height_steps = 0.02,
    #height_betw_steps = 0.04,
    step_labels = c(
      "beta1" = "Beta1 = {0.5, 1}",
      "X_level" = "X type = {contin., binary}", 
      "eta1" = "Mech. strength = {weak, strong}"
    )
  ) +
  facet_grid(var ~ miss_mech)


threecomp_summary |>#[var == "X.1"] |> 
  ggplot_nlp(
    estim = "bias", 
    method_var = "analy", 
    true = 0, 
    step_factors = c("miss_mech", "beta1", "X_level", "eta1"),
    point_dodge = 0.7,
    text_size = 4, 
    #pointsize = ptsize,
    #top_step = -0.225,
    #height_steps = 0.02,
    #height_betw_steps = 0.04,
    step_labels = c(
      "miss_mech" = "Missing mech = {MAR, MAR-T}",
      "beta1" = "Beta1 = {0.5, 1}",
      "X_level" = "X type = {contin., binary}", 
      "eta1" = "Mech. strength = {weak, strong}"
    )
  ) +
  facet_grid(var ~ .)

threecomp_summary |> 
  ggplot(aes(bias_mcse)) +
  geom_histogram() +
  facet_wrap(X_level ~ var)


rsimsum::multisimsum(
  data = breslow_all,
  par = "var",
  estvarname = "estimate", 
  se = "std.error",
  by = "beta1",
  ref = "ref",
  methodvar = "analy",
  true = "true"
)




ggplot_lolly(
  dat = threecomp_summary,
  estim = "bias",
  method_var = "analy",
  true = 0,
  mcarlo_se = "bias_mcse", 
  facet = "var ~ beta1"
) 

ggplot_lolly(
  dat = threecomp_summary,
  estim = "rmse",
  method_var = "analy",
  true = 0,
  mcarlo_se = "rmse_mcse", 
  facet = "var ~ beta1"
) 

  
