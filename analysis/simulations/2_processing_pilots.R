##************************************##
## Inital plotting of pilot scenarios ##
##************************************##


# Load compendium + tidyverse for plotting
devtools::load_all()

library(tidyverse)
library(RColorBrewer) 
theme_set(theme_bw(base_size = 14))


# Read-in summarised scenarios
pilots_summ <- list.files(
  path = "./analysis/sim_results/summarised_reps/estimates/",
  pattern = "*summarised.rds", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Adjust for eta minus label
  .[, scen_summary := gsub(
    "eta1=-", 
    "eta1=min", 
    scen_summary
  )] %>% 
  
  # Separate the scenario columns
  tidyr::separate(col = "scen_summary",
                  into = c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                           "rho", "eta1"),
                  sep = "-") %>% 
  mutate(analy = factor(
    analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
  ))


# Does SE / emp_SE decrease with more imputations? ------------------------


# Model-based SE
pilots_summ %>% 
  filter(m > 0,
         !(scen_num %in% c("28", "57")),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         beta1 == "beta1=1") %>% 
  ggplot(aes(m, se,
             col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")

# Empirical SEs
pilots_summ %>% 
  filter(m > 0,
         eta1 %in% c("eta1=min1", "eta1=NA"),
         beta1 == "beta1=1") %>% 
  ggplot(aes(m, emp_se,
             col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")




# Estimates ---------------------------------------------------------------


# With beta1 = 1
pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    !(scen_num %in% c("28", "57")),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
  ggplot(aes(analy, est, col = analy, shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = est - qnorm(0.975) * se,
                    ymax = est + qnorm(0.975) * se,
                    colour = analy), 
                width = 0.2,
                size = 1) +
  facet_grid(miss_mech ~ var) +
  guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")


# Using pointrage
pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    !(scen_num %in% c("28", "57")),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
  ggplot(aes(analy, est, col = analy)) +
             #shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_segment(aes(y = est - qnorm(0.975) * se,
                   yend = est + qnorm(0.975) * se,
                   x = analy,
                   xend = analy), 
               size = 2,
               alpha = .5) + 
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed")  +
  facet_grid(miss_mech ~ var) +
  #guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")

# With beta1 = 0
pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    beta1 == "beta1=0",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
  ggplot(aes(analy, est, col = analy, shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = est - qnorm(0.975) * se,
                    ymax = est + qnorm(0.975) * se,
                    colour = analy), 
                width = 0.2,
                size = 1) +
  facet_grid(miss_mech ~ var) +
  guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")


# Bias --------------------------------------------------------------------



# Actual bias
pilots_summ %>% 
  filter(m %in% c(0, 100),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         beta1 == "beta1=1") %>% 
  ggplot(aes(analy, bias, col = analy)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = bias, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  geom_point(aes(x = analy, y = bias + qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 41,
             size = 1.5) + 
  geom_point(aes(x = analy, y = bias - qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 40, 
             size = 1.5) + 
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")


# Absolute bias
pilots_summ %>% 
  filter(m %in% c(0, 100),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         beta1 == "beta1=1") %>%
  mutate(bias = abs(bias)) %>% 
  ggplot(aes(analy, bias, col = analy)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = bias, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  geom_point(aes(x = analy, y = bias + qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 41,
             size = 1.5) + 
  geom_point(aes(x = analy, y = bias - qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 40, 
             size = 1.5) + 
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")



# What about "strength" of the mechanism
pilots_summ %>% 
  filter(m %in% c(0, 100),
         miss_mech != "miss_mech=MCAR",
         beta1 == "beta1=1") %>%
  mutate(bias = abs(bias)) %>% 
  ggplot(aes(analy, bias, col = analy, 
             group = eta1, shape = eta1,
             linetype = eta1)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = bias, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  #geom_point(aes(x = analy, y = bias + qnorm(0.975) * mcarlo_se_bias),
  #           position = position_dodge(0.5), shape = 41,
  #           size = 1.5) + 
  #geom_point(aes(x = analy, y = bias - qnorm(0.975) * mcarlo_se_bias),
  #           position = position_dodge(0.5), shape = 40, 
  #           size = 1.5) + 
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")




# Coverage ----------------------------------------------------------------


# Plot coverage with monte carlo SE of coverage
pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    eta1 %in% c("eta1=min1", "eta1=NA"),
    beta1 == "beta1=1"
  ) %>% 
  mutate(mcarlo_se_cover = sqrt((cover * (1 - cover)) / 160)) %>% 
  ggplot(aes(y = cover, x = analy, col = analy)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_segment(aes(yend = cover + qnorm(0.975) * mcarlo_se_cover,
                      y = cover - qnorm(0.975) * mcarlo_se_cover,
                   xend = analy, x = analy),
                  size = 2, col = "grey", alpha = .75) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 2) + 
  theme(legend.position = 'none') +
  facet_grid(miss_mech ~ var) +
  ggtitle("Coverage")





# Predicitons -------------------------------------------------------------


# Read-in summarised scenarios
pilots_summ_preds <- list.files(
  path = "./analysis/sim_results/summarised_reps/predictions/",
  pattern = "*summarised.rds", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Adjust for eta minus label
  .[, scen_summary := gsub(
    "eta1=-", 
    "eta1=min", 
    scen_summary
  )] %>% 
  
  # Separate the scenario columns
  tidyr::separate(col = "scen_summary",
                  into = c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                           "rho", "eta1"),
                  sep = "-")



# Plot emp_se vs m
pilots_summ_preds %>% 
  filter(m > 0,
         eta1 %in% c("eta1=min1", "eta1=NA"),
         beta1 == "beta1=1") %>% 
  group_by(analy, m, miss_mech, times) %>% 
  summarise(emp_se = mean(emp_se)) %>% 
  ggplot(aes(m, emp_se,
             col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_grid(miss_mech ~ times) +
  scale_color_brewer(palette="Dark2")



# Try predictions
pilots_summ_preds %>% 
  filter(
    m %in% c(0, 100),
    eta1 %in% c("eta1=min1", "eta1=NA"),
    beta1 == "beta1=1",
    times == "5"
  ) %>% 
  ggplot(aes(x = analy, y = rmse,
             col = `combo-X_Z`, group = `combo-X_Z`)) +
  geom_point(size = 2) +
  geom_line(alpha = .5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_grid(miss_mech ~ state) +
  xlab("Analysis") +
  ylab("Root MSE State Probability")


dev.new()


# Grouped lolly, grouped by state, average patient
pilots_summ_preds %>% 
  filter(
    m %in% c(0, 100),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA"),
    `combo-X_Z` == "mean_X-Z_mean"
    #`combo-X_Z` == "-1SD_X-Z_-1SD"
  ) %>%
  ggplot(aes(analy, rmse, col = analy, 
             group = state, shape = state,
             linetype = state)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = rmse, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  facet_grid(miss_mech ~ times) +
  scale_color_brewer(palette="Dark2") +
  xlab("Analysis") +
  ylab("Root MSE State Probability") +
  ggtitle("For combo mean_X-Z_mean")

dev.new()

# Grouped lolly, grouped by state for another covar combo
pilots_summ_preds %>% 
  filter(
    m %in% c(0, 100),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA"),
    #`combo-X_Z` == "mean_X-Z_mean"
    `combo-X_Z` == "-1SD_X-Z_-1SD"
  ) %>%
  ggplot(aes(analy, rmse, col = analy, 
             group = state, shape = state,
             linetype = state)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = rmse, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  facet_grid(miss_mech ~ times) +
  scale_color_brewer(palette="Dark2") +
  xlab("Analysis") +
  ylab("Root MSE State Probability") +
  ggtitle("For combo -1SD_X-Z_-1SD")


# All reps ----------------------------------------------------------------




pilots_allreps <- list.files(
  path = "./analysis/sim_results/summarised_reps/estimates/",
  pattern = "*allreps.rds", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Implement this directly in the summarising file...
  
  # Remove redundant labelling
  .[, scen_summary := gsub(
    pattern = "(-scen).*$", 
    replacement = "", 
    x = scen_summary
  )] %>% 
  
  # Adjust for eta minus label
  .[, scen_summary := gsub(
    "eta1=-", 
    "eta1=min", 
    scen_summary
  )] %>% 
  .[, bias := estimate - true] %>% 
  
  # Separate the scenario columns
  tidyr::separate(col = "scen_summary",
                  into = c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                           "rho", "eta1"),
                  sep = "-")  
  



# Play with rsimsum -------------------------------------------------------



library(rsimsum)

simsum_scen2 <- simsum(
  data = pilots_allreps %>% 
    filter(m %in% c(0, 100),
           eta1 %in% c("eta1=min1", "eta1=NA"),
           beta1 == "beta1=1"),
  estvarname = "estimate", 
  se = "std.error",
  true = "true", 
  methodvar = "analy",
  ref = "ref",
  by = c("miss_mech", "var"),
  x = TRUE
)

autoplot(simsum_scen2, type = "nlp", stats = "bias")
autoplot(simsum_scen2, type = "est_ridge")
autoplot(summary(simsum_scen2), type = "lolly", stats = "bias")
autoplot(summary(simsum_scen2), type = "forest", stats = "cover")



# New plots ---------------------------------------------------------------

pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    !(scen_num %in% c("28", "57")),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
  mutate(analy = factor(
    analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
  )) %>% 
  ggplot(aes(analy, est, col = analy)) +
  #shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_segment(aes(y = est - qnorm(0.975) * se,
                   yend = est + qnorm(0.975) * se,
                   x = analy,
                   xend = analy), 
               size = 2,
               alpha = .5) + 
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed")  +
  facet_grid(miss_mech ~ var) +
  #guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")


# Stength of mechanism
pilots_summ %>% 
  filter(
    m %in% c(0, 100),
    !(scen_num %in% c("28", "57")),
    beta1 == "beta1=1",
    eta1 != "eta1=NA"
  ) %>% 
  ggplot(aes(analy, est, col = analy, group = eta1,
             shape = eta1)) +
  #shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_linerange(aes(ymin = est - qnorm(0.975) * se,
                   ymax= est + qnorm(0.975) * se,
                   xmin = analy,
                   xmax = analy), 
               size = 2,
               alpha = .5,  
               position = position_dodge(0.75)) + 
  geom_point(size = 2.5, position = position_dodge(0.75)) +
  geom_hline(aes(yintercept = true), linetype = "dashed")  +
  facet_grid(miss_mech ~ var) +
  scale_shape_manual(values = c(16, 15)) +
  #guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")


# Look at bias plot for small scenarios

# This is the old one
pilots_summ %>% 
  filter(m %in% c(0, 100),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         !(scen_num %in% c("28", "57")),
         beta1 == "beta1=1") %>%
  mutate(bias = abs(bias)) %>% 
  ggplot(aes(analy, bias, col = analy)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = bias, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  geom_point(aes(x = analy, y = bias + qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 41,
             size = 1.5) + 
  geom_point(aes(x = analy, y = bias - qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 40, 
             size = 1.5) + 
  facet_grid(miss_mech ~ var) +
  scale_color_brewer(palette="Dark2")


# New one

# For
pilots_summ %>% 
  filter(m %in% c(0, 100),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         (scen_num %in% c("28", "57")),
         beta1 == "beta1=1") %>%
  mutate(bias = bias) %>% 
  ggplot(aes(analy, bias, col = analy)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(position = position_dodge(0.75), size = 2) + 
  geom_linerange(aes(ymin = 0, ymax = bias, 
                     xmin = analy, xmax = analy),
                 position = position_dodge(0.75)) +
  coord_flip() + 
  geom_point(aes(x = analy, y = bias + qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 41,
             size = 1.5) + 
  geom_point(aes(x = analy, y = bias - qnorm(0.975) * mcarlo_se_bias),
             position = position_dodge(0.5), shape = 40, 
             size = 1.5) + 
  facet_grid(var ~ X_level) +
  scale_color_brewer(palette="Dark2") +
  ggtitle("For normal MAR")


pilots_summ %>% 
  filter(m %in% c(0, 100),
         eta1 %in% c("eta1=min1", "eta1=NA"),
         (scen_num %in% c("28", "57")),
         beta1 == "beta1=1") %>% 
  View()

