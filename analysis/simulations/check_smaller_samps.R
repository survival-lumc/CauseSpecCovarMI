estims_scen28_seed4481 %>%
  filter(m > 0) %>% 
  ggplot(aes(m, std.error, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(. ~ var)

estims_scen57_seed9121 %>%
  filter(m > 0) %>% 
  ggplot(aes(m, std.error, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(. ~ var)


preds_scen28_seed4481 %>%
  filter(m > 0,
         `combo-X_Z` == "mean_X-Z_mean") %>% 
  group_by(state, times, analy, m) %>% 
  summarise(err = sqrt(mean(sq_err))) %>% 
  ggplot(aes(m, err, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(state ~ times)




# Read
estims_scen28_summarised %>% 
  filter(m > 0) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(. ~ var, scales = "free") +
  ylab("Pooled model SE")

estims_scen57_summarised %>% 
  filter(m > 0) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(. ~ var, scales = "free")


preds_scen57_summarised %>%
  filter(m > 0,
         `combo-X_Z` == "0_X-Z_mean",
         times == 10) %>% 
  group_by(state, times, analy, m) %>% 
  #summarise(err = sqrt(mean(sq_err))) %>% 
  ggplot(aes(m, emp_se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(state ~ ., scales = "free")
#facet_grid(state ~ times, scales = "free")


preds_scen28_summarised %>%
  filter(m > 0) %>% 
         #`combo-X_Z` == "mean_X-Z_mean",
         #times == 10) %>% 
  mutate(
    state = factor(state,
                   levels = c("1", "2", "3"),
                   labels = c("EFS", "REL", "NRM")),
    times = factor(times,
                   levels = c(0.5, 5, 10),
                   labels = c("6 mo.", "5y", "10y"))
  ) %>% 
  group_by(state, times, analy, m) %>% 
  summarise(se = mean(emp_se))  %>% 
  #summarise(err = sqrt(mean(sq_err))) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(state ~ times, scales = "free")


preds_scen28_summarised %>%
  filter(m > 0) %>% 
  #`combo-X_Z` == "mean_X-Z_mean",
  #times == 10) %>% 
  mutate(
    state = factor(state,
                   levels = c("1", "2", "3"),
                   labels = c("EFS", "REL", "NRM")),
    times = factor(times,
                   levels = c(0.5, 5, 10),
                   labels = c("6 mo.", "5y", "10y"))
  ) %>% 
  group_by(state, times, analy, m) %>% 
  summarise(se = mean(emp_se))  %>% 
  #summarise(err = sqrt(mean(sq_err))) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(state ~ times, scales = "free") +
  ylab("Empirical SE of pooled probabilities \n (all covariate combos)")

theme_set(theme_light(base_size = 14))


# Lets look at all reps
estims_scen28_allreps %>% 
  group_by(m, analy, var) %>% 
  summarise(n = n()) %>% 
  View()

estims_scen28_allreps %>% 
  filter(m > 0) %>% 
  #group_by(m, analy, var) %>% 
  ggplot(aes(std.error, fill = factor(m))) +
  geom_density(alpha = .5) +
  facet_grid(analy ~ var)

estims_scen57_allreps %>% 
  filter(m > 0) %>% 
  #group_by(m, analy, var) %>% 
  ggplot(aes(std.error, fill = factor(m))) +
  geom_density(alpha = .75) +
  facet_grid(analy ~ var)


preds_scen28_allreps %>% 
  filter(m > 0,
         times == 5,
         `combo-X_Z` == "mean_X-Z_mean") %>% 
  ggplot(aes(p_pool, fill = factor(m))) +
  geom_density(alpha = .75) +
  geom_vline(aes(xintercept = true), linetype = "dashed") +
  facet_grid(analy ~ state)

preds_scen28_allreps %>% 
  filter(m > 0) %>% 
  mutate(D = p_pool - true) %>% 
  group_by(state, m, analy, times) %>% 
  summarise(n = n(),
            Q10 = quantile(D, 0.1),
            Q90 = quantile(D, 0.9),
            diff = (Q90 - Q10)) %>% 
  ggplot(aes(m, diff, col = analy)) +
  geom_point() +
  geom_line() +
  facet_wrap(state ~ times, scales = "free")





# Plot H&L
estims_scen28_summarised %>% 
  filter(m > 0) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(. ~ var, scales = "free") +
  ylab("Pooled model SE")


estims_scen28_allreps %>% 
  filter(m > 0,
         var %in% c("X.1", "X.2")) %>% 
  #group_by(m, analy, var) %>% 
  ggplot(aes(std.error, fill = factor(m))) +
  geom_density(alpha = .5) +
  facet_grid(analy ~ var)


preds_scen28_summarised %>%
  filter(m > 0) %>% 
  #`combo-X_Z` == "mean_X-Z_mean",
  #times == 10) %>% 
  mutate(
    state = factor(state,
                   levels = c("1", "2", "3"),
                   labels = c("EFS", "REL", "NRM")),
    times = factor(times,
                   levels = c(0.5, 5, 10),
                   labels = c("6 mo.", "5y", "10y"))
  ) %>% 
  group_by(state, times, analy, m) %>% 
  summarise(se = mean(emp_se))  %>% 
  #summarise(err = sqrt(mean(sq_err))) %>% 
  ggplot(aes(m, se, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(state ~ times, scales = "free") +
  ylab("Empirical SE of pooled probabilities \n (all covariate combos)")


estims_scen57_summarised %>% 
  
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
  filter(
    m %in% c(0, 100),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
  ggplot(aes(analy, est, col = analy, shape = factor(m))) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = est - qnorm(0.975) * emp_se,
                    ymax = est + qnorm(0.975) * emp_se,
                    colour = analy), 
                width = 0.2,
                size = 1) +
  facet_grid(. ~ var) +
  guides(shape = FALSE) +
  scale_color_brewer(palette="Dark2")





estims_scen57_summarised %>% 
  
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
  filter(
    m %in% c(0, 100),
    beta1 == "beta1=1",
    eta1 %in% c("eta1=min1", "eta1=NA")
  ) %>% 
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
  facet_grid(. ~ var) +
  scale_color_brewer(palette="Dark2")
