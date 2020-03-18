library(tidyverse)
theme_set(theme_bw(base_size = 18))

summary(dat)

plot(dat$t, dat$H1, type = "l")
lines(dat$t, dat$H2, col = "blue")

# Plot the parameter estimates
estimates %>%
  filter(m %in% c(10, 0)) %>% 
  ggplot(aes(analy, estimate)) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error,
                    colour = var), 
                width = 0.2,
                size = 1) +
  geom_point(size = 2) +
  coord_flip() +
  #ggtitle(scen_labels[which(scen_names == input$select_scen)]) +
  facet_wrap(. ~ var) +
  #scale_y_continuous(limits = c(-.75, 1.25), breaks = c(-.5, 0, .5, 1)) +
  xlab("Analysis") + ylab("Coefficient value") +
  theme(legend.position = "none")




pooled_preds %>% 
  group_by(`combo-X_Z`, m, state, times, analy) %>% 
  summarise(mean_sq_err = sqrt(mean(sq_err))) %>% 
  filter(m == 10) %>% 
  filter(`combo-X_Z` == "-1SD_X-Z_mean") %>% 
  ggplot(aes(x = analy, y = mean_sq_err, col = analy, group = 1)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_grid(state ~ times)



pooled_preds %>% 
  group_by(`combo-X_Z`, m, state, times, analy) %>% 
  summarise(mean_sq_err = sqrt(mean(sq_err))) %>% 
  filter(m == 10) %>% 
  #filter(`combo-X_Z` == "-1SD_X-Z_mean") %>% 
  filter(`combo-X_Z` == "mean_X-Z_mean" |
           `combo-X_Z` == "-1SD_X-Z_-1SD" |
           `combo-X_Z` == "-1SD_X-Z_+1SD" |
           `combo-X_Z` == "+1SD_X-Z_-1SD" |
           `combo-X_Z` == "+1SD_X-Z_+1SD" ) %>% 
  ggplot(aes(x = analy, y = mean_sq_err, col = `combo-X_Z`, group = `combo-X_Z`)) +
  geom_point(size = 2) +
  geom_line(alpha = .5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_grid(state ~ times) +
  xlab("Analysis") +
  ylab("Root MSE State Probability")


# What do we do with a scen_summary, if we would like to group by other factors
# do this after summarising by group_by(scen_summary)
lol <- estims_scen1_seed161 %>% 
  tidyr::separate(col = "scen_summary",
                  into = c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                           "rho", "eta1", "scen_num", "rep", "seed"),
                  sep = "-") %>% 
  mutate_if(is.character, ~ gsub(pattern = ".*=", # Remove label before =
                                 replacement = "", # Empty replace
                                 x = .))

