library(tidyverse)
theme_set(theme_bw(base_size = 18))



# Plotting scenario 1 (finished) ------------------------------------------


# Checks
estims_scen1 <- readRDS("analysis/results/estimates_summarised/estims_scen1_summarised.rds")
preds_scen1 <- readRDS("analysis/results/predictions_summarised/preds_scen1_summarised.rds")

# Warnings


# What do we do with a scen_summary, if we would like to group by other factors
# do this after summarising by group_by(scen_summary)
estims_scen1 %>% 
  tidyr::separate(col = "scen_summary",
                  into = c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                           "rho", "eta1", "scen_num_dupl"),
                  sep = "-") %>% 
  dplyr::select(-scen_num_dupl)



# Exploring estimates -----------------------------------------------------



# Estimates first
estims_scen1 %>%
  filter(m %in% c(0, 100)) %>% 
  ggplot(aes(analy, est)) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = est - 1.96 * se,
                    ymax = est + 1.96 * se,
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



# Does SE decrease with more imps? ----------------------------------------


# Plot empirical SE vs m
estims_scen1 %>% 
  filter(m > 0) %>% 
  pivot_longer(se:emp_se, 
               names_to = "se_type",
               values_to = "se") %>%
  group_by(m, se_type, var) %>% 
  summarise(se = mean(se)) %>% 
  ggplot(aes(m, se, col = se_type)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_wrap(. ~ var)


# Check out coverage
estims_scen1 %>% 
  filter(m %in% c(0, 100)) %>% 
  ggplot(aes(y = cover, x = analy, fill = analy)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_wrap(. ~ var) + 
  geom_bar(stat = "identity", alpha = .75) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  theme(legend.position = 'none') +
  ggtitle("Coverage")

# Plot rmse
estims_scen1 %>% 
  filter(m %in% c(0, 100)) %>% 
  ggplot(aes(x = analy, y = rmse, col = analy, group = 1)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_wrap(. ~ var)


# Trellis for predictions -------------------------------------------------



# Try predictions
preds_scen1 %>% 
  filter(m %in% c(0, 100)) %>% 
  ggplot(aes(x = analy, y = rmse,
             col = `combo-X_Z`, group = `combo-X_Z`)) +
  geom_point(size = 2) +
  geom_line(alpha = .5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_grid(state ~ times) +
  xlab("Analysis") +
  ylab("Root MSE State Probability")


# Use time as x axis
preds_scen1 %>% 
  filter(m %in% c(0, 100)) %>% 
  ggplot(aes(x = times, y = rmse, col = analy, group = analy)) +
  geom_point(size = 2) +
  geom_line(alpha = .5) +
  scale_x_continuous(guide = guide_axis(n.dodge = 2)) +
  facet_grid(`combo-X_Z` ~ state) 

