summary(dat)

plot(dat$t, dat$H1, type = "l")
lines(dat$t, dat$H2, col = "blue")

library(tidyverse)
theme_set(theme_bw(base_size = 18))

pooled_preds %>% 
  group_by(`combo-X_Z`, m, state, times, analy) %>% 
  summarise(mean_sq_err = sqrt(mean(sq_err))) %>% 
  filter(m == 4) %>% 
  filter(`combo-X_Z` == "-1SD_X-Z_mean") %>% 
  ggplot(aes(x = analy, y = mean_sq_err, col = analy, group = 1)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_grid(state ~ times)
