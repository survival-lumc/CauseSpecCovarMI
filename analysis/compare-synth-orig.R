
library(data.table)
library(ggplot2)

dat_mds_merged <- rbind(dat_mds, dat_mds_synth, idcol = "df")
dat_mds_merged[, df := factor(df, levels = 1:2, labels = c("orig", "synth"))]

dat_mds_merged$df

dat_mds_merged %>% 
  ggplot(aes(x = cytog_threecat)) +
  geom_bar(aes(y = ..prop.., fill = df, group = df), position = "dodge")
# y=..prop..

dat_mds_merged %>% 
  ggplot(aes(x = age_allo1_decades)) +
  geom_histogram(aes(fill = df, group = df), position = "dodge", bins = 20)

# Better to cut into factor...
dat_mds_merged %>% 
  ggplot(aes(x = cut_interval(agedonor_allo1_decades, n = 15))) +
  geom_bar(aes(fill = df, group = df), position = "dodge")

# compare bivariate too?
plots_ls <- purrr::map(
  .x = names(dat_mds_merged)[-1],
  .f = ~ {
    is_fact <- is.factor(dat_mds_merged[[.x]])

    if (is_fact) {
      p <- dat_mds_merged %>% 
        ggplot(aes(x = .data[[.x]])) +
        geom_bar(aes(fill = df, group = df), position = "dodge")
    } else {
      p <- dat_mds_merged %>% 
        ggplot(aes(x = cut_interval(.data[[.x]], n = 15))) +
        geom_bar(aes(fill = df, group = df), position = "dodge") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2))
    }
    
    return(p)
  }
)

ggpubr::ggarrange(plotlist = plots_ls,
                  nrow = 7, common.legend = TRUE, ncol = 2)

# Make function for categorising - based on cut.. (with quantiles/labels, option
# for ordered etc.)
# or just ggplot2 cut_number?

