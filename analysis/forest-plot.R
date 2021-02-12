##***********************##
## Run analysis mds data ##
##***********************##


# Read-in
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% 
  data.table::setDT()

#

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 


# Prep formulas -----------------------------------------------------------


outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)]) 

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", rhs))
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", rhs))


# Load imputations --------------------------------------------------------


# Load in objects
imps_mice <- readRDS("data/imps-mds-mice.rds")
imps_smcfcs <- readRDS("data/imps-mds-smcfcs.rds")

# Prepare lists 
impdats_mice <- mice::complete(imps_mice, action = "all")
impdats_smcfcs <- imps_smcfcs$impDatasets

# Relapse models ----------------------------------------------------------

mice_rel <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)

smcfcs_rel <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)


# NRM models --------------------------------------------------------------


mice_nrm <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)

# smcfcs
smcfcs_nrm <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)


# CCA ---------------------------------------------------------------------


mod_rel <- survival::coxph(form_rel, data = dat_mds) %>% 
  summary(conf.int = 0.95) %$% 
  conf.int %>% 
  data.table::data.table(keep.rownames = "term") %>% 
  data.table::setnames(
    old = c("exp(coef)", "lower .95", "upper .95"), 
    new = c("estimate", "2.5 %", "97.5 %")
  )

mod_nrm <- survival::coxph(form_nrm, data = dat_mds) %>% 
  summary(conf.int = 0.95) %$% 
  conf.int %>% 
  data.table::data.table(keep.rownames = "term") %>% 
  data.table::setnames(
    old = c("exp(coef)", "lower .95", "upper .95"), 
    new = c("estimate", "2.5 %", "97.5 %")
  )






# Prepare forest plot data ------------------------------------------------


res_rel <- rbind(
  mod_rel %>% reflevels_add_summary(dat_mds, form_rel), 
  mice_rel %>% reflevels_add_summary(dat_mds, form_rel), 
  smcfcs_rel %>% reflevels_add_summary(dat_mds, form_rel), 
  fill = TRUE, 
  idcol = "method"
)

res_rel[, ':=' (
  method = factor(method, levels = 1:3, labels = c("CCA", "MICE", "SMC-FCS")),
  `exp(-coef)` = NULL,
  std.error = NULL,    
  statistic = NULL,       
  df = NULL,     
  p.value = NULL
)]

res_nrm <- rbind(
  mod_nrm %>% reflevels_add_summary(dat_mds, form_nrm), 
  mice_nrm %>% reflevels_add_summary(dat_mds, form_nrm), 
  smcfcs_nrm %>% reflevels_add_summary(dat_mds, form_nrm), 
  fill = TRUE, 
  idcol = "method"
)

res_nrm[, ':=' (
  method = factor(method, levels = 1:3, labels = c("CCA", "MICE", "SMC-FCS")),
  `exp(-coef)` = NULL,
  std.error = NULL,    
  statistic = NULL,       
  df = NULL,     
  p.value = NULL
)]

res_full <- rbind(res_rel, res_nrm, idcol = "state") 
res_full[, state := factor(state, levels = 1:2, labels = c("REL", "NRM"))]
res_full[, method := factor(method, levels = c("SMC-FCS", "MICE", "CCA"))]

# Add separate variables 
predictors_patt <- paste0("(", paste(predictors, collapse = "|"), ")")
res_full[, ':=' (
  var_name = regmatches(x = term, m = regexpr(pattern = predictors_patt, text = term)),
  levels = gsub(pattern = predictors_patt, replacement = "", x = term)
)]

#

data_dict <- readRDS("data/data-dictionary.rds")

# Add reference values
res_full[is.na(estimate), c("estimate", "2.5 %", "97.5 %") := NA]
res_full[levels == "", levels := NA_character_]


res_full <- merge(res_full, data_dict, by = c("var_name", "levels"))



# Forest ------------------------------------------------------------------

# Copy for testing
test_forest <- data.table::copy(res_full)

# Add rows for variable names of factors
new_rows <- test_forest[!is.na(levels_lab), .(var_name = unique(var_name)), 
                        by = .(var_label, method, state)]
test_forest <- rbind(test_forest, new_rows, fill = TRUE)
test_forest[is.na(level_num), level_num := 0]
test_forest[, levels_lab := ifelse(is.na(levels_lab), var_label, paste0("    ", levels_lab))]
data.table::setorder(test_forest, var_label, level_num, method, state)
test_forest[, levels_lab := factor(levels_lab, levels = unique(levels_lab))]

# test_forest the row colours 
test_forest[, graph_x := as.numeric(levels_lab)]
test_forest[order(levels_lab), colour_row := rep(
  c("white", "gray"), length.out = .N
), by = .(method, state)]

theme_set(theme_void(base_size = 14))

# Make the table for the left side
dat_table <- test_forest[, c("levels_lab", "count", "colour_row",
                             "graph_x")] %>% unique()

plot_table <- dat_table %>% 
  ggplot(aes(y = levels_lab)) + 
  geom_rect(
    aes(fill = colour_row), 
    ymin = dat_table$graph_x - 0.5, 
    ymax = dat_table$graph_x + 0.5,
    xmin = -Inf, 
    xmax = Inf
  ) +
  geom_text(aes(label = levels_lab), x= 0, hjust = 0) +
  geom_text(aes(label = count), x= 0.6, hjust = 0) +
  annotate("text", x = 0, y = 37, label = "Variable", hjust = 0, fontface = "bold") +
  annotate("text", x = 0.6, y = 37, label = "n", hjust = 0, fontface = "bold") +
  coord_cartesian(xlim = c(0, 0.75), ylim = c(0, 40)) +
  scale_fill_manual(values = c("gray90", "white"), guide = "none") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_y_discrete(limits = rev(levels(dat_table$levels_lab)))

# Plot for relapse
plot_dat_rel <- test_forest[state == "REL"]

plot_rel <- plot_dat_rel %>%
  ggplot(aes(x = levels_lab, y = estimate, group = method)) +
  scale_y_continuous("Hazard ratio (95% CI)",
                     trans = "log", breaks = c(0.5, 1, 1.5, 2, 3)) + 
  geom_rect(
    aes(fill = colour_row), 
    xmin = plot_dat_rel$graph_x - 0.5, 
    xmax = plot_dat_rel$graph_x + 0.5,
    ymin = -Inf, 
    ymax = +Inf
  ) +
  geom_linerange(
    aes(ymin = `2.5 %`, 
        ymax = `97.5 %`, 
        xmin = levels_lab, 
        xmax = levels_lab, 
        col = method),
    position = position_dodge(width = 0.75), 
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_point(
    aes(col = method, shape = method), 
    position = position_dodge(width = 0.75), 
    size = 1.25,
    na.rm = TRUE
  ) +
  coord_flip(ylim = c(0.65, 4), xlim = c(0, 40)) +
  geom_segment(x = 0, y = log(1), xend = 36, yend = log(1), linetype = "dashed") +
  geom_segment(x = 0, y = log(2), xend = 36, yend = log(2), linetype = "dotted") +
  annotate("text", x = 37, y = 1, label = "Relapse", hjust = 0.5, fontface = "bold") +
  theme(
    legend.position = "bottom",
    #axis.text.y = element_blank(), 
    #axis.ticks.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.ticks.x = element_line(colour = "grey10"),
    axis.ticks.length.x = unit(.15, "cm"),
    axis.line.x = element_line(colour = "grey10"),
    axis.text.x = element_text(colour = "black", size = 10),
  ) +
  scale_x_discrete(limits = rev(levels(test_forest$levels_lab))) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  scale_shape_discrete(na.translate = FALSE) +
  scale_fill_manual(values = c("gray90", "white"), guide = "none") +
  guides(
    colour=guide_legend("Method", byrow=TRUE, reverse = T), # nrow=2
    shape=guide_legend("Method",byrow=TRUE, reverse = T)
  )


plot_rel

# Plot for nrm
plot_dat_nrm <- test_forest[state == "NRM"]


plot_nrm <- plot_dat_nrm %>%
  ggplot(aes(x = levels_lab, y = estimate, group = method)) +
  #geom_segment(x = 0, y = 2, xend = 36, yend = 2) +
  scale_y_continuous(trans = "log", breaks = c(0.5, 1, 1.5, 2, 3)) +
  geom_rect(
    aes(fill = colour_row), 
    xmin = plot_dat_nrm$graph_x - 0.5, 
    xmax = plot_dat_nrm$graph_x + 0.5,
    ymin = -Inf, 
    ymax = +Inf
  ) +
  geom_linerange(
    aes(ymin = `2.5 %`, ymax = `97.5 %`, xmin = levels_lab, xmax = levels_lab, col = method),
    position = position_dodge(width = 0.75), 
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_point(
    aes(col = method, shape = method), 
    position = position_dodge(width = 0.75), 
    size = 1.25,
    na.rm = TRUE
  ) +
  coord_flip(ylim = c(0.5, 4), xlim = c(0, 40)) +
  geom_segment(x = 0, y = log(1), xend = 36, yend = log(1), linetype = "dashed") +
  geom_segment(x = 0, y = log(2), xend = 36, yend = log(2), linetype = "dotted") +
  annotate("text", x = 37, y = 1, label = "Non-relapse mortality", hjust = 0.5, fontface = "bold") +
  theme(
    legend.position = "bottom",
    #axis.text.y = element_blank(), 
    #axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = "grey10"),
    axis.ticks.length.x = unit(.15, "cm"),
    axis.line.x = element_line(colour = "grey10"),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.title = element_blank()
  ) +
  scale_x_discrete(limits = rev(levels(test_forest$levels_lab))) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  scale_shape_discrete(na.translate = FALSE) +
  scale_fill_manual(values = c("gray90", "white"), guide = "none") +
  guides(
    colour=guide_legend("Method", byrow=TRUE, reverse = T), # nrow=2
    shape=guide_legend("Method",byrow=TRUE, reverse = T)
  )

library(patchwork)
forests <- plot_table + 
  plot_rel + 
  plot_nrm + 
  plot_layout(guides = "collect", widths = c(1, 1.5, 1.5)) & 
  theme(legend.position = "top", plot.margin = margin(0, 0, 0, 2))#,
        #panel.grid = element_blank()) 

ggsave(filename = "test.png", plot = forests, width = 10, height = 11)




# New plot just REL -------------------------------------------------------


test_forest2 <- data.table::copy(test_forest)
# Make the table for the left side
dat_table <- test_forest2[, c("levels_lab", "count", "count_REL", "colour_row",
                             "graph_x")] %>% unique()

plot_table <- dat_table %>% 
  ggplot(aes(y = levels_lab)) + 
  geom_rect(
    aes(fill = colour_row), 
    ymin = dat_table$graph_x - 0.5, 
    ymax = dat_table$graph_x + 0.5,
    xmin = -Inf, 
    xmax = Inf
  ) +
  geom_text(aes(label = levels_lab), x= 0, hjust = 0) +
  geom_text(aes(label = count), x= 1, hjust = 1) +
  geom_text(aes(label = count_REL), x= 1.3, hjust = 1) +
  annotate("text", x = 0, y = 37, label = "Variable", hjust = 0, fontface = "bold") +
  annotate("text", x = 1, y = 37, label = "n", hjust = 1, fontface = "bold") +
  annotate("text", x = 1.3, y = 37, label = "# Events", hjust = 1, fontface = "bold") +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 40)) +
  scale_fill_manual(values = c("gray90", "white"), guide = "none") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_y_discrete(limits = rev(levels(dat_table$levels_lab)))

plot_table

# Plot for relapse
plot_dat_rel <- test_forest[state == "REL"]

plot_rel <- plot_dat_rel %>%
  ggplot(aes(x = levels_lab, y = estimate, group = method)) +
  scale_y_continuous("Hazard ratio (95% CI)",
                     trans = "log", breaks = c(0.5, 1, 1.5, 2, 3)) + 
  geom_rect(
    aes(fill = colour_row), 
    xmin = plot_dat_rel$graph_x - 0.5, 
    xmax = plot_dat_rel$graph_x + 0.5,
    ymin = -Inf, 
    ymax = +Inf
  ) +
  geom_linerange(
    aes(ymin = `2.5 %`, 
        ymax = `97.5 %`, 
        xmin = levels_lab, 
        xmax = levels_lab, 
        col = method),
    position = position_dodge(width = 0.75), 
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_point(
    aes(col = method, shape = method), 
    position = position_dodge(width = 0.75), 
    size = 1.25,
    na.rm = TRUE
  ) +
  coord_flip(ylim = c(0.65, 4.5), xlim = c(0, 40)) +
  geom_segment(x = 0, y = log(1), xend = 36, yend = log(1), linetype = "dashed") +
  geom_segment(x = 0, y = log(2), xend = 36, yend = log(2), linetype = "dotted") +
  annotate("text", x = 37, y = 1, label = "Relapse", hjust = 0.5, fontface = "bold") +
  theme(
    legend.position = "bottom",
    #axis.text.y = element_blank(), 
    #axis.ticks.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.ticks.x = element_line(colour = "grey10"),
    axis.ticks.length.x = unit(.15, "cm"),
    axis.line.x = element_line(colour = "grey10"),
    axis.text.x = element_text(colour = "black", size = 10),
  ) +
  scale_x_discrete(limits = rev(levels(test_forest$levels_lab))) +
  scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  scale_shape_discrete(na.translate = FALSE) +
  scale_fill_manual(values = c("gray90", "white"), guide = "none") +
  guides(
    colour=guide_legend("Method", byrow=TRUE, reverse = T), # nrow=2
    shape=guide_legend("Method",byrow=TRUE, reverse = T)
  )


plot_rel



library(patchwork)
forests <- plot_table + 
  plot_rel + 
  plot_layout(guides = "collect", widths = c(1, 1)) & 
  theme(legend.position = "top", plot.margin = margin(0, 0, 0, 2))#,
#panel.grid = element_blank()) 

forests

ggsave(filename = "forest_new.png", plot = forests, width = 10, height = 11)

