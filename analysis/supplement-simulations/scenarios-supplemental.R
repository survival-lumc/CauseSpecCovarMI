
scenarios_raw <- expand.grid(
  # Varied
  "beta1" = c(0.5, 1),
  "X_level" = c("continuous", "binary"),
  "miss_mech" = c("MAR", "MAR_GEN"),
  "eta1" = c(-1, -2),
  
  # Fixed
  "n" = 2000,
  "haz_shape" = "different",
  "prop_miss" = 0.5,
  "n_sim" = 160
)
