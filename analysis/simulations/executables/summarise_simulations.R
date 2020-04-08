##*************************************##
## Summarising simulation replications ##
##*************************************##


devtools::load_all()


# Estimates ---------------------------------------------------------------


# Read individual rds files and collapse into large DT
estims_allreps <- list.files(
  path = "./analysis/simulations/sim-reps_individual/",
  pattern = "estims*", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Extract scen_num, to split
  .[, scen_num := gsub(
    pattern = ".*(scen_num=)|(-rep).*$", 
    replacement = "", 
    x = scen_summary
  )] 


# Save each scenario all reps
estims_allreps %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>%  
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/simulations/sim-reps_all/estimates/", 
                    "estims_scen", .y, "_allreps.rds")
    ) 
  )
  

# Summarise estimates
estims_allreps %>% 
  
  #Compute bias, covarage, rejection rate + remove replicate/seed index
  .[, ':=' (
    bias = estimate - true,
    cover = `2.5 %` < true & true < `97.5 %`,
    rej = `p.value` < 0.05,
    scen_summary = gsub(
      pattern = "(-scen_num).*$", 
      replacement = "", 
      x = scen_summary
    )
  )] %>% 
  
  # Summarise quantities of interest 
  .[, .(
    n = .N,
    est = mean(estimate),
    se = mean(std.error),
    emp_se = sd(estimate),
    cover = mean(cover),
    bias = mean(bias),
    rej = mean(rej),
    warns = mean(warns)
  ), by = .(var, m, analy, true, scen_summary, scen_num)] %>% 
  
  # Compute monte carlo error for bias and coverage
  .[, ':=' (
    mcarlo_se_bias = emp_se / sqrt(n),
    mcarlo_se_cover = sqrt((cover * (1 - cover)) / n)
  )] %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>%  
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/simulations/sim-reps_summarised/estimates/", 
                    "estims_scen", .y, "_summarised.rds")
    ) 
  )
  
 
# Predictions -------------------------------------------------------------


# Read individual rds files and collapse into large DT
preds_allreps <- list.files(
  path = "./analysis/simulations/sim-reps_individual/",
  pattern = "preds*", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Extract scen_num, to split
  .[, scen_num := gsub(
    pattern = ".*(scen_num=)|(-rep).*$", 
    replacement = "", 
    x = scen_summary
  )] 
  

# Save allreps 
preds_allreps %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>% 
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/simulations/sim-reps_all/predictions/", 
                    "preds_scen", .y, "_allreps.rds")
    ) 
  )


# Saved summarised file
preds_allreps %>% 
  
  # Compute bias + remove redundant replicate/seed index
  .[, ':=' (
    bias = p_pool - true,
    scen_summary = gsub(
      pattern = "(-scen_num).*$", 
      replacement = "", 
      x = scen_summary
    )
  )] %>% 
  
  # Summarise quantities of interest
  .[, .(
    n = .N,
    prob = mean(p_pool),
    emp_se = sd(p_pool),
    bias = mean(bias), # maybe add relative bias?
    rmse = sqrt(mean(sq_err)),
    rmse_se = sd(sqrt(sq_err))
  ), by = .(
    analy, m,`combo-X_Z`, 
    times, state, 
    true, scen_summary, scen_num
  )] %>% 
  
  # Add monte carlo se of bias and rmse, keep scenario number
  .[, ':=' (
    mcarlo_se_bias = emp_se / sqrt(n),
    mcarlo_se_rmse = rmse_se / sqrt(n)
  )] %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>%  
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/simulations/sim-reps_summarised/predictions/", 
                    "preds_scen", .y, "_summarised.rds")
    ) 
  )

