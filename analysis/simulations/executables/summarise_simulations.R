##*************************************##
## Summarising simulation replications ##
##*************************************##

devtools::load_all()

# Estimates ---------------------------------------------------------------


# Read individual rds files and collapse into large DT
list.files(
  path = "./analysis/results/estimates/",
  pattern = "*.rds", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist()  %>%
  
  # Compute bias, covarage, rejection rate 
  # + remove redundantreplicate/seed index
  .[, ':=' (
    bias = estimate - true,
    cover = `2.5 %` < true & true < `97.5 %`,
    rej = `p.value` < 0.05,
    scen_summary = gsub(
      pattern = "(-rep).*$", 
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
    rmse = rmse(estimate, true, .N),
    warns = mean(warns)
  ), by = .(var, m, analy, true, scen_summary)] %>% 
  
  # Compute monte carlo error for bias
  # Keep scenario number for file saving
  .[, ':=' (
    mcarlo_se_bias = emp_se / sqrt(n),
    scen_num = gsub(
      pattern = ".*(-scen_num=)", 
      replacement = "", 
      x = scen_summary
    )
  )] %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>%  # Replace $m by $scen_num
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/results/estimates_summarised/", 
                    "estims_scen", .y, "_summarised.rds")
    ) 
  )
  
 
# Predictions -------------------------------------------------------------


# Read individual rds files and collapse into large DT
list.files(
  path = "./analysis/results/predictions/",
  pattern = "*.rds", 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Compute bias + remove redundantreplicate/seed index
  .[, ':=' (
    bias = p_pool - true,
    scen_summary = gsub(
      pattern = "(-rep).*$", 
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
    rmse = mean(sqrt(sq_err))
  ), by = .(
    analy, m,`combo-X_Z`, 
    times, state, 
    true, scen_summary
  )] %>% 
  
  # Add monte carlo se of bias, keep scenario number
  .[, ':=' (
    mcarlo_se_bias = emp_se / sqrt(n),
    scen_num = gsub(
      pattern = ".*(-scen_num=)", 
      replacement = "", 
      x = scen_summary
    )
  )] %>% 
  
  # Split into list of dataframes based on scenario
  split(., .$scen_num) %>%  # Replace $m by $scen_num
  
  # Save each summarised scenario into a separate .rds file
  purrr::imap(
    ~ saveRDS(
      object = .x,
      file = paste0("./analysis/results/predictions_summarised/", 
                    "preds_scen", .y, "_summarised.rds")
    ) 
  )

