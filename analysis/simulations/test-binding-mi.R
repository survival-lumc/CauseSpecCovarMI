library(data.table)

indiv_files <- list.files("all_estimates/estimates", full.names = T)

big_dat <- lapply(indiv_files, readRDS)

bigo <- rbindlist(big_dat)

library(pbapply)



# Preds -------------------------------------------------------------------


indiv_files <- list.files("all_preds/predictions/", full.names = T)

pboptions(type = "timer")
big_dat <- pblapply(indiv_files, readRDS)

bigo <- rbindlist(big_dat)

