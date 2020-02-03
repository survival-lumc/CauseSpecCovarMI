###############################
## Additional things to load ##
###############################


# Library
library(tidyverse)

# Set up scenario names
load("final_dat.RData")

# Throw away seed
final_dat <- final_dat %>% 
  select(-seed) %>%
  mutate_if(is.numeric, ~ round(., 3))

# Set up labels
scen_names <- unique(final_dat$scen_name) # extract names
betas <- as.numeric(str_extract(scen_names, "\\-*\\d+\\.*\\d*")) # get betas values
miss <- as.character(str_extract(scen_names, "high|low")) # high or low missing
mechs <- as.character(str_extract(scen_names, "MCAR|MAR_W|MAR_MRR")) # missingness mech

# Bind and paste
scen_labels <- apply(cbind(betas, miss, mechs), 1, function(x) {
  paste0("Beta = ", x[1], ", miss = ", x[2], ", mech = ", x[3])
})

# Assign names
scen_names <- setNames(scen_names, scen_labels)

rm(betas, mechs, miss)





# set theme
# Make general plotting theme
theme_MIsim_app <- function(base_family = "") {
  
  theme_light(base_family = base_family) %+replace%
    
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                    margin = margin(0, 0, 10, 0)),
          legend.position = "none",
          panel.spacing.x=unit(2, "lines"),
          panel.spacing.y=unit(1, "lines"), 
          strip.text = element_text(size = 16, colour = "black", face = "bold"),
          axis.title = element_text(size = 14),
          axis.text =  element_text(size = 12)) 
}
