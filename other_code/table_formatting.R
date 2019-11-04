test <- final_dat %>% 
  
  # For one scenario first
  filter(scen_name == "0_high_MCAR") %>% 
  select(-seed, -scen_name) %>% 
  gather(key = measure, value, 4:11) %>% 
  unite(var_un, m, measure, sep = "_") %>% 
  spread(key = var_un, value) 

test
  
  
  