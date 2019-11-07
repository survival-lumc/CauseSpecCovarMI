# Multiple imputation for competing risks

The is the GitLab directory for the 'Missing data for competing risks' project. For more information, consult the [wiki](https://git.lumc.nl/hp-ldw-eb/mi-competing-risks/wikis/home).

## Navigating the R directory

To begin with, clone (using SSH or HTTPS) the project in a directory of your choice. Thereafter, double-clicking on the `code_MI_comprisks.Rproj` file will open an R-studio session with the project directory set as current working directory. 

Open `simulations_main.R` , and load all the packages needed throughout the directory with the following code at the top of the file:

``` R
if (suppressWarnings(!require("pacman", character.only = T))) install.packages("pacman", dep = T); library(pacman)

pacman::p_load("MASS", "tidyverse", "mice", "smcfcs",
               "mitools", "survival", "foreign",
               "naniar", "VIM", "smcfcs", "cmprsk")
```

The following part of the R script allows you to navigate between sub-directories using relative paths. For example, let us assume you would like to run the shiny app allowing you to visualise the effect of different Weibull parameters on the simulated event times.

Step 1: Run `file.edit("shiny/shiny_weibull_params/parameter_cuminc.R")` . This will automatically open the main R file from the sub-directory.

Step 2: At the top of the opened file, run  `setwd("shiny/shiny_weibull_params/")` . This will set your working directory to that of the shiny app, and allow you to appropriately call any files sourced by the script.

Step 3: Run the shiny app, edit the file, commit any changes.

Step 4: Set working directory back to head directory, and switch back to main file using the following piece of code - 

```R
setwd("../../../code_MI_comprisks/")
file.edit("simulations_main.R")
```

Once back in the head directory, you can navigate again towards other sub-directories in a similar manner.