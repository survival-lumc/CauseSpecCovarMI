# Missing data for cause-specific models

The is the GitLab directory for the 'Missing data for cause-specific models' project, in collaboration with EBMT staff and Paris Diderot University.

## Set-up

The code directory for this project is set-up as an R-package 'compendium'. Current stucture:

```powershel
MI-cause-specific_R
    ├───analysis
    │   ├───data
    │   │   ├───derived_data
    │   │   └───raw_data
    │   ├───figures
    │   ├───paper
    │   └───visualisation
    │       ├───shiny_sims_visuals
    │       └───shiny_weibull_params
    ├───doc
    ├───inst
    │   └───testdata
    ├───man
    ├───Meta
    ├───R
    ├───temp
    │   ├───notebooks
    │   │   ├───MI-cytoscoring_EBMT
    │   │   └───MI_mds-long-term_EBMT
    │   │       └───figures
    │   ├───old_helpers
    │   └───shark
    │       ├───data_generation
    │       └───support_functions
    └───vignettes
```



You need to install the `devtools` upon cloning the directory to run all code. To do so, call `devtools::load_all()`.

## Future

The current repository will eventually be available of github, ideally also with a docker container.

## Contributors

| Name                  | Institute     | Role         |
| --------------------- | ------------- | ------------ |
| Edouard Bonneville    | LUMC          | Maintainer   |
| Hein Putter           | LUMC          | Supervisor   |
| Liesbeth de Wreede    | LUMC          | Supervisor   |
| Matthieu Resche-Rigon | Paris Diderot | Collaborator |
