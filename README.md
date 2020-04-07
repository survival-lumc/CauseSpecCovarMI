# Missing data for cause-specific models

The is the GitLab directory for the 'Missing data for cause-specific models' project, in collaboration with EBMT staff and Paris Diderot University.

## Set-up

The code directory for this project is set-up as an R-package 'compendium'. Current stucture:

```
.
├── analysis
│   ├── data
│   │   ├── derived_data
│   │   └── raw_data
│   ├── notebooks
│   │   └── MI_mds-long-term_EBMT
│   ├── other
│   │   └── out_err_messages
│   ├── simulations
│   │   ├── executables
│   │   ├── sim-reps_all
│   │   │   ├── estimates
│   │   │   └── predictions
│   │   ├── sim-reps_individual
│   │   └── sim-reps_summarised
│   │       ├── estimates
│   │       └── predictions
│   └── visualisation
│       ├── shiny_sims_visuals
│       └── shiny_weibull_params
├── inst
│   └── testdata
├── man
├── R
└── vignettes
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
