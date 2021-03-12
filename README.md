# Multiple imputation for cause-specific Cox models: assessing methods for estimation and prediction

(To be submitted to Statistical Methods in Medical Research, SMMR)

- Abstract here with link to paper/pre-print

## Supplementary material

- [Generating simulated data](articles/data-generation.html)
- [Performance measures](articles/performance-measures.html)
- [Full results: regression coefficients](articles/regr-results-full.html)
- [Full results: predictions](articles/preds-results-full.html)

## Usage 

### General

This project is set-up as an [R-package compendium](https://github.com/ropensci/rrrpkg), with the following structure:

```
.
├── R   # User-made functions
├── analysis    # Scripts for illustrative analysis
│   ├── figures   # Manuscript figures
│   └── simulations   # Scripts to run simulation study
├── data    # Full and processed simulation data/synthetic MDS data
│   └── sim-reps_indiv    # Individual simulation replications
│       ├── preds
│       └── regr
├── data-raw    # Scripts for processing raw 
├── man   # Documentation user-made functions
└── vignettes   # Supplementary manuscript results
```

The `vignettes` directory contains the supplementary material from the manuscript, and is rendered here (add link) using [`pkgdown`](https://pkgdown.r-lib.org/).

To reproduce results from the manuscript, first install the package and its dependencies using

```R
remotes::install_github("survival-lumc/CauseSpecCovarMI", dependencies = TRUE)
```

You can then clone the directory using `git clone https://github.com/survival-lumc/CauseSpecCovarMI.git` , and set your RStudio working directory to the home of this repository (you can make use the included R project file).

### Reproducing the simulation study

The simulation study was run on a computer cluster using a [SLURM](https://slurm.schedmd.com/documentation.html) scheduler. If you have access to such a scheduler, you should install the package and clone the directory on the cluster first, and then run (from the home directory)

```bash
sbatch analysis/simulations/all-simulations.sh 
```

The above will run each replication of each scenario, and for each will save two `.rds` files (one for regression coefficients, another for the predictions) containing the results to the `data/sim-reps_indiv/`.  To summarise all of these, you can run

```bash
sbatch analysis/simulations/summarise-simulations.sh
```

This will produce four files, which represent the full and summarised results for both the regression coefficient results and those for the predictions:

- `data/sims_regr_full.fst` 
- `data/sims_regr_summary.fst` 
- `data/sims_preds_full.fst` 
- `data/sims_preds_summary.fst` 

These are included in the package/on Zenodo, such that one does not need to run the entire simulation study again (which would take on the order of days).

Should you want to replicate one replication of one scenario locally (for example the 10th replication of scenario 3), you can run

```bash
Rscript analysis/simulations/run-simulation.R 3 10
```

The scenarios (and their numbers/parameters) are contained in `data/scenarios.rds`.

### Reproducing manuscript figures (simulation study sections)

The figures can be reproduced locally by running

```sh
Rscript analysis/simulations/manuscript-figures-simulations.R
```

### Reproducing illustrative analysis

The original dataset on which the illustrative analysis was performed could not be shared, however a synthetic version is provided. This will yield similar results after the imputation procedure, but not comparable for complete-case analysis. On a cluster using SLURM, you can run

```bash
sbatch analysis/run-illustrative-analysis.sh
```

This will make use of 10 cores to run the 100 imputations (with 25 iterations) from both `mice()` and `smcfcs()` in parallel.

It can also be run locally by first changing the number of cores in `analysis/illustrative-analysis.R` , and then

```R
Rscript analysis/illustrative-analysis.R
```



## Authors

| Name                                                         | Affiliation                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Edouard F. Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807?setlanguage=English&setcountry=en) | Leiden University Medical Center (Leiden, NL)                |
| [Matthieu Resche-Rigon](https://www.researchgate.net/scientific-contributions/Matthieu-Resche-Rigon-56101026) | Paris Diderot University / Saint Louis Hospital (Paris, FR)  |
| [Johannes Schetelig](https://www.researchgate.net/scientific-contributions/Johannes-Schetelig-38769437) | Universitätsklinikum Dresden / DKMS Clinical Trials Unit (Dresden, DE) |
| [Hein Putter](https://www.lumc.nl/org/bds/medewerkers/hputter?setlanguage=English&setcountry=en) | Leiden University Medical Center (Leiden, NL)                |
| [Liesbeth C. de Wreede](https://www.lumc.nl/org/bds/medewerkers/lcdewreede?setlanguage=English&setcountry=en) | Leiden University Medical Center (Leiden, NL)                |

