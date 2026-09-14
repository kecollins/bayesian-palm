This repository contains the code used to produce the paper "Efficient Bayesian inference for spatial point patterns using the Palm likelihood". Please contact Kevin M. Collins at kmcolli9@ncsu.edu for any help or questions.

All computation was done using R version 4.4.1 running under macOS Sonoma 14.8.9. We include a `renv.lock` file to manage R package dependencies for reproducibility. Please install the package `renv` and run the command `renv::restore()` before any of the scripts below.

Note that all posterior sampling is done in Stan using the `rstan` package.

We recommend running the script `example.R` as the actual simulation studies below require days of computation time, even when run in parallel. The script `example.R` will generate a point pattern and fit a model for each of the processes discussed in the manuscript. It will output three `.pdf` files with traceplots from each model for all parameters. The actual code for model fitting can be found in `example_wrappers.R`. These examples require the package `doParallel`. We used 10 cores, but the number can be modified in `example.R` to suit the capabilities of your machine.

The folders `lgcp`, `thomas_process`, and `dpp` contain the files necessary to reproduce each of the simulation studies.

`lgcp`
- `lgcp_sim.R` runs the entire simulation study
- `lgcp_sim_results.R` will output the results found in tables in the manuscript and in the appendix

`thomas_process`
- `thomas_sim.R` runs the entire simulation study
- `thomas_sim_results.R` will output the results found in tables in the manuscript

`dpp`
- `dpp_sim.R` runs the entire simulation study
- `dpp_sim_results.R` will output the results found in tables in the manuscript

The folder `data-analysis` has the code necessary to reproduce the real data analysis. The dataset itself can be found in `spatstat::bei` and `spatstat::bei.extra`.

`data-analysis`
- `bei_analysis.R` will run the analysis and output the results found in the table in Section 5