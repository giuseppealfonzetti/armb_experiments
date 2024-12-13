
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Intro

This repository contains all the code needed to reproduce all the
simulations and real data experiments presented in the paper *Scalable
inference via averaged Robbins-Monro bootstrap*.

## How to use this repository

Clone the repository locally with

`git clone https://github.com/giuseppealfonzetti/armb_experiments`

Use [renv](https://rstudio.github.io/renv/articles/renv.html) to restore
the `R` environment:

``` r
renv::restore()
```

Run `make` to run the experiments pipeline from the console or use the
`Build All` button from Rstudio.

As can be seen in the `Makefile`, running the pipeline creates an
`output/` folder where all results will be stored, and a `data/` folder
where the data for the real data application will be downloaded.

## Description

The `R/` folder contains the scripts:

- `mle.R`;
- `ground_truth.R`;
- `bootstrap.R`;
- `blb.R`;
- `armb.R`;

Their names point to their aims. Note that the scripts are parameterised
in order to be used in all settings investigated in the paper. Data from
different settings is generated using the content of `R/sett/`.

The `R/` folder also contains the scripts:

- `sims_plots.R`, which reproduces the plots presented in the paper;
- `aps.R`, which reproduces the results of the real data application.

Finally, note that you can specify the number of cores to be used by
changing the object `num_cores` in `Makefile`.
