# Simulations

This folder contains supplementary code for "Targeting relative risk heterogeneity with causal forests". The code is organized into five stages: setup, minimal examples (optional), data generation, simulation runs, and diagnostics.

### 1. Setup

Run this script once at the start of each session before running the scripts below.

- **`Setup.R`**: Installs missing package dependencies (CRAN and GitHub) and loads required libraries.

### 2. Minimal Examples (Optional)

- **`Examples.R`**: Demonstrates an end-to-end workflow: generate example data, fit the model, predict relative conditional treatment effects, and run an omnibus heterogeneity test. Contains implementations for all three methods discussed in the paper: relative risk causal forests (RRCF), generalized random forests (GRF), and model-based forests (MOB).

### 3. Data Generation

These scripts rely on functionality from the <a href="https://github.com/rje42/causl">causl</a> package.

- **`CopulaRCT.R`**: Simulates randomized controlled trial (RCT) data using the frugal parameterization.

- **`CopulaObs.R`**: Simulates observational data using the frugal parameterization.

### 4. Simulation Runs

These scripts are the main entry points for reproducing the simulation results in the paper.

- **`SimRCT.R`**: Runs the full RCT simulation study across 100 trials. The script logs summary statistics, including ANOVA calibration p-values, variable importance measures, and related metrics.

- **`SimObs.R`**: Runs the corresponding observational data simulation study across 100 trials.

### 5. Diagnostics

Figures reported in the paper are produced using this script.

- **`Diagnostics.R`**: Aggregates simulation outputs and generates diagnostic figures, including empirical power and type I error rates; distributions of ANOVA calibration p-values; and variable importance summaries.

