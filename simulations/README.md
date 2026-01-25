# Simulations

This folder contains all code required to reproduce the simulation studies in  
**Shirvaikar, Storås, Lin and Holmes (2024)**.

The simulations are organized into three stages: data generation, simulation runs, and diagnostics.

### Data Generation

- **`CopulaRCT.R`**  
  Simulates randomized controlled trial (RCT) data using the *frugal parameterization*.

- **`CopulaObs.R`**  
  Simulates observational data using the frugal parameterization.

Both scripts rely on functionality from the <a href="https://github.com/rje42/causl">causl</a> package.

### Simulation Runs

- **`SimRCT.R`**  
  Runs the full RCT simulation study across 100 trials. The script logs summary statistics
  including ANOVA calibration p-values, variable importance measures, and related metrics.

- **`SimObs.R`**  
  Runs the corresponding observational data simulation study across 100 trials.

These scripts are the main entry points for reproducing the simulation results in the paper.

### Diagnostics and Visualization

- **`Diagnostics.R`**  
  Aggregates simulation outputs and generates diagnostic figures, including:
  - Empirical power and type I error rates  
  - Distributions of ANOVA calibration p-values  
  - Variable importance summaries  

All figures reported in the paper are produced using this script.

## Recommended Workflow

1. Run `CopulaRCT.R` and/or `CopulaObs.R` to verify and inspect data generation.
2. Execute `SimRCT.R` and/or `SimObs.R` to run the full simulation studies.
3. Run `Diagnostics.R` to aggregate results and produce summary plots and tables.
