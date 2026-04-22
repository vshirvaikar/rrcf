# Simulations

This folder contains supplementary code for **Shirvaikar, Storås, Lin and Holmes (2024)**.

**`ExampleRRCF.R`** demonstrates how to generate example data, fit a relative risk causal forest, predict relative conditional treatment effects, and conduct an omnibus test for the overall detection of heterogeneity. To aid interested readers, **`ExampleMOB.R`** mirrors the same example data setup and workflow for model-based forests (MOB).

Code to reproduce the simulation study is organized into three stages: data generation, simulation runs, and diagnostics.

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

### Recommended Workflow

1. Run `ExampleRRCF.R` and/or `ExampleMOB.R` for a quick end-to-end example.
2. Run `CopulaRCT.R` and/or `CopulaObs.R` to verify and inspect data generation.
3. Execute `SimRCT.R` and/or `SimObs.R` to run the full simulation studies.
4. Run `Diagnostics.R` to aggregate results and produce summary plots and tables.
