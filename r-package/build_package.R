# Builds and tests the rrcf package.
#
# To build the package for development:
#   `Rscript build_package.R`
#
# To prepare a CRAN build:
#   `Rscript build_package.R --as-cran`

args <- commandArgs(TRUE)
library(Rcpp)
library(devtools)
library(testthat)
library(roxygen2)

package.name <- "rrcf"

# If built for CRAN, exlude all test except ones with "cran" in the filename
# by adding the following regex to .Rbuildignore.
if (!is.na(args[1]) && args[1] == "--as-cran") {
  write_union("rrcf/.Rbuildignore", "^tests/testthat/test_((?!cran).).*")
  write_union("rrcf/.Rbuildignore", "^tests/testthat/data")
  write_union("rrcf/.Rbuildignore", "^tests/testthat/Rplots.pdf")
}

# Auto-generate documentation files
roxygen2::roxygenise(package.name)
## ISSUE: This does not add useDynLib(rrcf) to NAMESPACE - needs to be included

# Run Rcpp and build the package.
# Symlinks in `rrcf/src` point to the Rcpp bindings (`rrcf/bindings`) and core C++ (`core/src`).
# Note: we don't link in third_party/Eigen, because for the R package build we provide
# access to the library through RcppEigen.
compileAttributes(package.name)
clean_dll(package.name)
build(package.name)

# Test installation and run some smoke tests.
install(package.name)
library(package.name, character.only = TRUE)
# Treat warnings as errors.
options(warn = 2)
test_package(package.name)
