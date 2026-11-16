cran_packages <- c("data.table", "grf", "ggplot2", "dplyr", 
                   "tidyr", "model4you", "remotes")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("causl", quietly = TRUE)) {
  remotes::install_github("rje42/causl")
}

if (!requireNamespace("rrcf", quietly = TRUE)) {
  remotes::install_github("vshirvaikar/rrcf", subdir = "r-package/rrcf")
}

all_packages <- c(cran_packages, "causl", "rrcf")

for (pkg in all_packages) {
  library(pkg, character.only = TRUE)
}
