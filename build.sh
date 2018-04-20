#!/bin/bash

Rscript -e "Rcpp::compileAttributes(pkgdir='.');"
Rscript -e "roxygen2::roxygenise();"
Rscript -e "devtools::build_vignettes();"
ls | grep -P "^marathonR_[0-9\.]+tar\.gz" | xargs rm -f
R CMD build .

ls | grep -P "^marathonR_[0-9\.]+tar\.gz" | xargs R CMD INSTALL

R CMD BATCH devtest.R
cat devtest.Rout
rm -f devtest.Rout
rm -f .RData
