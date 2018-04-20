## marathonR

This R package brings the highly-optimized sampling routines from the
[marathon](https://github.com/srechner/marathon) library to R.

It provides methods for the uniform sampling of binary matrices with restricted row and column sums.
In particular, it contains functions for the uniform sampling of 
  a) binary matrices with fixed row and column sums.
  b) binary matrices whose row and column sums are restricted by lower and upper bounds.

When using these methods in scientific publications, please cite

Steffen Rechner and Annabell Berger: marathon: An open source software library
for the analysis of Markov-Chain Monte Carlo algorithms , PLoS ONE 2016.

### Installation

#### ubuntu

Tested on: Ubuntu 16.04, Xubuntu 17.04

The package requires the following dependencies:

-   Ubuntu
    -   libgmp-dev

-   R
    -   Rcpp
    -   BH

```{sh}
sudo apt-get install libgmp-dev
```

When all dependencies are installed one can install marathonR:

```{r}
install.packages("marathonR")
```

#### windows

On windows the package is delivered with precompiled binaries. Therefore no
other dependencies are required.

```{r}
install.packages("marathonR")
```