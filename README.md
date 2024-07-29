# jdcbioinfo
Miscellaneous R functions helpful for analysis. Not built for public consumption.

<!-- badges: start -->
[![R-CMD-check](https://github.com/jdreyf/jdcbioinfo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jdreyf/jdcbioinfo/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/jdcbioinfo/master.svg)](https://codecov.io/github/jdreyf/jdcbioinfo?branch=master)
<!-- badges: end -->

## Install
library(remotes)  
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"))
remotes::install_github(repo="jdreyf/ezlimmaplot", build_opts = c("--no-resave-data", "--no-manual"))
remotes::install_github(repo="jdreyf/Hitman", build_opts = c("--no-resave-data", "--no-manual"))
remotes::install_github(repo="jdreyf/jdcbioinfo", build_opts = c("--no-resave-data", "--no-manual"))