# Shared setup for thesis scripts.
# Loads packages, sets reproducibility seed, and defines shared paths.

library(readxl)
library(quantmod)
library(ggplot2)
library(ggthemes)
library(forecast)
library(MASS)
library(rugarch)
library(rmgarch)
library(reshape2)
library(fGarch)
library(evir)
library(xts)
library(tseries)
library(gridExtra)
library(zoo)
library(copula)
library(matrixcalc)
library(pracma)
library(MCMCpack)
library(viridis)

set.seed(191023)
bond_data_path <- "data/raw/final_bond_data.xlsx"
