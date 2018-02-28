# pestim
Population estimation using mobile phone data

`pestim` R package implements a hierarchical model to estimate the population counts of different
territorial cells combining the information from aggregated mobile phone data and a population
register or survey data, both at a given time instant and along a sequence of time periods.

`pestim` package is freely available under the GPL3 and EUPL licenses at the following address:
[https://github.com/MobilePhoneESSnetBigData/pestim]. It requires at least R version
3.3.0, but upgrading R to the newest version is highly recommended. It can be installed
using `install_github()` function from devtools package:
```r
library(devtools)
install_github("MobilePhoneESSnetBigData/pestim", build_vignettes=TRUE)
```

The Reference Manual is available at the following address: [https://github.com/MobilePhoneESSnetBigData/pestim/raw/master/doc/pestim_Reference_Manual.pdf].
