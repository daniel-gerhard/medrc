medrc
=====

The medrc package moved to https://github.com/DoseResponse/medrc


An R package for mixed effect dose-response modeling, combining the packages nlme and drc.

The package provides functions for fitting hierarchical dose-response models, i.e. 5-parameter log-logistic models, and automated inference for derived parameters, like effective doses, relative potency, and benchmark doses.
By numerical integration of the random effect distribution, conditional on the estimated variance components, population-averaged estimates are available.

The package can easily be installed, using the devtools package:
```
library(devtools)
install_github("daniel-gerhard/medrc")
```

