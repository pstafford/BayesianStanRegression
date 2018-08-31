# BayesianStanRegression

This repository provides the `stan` files used in the article:

Stafford, P.J. (2018), Continuous integration of data into ground-motion models using Bayesian updating. _Journal of Seismology_, _under review_

The repository contains just two `*.stan` files:
- `PJSstanModel_reES_AS.stan` is the file for performing  regression analysis on a traditional dataset. The data passed to this code contains observations from a number of earthquake events, and the code specifies informative univariate priors for the fixed effects. Note that the priors are defined for the peak ground acceleration as used in the above article.
- `PJSstanModel_reES_update_AS.stan` is the file used for performing Bayesian updating. The data passed to this code contains just observations from a single new event. The priors are specified as a multivariate prior for the fixed effects and univariate priors for the variance components and random effects.

Note that in both cases the file name contains the string `_reES_` to indicate that random effects for _Event_ and _Station_ are considered in these scripts.

To actually run the regressions, ground-motion data must be prepared and passed through calls to the Stan program. While the `rstan` interface was used in this particular study, it is possible to prepare the data and interact with `Stan` in a number of ways. See http://mc-stan.org for more details.
