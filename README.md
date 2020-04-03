# SightabilityModel
Contains code for the SightabilityModel Package

## Versions and installation

  * **Github** To install the latest development version from Github, 
    install the newest version of the **devtools** package; then run
```
devtools::install_github("jfieberg/SightabilityModel", dependencies = TRUE,
                        build_vignettes = TRUE)
```
## Features
Estimates abundance, density, and ratio of variables in simple random 
sample or stratified random sample block design often used in studying
ungulate populations.

In a typical study, the survey area is divided into strata and each
stratum divided into blocks. Not all blocks need be the same area.
Aerial surveys of selected blocks find groups of animals and the number
of animals in the group is recorded, often divided by sex (bulls and cows)
and age (calve and mature).

A key problem is that not all animals in a group are seen
because sightability is < 1. A separate study of a known
number of animal (e.g. radio collared) is used to estimate
a sightability model. The actual sightability study data or the fitted
mode from a sightability model using logistic regression is used to
correct the observed counts in the study. The corrected (for sightability)
counts are used to estimate density, abundance or ratios (e.g. bulls to cows).

## References

Fieberg, J. 2012. Estimating Population Abundance Using Sightability Models: R SightabilityModel Package. Journal of Statistical Software, 51(9), 1-20. URL http://www.jstatsoft.org/v51/i09/.

Steinhorst, Kirk R. and Samuel, Michael D. 1989. Sightability Adjustment Methods for Aerial Surveys of Wildlife Populations. Biometrics 45:415–425.

Wong, C. 1996. Population size estimation using the modified Horvitz-Thompson estimator with estimated sighting probabilities. Dissertation, Colorado State University, Fort Collins, USA.

