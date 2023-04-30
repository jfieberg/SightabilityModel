# SightabilityModel 1.5.3

* Substantial Improvement in computational time for Wong.est() function.

# SightabilityModel 1.5.2

* Resubmission to remove Imports/Suggests: tidyverse
* Typographical errors corrected
* No new functionality.

# SightabilityModel 1.5.1

* Fixup vignette index to sort properly
* Misc small changes to vignettes
* No new functionality added.

# SightabilityModel 1.5

* Added new vignette dealing with Domain stratification (Heard et al 2008)
* Fixed MoosePopR() and SightabilityPopR() to deal with 100% census
* Added new functions
  + compute.SCF() to compute sightability function given data and model
  + MoosePopR_DomStrat() for domain stratification. Note that SE have not been
  adjusted for multiple domains measured on the same sampling unit. Bootstrapping may 
  be required.
  + SightabilityPopR_DomStrat() with similar functionality for sightability models
  + MoosePopR_DomStrat_bootrep() generates a bootstrap replicate and analysis. Allows
  for adjustments for sightability. Can be used for sightability models (mean-per-unit estimates) as well by
  setting block area to 1 and total block area to number of blocks.
* Numerous typographical corrections in code and vignettes

# SightabilityModel 1.4.2

* Updated URLs to Journal of Statistical Software to DOI as recommended by JSS editors.
* No change in functionality



