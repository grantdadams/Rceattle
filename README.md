# Rceattle

Rceattle: an R package for estimation of CEATTLE using template model builder



CEATTLE is short for Climate-Enhanced, Age-based model with Temperature-specific Trophic Linkages and Energetics, which is a multi-species age-structured assessment model developed for groundfish in the Bering Sea, USA by Holsman et al. (2015). To incorporate the impacts of climate, the model includes temperature-dependent von Bertalanffy weight-at-age functions (VBGF) and temperature-specific, bioenergetics-based predation interactions. Inputs of the model include U.S. National Marine Fisheries Service Alaska Fisheries Science Center (AFSC) survey and fishery data. Outputs include historical estimates of predation mortality, fishing mortality, biomass, recruitment, etc.



'Rceattle' is an 'R' package designed to implement the CEATTLE model using Template Model Builder ('TMB'; Kristensen et al. 2015), which can be installed using following https://github.com/kaskr/adcomp/wiki/Download. Rceattle is structured similar to the original manuscript in terms of modularization. Seperate functions (i.e. modules) estimate retrospective temperature- and size-specific predator rations, prey preference, and weight-at- age. These are then used as inputs to the CEATTLE model to evaluate how predation mortality, recruitment, and survival of three target species change under historical climate conditions and harvest rates.


References

Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2016). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep Sea Research Part II: Topical Studies in Oceanography, 134, 360-378.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. (2015). TMB: automatic differentiation and Laplace approximation. arXiv preprint arXiv:1509.00660.

Adams G., Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E. Dissertation. University of Washington.
