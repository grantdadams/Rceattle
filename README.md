# Rceattle

Rceattle: an R package for estimation of CEATTLE using template model builder. 

NOTE: this branch is the most up to date (including stock-recruit environmental linkages), but has some model outputs commented out for quicker MSE runs.


CEATTLE is short for Climate-Enhanced, Age-based model with Temperature-specific Trophic Linkages and Energetics, which is a multi-species age-structured assessment model developed for groundfish in the Bering Sea, USA by Holsman et al. (2015) and Gulf of Alaska, USA (Adams et al, 2022). Essentially, it is a multispecies statistical catch-at-age model. To incorporate the impacts of climate, the model includes temperature-specific, bioenergetics-based predation interactions. Inputs of the model include survey/fishery index, catch, and length/age composition data in addition to empirical weight-at-age, diet proportion, and ration data. Outputs include historical estimates of predation mortality, fishing mortality, biomass, recruitment, etc.



'Rceattle' is an 'R' package designed to implement the CEATTLE model using Template Model Builder ('TMB'; Kristensen et al. 2015), which can be installed using following https://github.com/kaskr/adcomp/wiki/Download. Rceattle is structured similar to Adams et al (2022). Data are read in via an excel document (see examples) for model fitting (see examples). Projections can be conducted under alternative harvest control rules, climate projections, and recruitment. Model diagnostic, validation, simulation, and closed loop simulation testing (management strategy evaluation) functions are included as well. The package supports one- or -two sex models with multiple fisheries and surveys with flexible catchability and selectivity parameterizations. See vignette (in progress) for model parameterizations. 

Code and function examples using data from the Bering Sea and Gulf of Alaska groundfish application can be found in the [examples](https://github.com/grantdadams/Rceattle/tree/master/examples) folder.

**Getting started**

See [onboarding document](https://github.com/grantdadams/Rceattle/wiki/Onboarding) and Wiki for getting started!

**References**

Holsman, K. K., Ianelli, J., Aydin, K., Punt, A. E., & Moffitt, E. A. (2016). A comparison of fisheries biological reference points estimated from temperature-specific multi-species and single-species climate-enhanced stock assessment models. Deep Sea Research Part II: Topical Studies in Oceanography, 134, 360-378.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. (2015). TMB: automatic differentiation and Laplace approximation. arXiv preprint arXiv:1509.00660.

Adams, G.D., Holsman, K.K., Barbeaux, S.J., Dorn, M.W., Ianelli, J.N., Spies, I., Stewart, I.J., Punt, A.E., 2022. An ensemble approach to understand predation mortality for groundfish in the Gulf of Alaska. Fish. Res. 251, 106303. doi:10.1016/j.fishres.2022.106303
