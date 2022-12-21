# Group SCR
Group SCR sampler from Emmett et al. (2021) implemented in Nimble.

https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3576

This version uses single data augmentation (for groups) rather than double data augmentation (for groups and individuals). Individual data augmentation not necessary without individual covariates and/or random effects. 

See "testscript GroupSCR.R"" and "testscript GroupSCR Marginal.R". Marginal version (marginalize over group site visits) is slower, but likely most efficient due to improved mixing.

Note, some custom functions share a name between these two samplers, so you must restart the R session when switching samplers so that Nimble will reload the correct versions. 

Update 12/21/22: There was a bug in the custom count update where I was simulating from a Poisson instead of zero-truncated Poisson. I've fixed this and started formal testing (finally). Also, I think I'll make a binary version.

Disclaimer: This code has not been tested extensively. I will remove this disclaimer when the code has been tested.
