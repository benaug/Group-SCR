# Group SCR
Group SCR sampler from Emmett et al. (2021) implemented in Nimble. Assumes individual and pack ID is known. Can possibly be extended to allow for unknown pack IDs.

https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3576

These nimble samplers use single data augmentation (for groups) rather than double data augmentation (for groups and individuals) as in the original paper. Individual data augmentation not necessary without individual covariates and/or random effects. 

There are 4 samplers here. Two assume group site visitation and individual detection given group site visitation are count variables. Two others assume they are binary variables. Two other possibilities are that one is a count variable and the other is binary. Versions with binary site visitation run the most quickly and were not provided in Emmett et al. 

For each of the sampler types described above, there are two versions, one that updates the latent group site use variables and a second that marginalizes over them. The marginal versions will typically be preferred, I think. They use less RAM, mix better, and at least for the binary samplers, run faster.

See "testscripts" for each model type. "testscript GroupSCR ..." are count versions, and "testscript GroupSCR Bernoulli" are binary versions. Marginal versions indicated with "marginal" in file name.

All 4 versions have been tested via simulation with the parameter values and sampling design indicated in the testscripts. The original R-based samplers that use most of the same MCMC methods have been tested much more extensively.

Note, some custom functions share a name between these two samplers, so you must restart the R session when switching samplers so that Nimble will reload the correct versions or deregister them.

Final note: This model can be very slow due to the group site visitation latent variables and becomes progressively slower as the number of groups, individuals, traps, and occasions increase. The occasion-level data is necessary to fit the model, so we cannot aggregate detections over occasions and work with the "i x j" data. This is why the testscripts are set up with K=2. K=1 could also be handled, but requires changes to custom updates, etc., because nimble won't let you use an array with one dimension size being 1. E.g., M x J x K works, M x J x 1 does not work as an array, you need to tell nimble it is an M x J matrix.