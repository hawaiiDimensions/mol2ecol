setwd('~/Dropbox/hawaiiDimensions/mol2ecol')
source('simMetaBar.R')
library(nimble)
library(socorro)

# set parameters for simulation
S <- 70
set.seed(1)
n <- sample(round(seq(1, 40, length.out = S)))
N <- sum(n)
Nreads <- 1e+06
sigCopy <- 0.1

# simulate phylogeny and number of reads resulting form metabarcoding
set.seed(1)
sim <- simMetaBar(abund = n, sigCopy = sigCopy, nreads = Nreads)

# extract needed objects from output
numberReads <- sim$reads
tre <- sim$tre
ordID <- sim$ordID

# model code
mod <- nimbleCode({
    # phylo correlated copy number
    PrecCopy[1:nOrd, 1:nOrd] <- tauCopy * Dinv[1:nOrd, 1:nOrd]
    logNu[1:nOrd] ~ dmnorm(mu0[1:nOrd], PrecCopy[1:nOrd, 1:nOrd])
    nu[1:nOrd] <- exp(logNu[1:nOrd])
    
    # define dirichlet-multinom params and the distribution of x_{reads}
    for(i in 1:S) {
        alpha[i] <- x[i] * nu[ordID[i]]
    }
    
    ## reads model
    xreads[1:S] ~ dmulti(alpha[1:S], Nreads)
    
    # priors
    x[1:S] ~ dmulti(p0[1:S], N) # p0 defined in constants 
    tauCopy ~ dgamma(0.01, 0.0001)
})

# model constants, data and inits
modConstants <- list(S = S, N = N, Nreads = Nreads, ordID = ordID, nOrd = max(ordID),
                     Dinv = solve(vcv(tre)), # vcv gives the var-cov mat of a phylo for sig=1
                     mu0 = rep(0, max(ordID)), p0 = rep(1/S, S))

modData <- list(xreads = numberReads)

modInits <- list( 
    # priors
    tauCopy = 1,
    x = {
        foo <- rep(round(N / S), S)
        foo[0:(N - sum(foo))] <- foo[0:(N - sum(foo))] + 1 
        foo
    },
    
    # deterministic relationships
    PrecCopy = modConstants$Dinv,
    nu = rep(1, modConstants$nOrd),
    
    # hyperdistribution
    logNu = rep(0, modConstants$nOrd)
)

# build model
mod <- nimbleModel(code = mod, 
                   constants = modConstants, data = modData, inits = modInits)

# compile model
Cmod <- compileNimble(mod)
modConf <- configureMCMC(mod)

# add monitor for logNu
modConf$addMonitors('logNu')

# set the thinning param
thin <- 10
modConf$setThin(thin)

# add block sampling for x
modConf$addSampler(target = sprintf('logNu[1:%s]', modConstants$nOrd),
                   type = 'RW_block')

# compile MCMC
modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC, project = mod)

# set MCMC params
mcmcN <- 5e+03
burn <- 4e+02
niter <- (mcmcN + burn) * modConf$thin
CmodMCMC$run(niter)

# the posterior sample    
samp <- as.matrix(CmodMCMC$mvSamples)[-(1:burn), ]
apply(exp(samp[, grep('Nu', colnames(samp))]), 2, range)

# write out figures
pdf('fig_realVsEstimatedAbund.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(sim$abund, colMeans(samp[, grep('x', colnames(samp))]), 
     xlab = 'Real abundance', ylab = 'Estimated abundance', 
     pch = 16, ylim = range(samp[, grep('x', colnames(samp))]), 
     panel.first = segments(x0 = sim$abund, 
                            y0 = apply(samp[, grep('x', colnames(samp))], 2, 
                                       quantile, prob = 0.025), 
                            y1 = apply(samp[, grep('x', colnames(samp))], 2, 
                                       quantile, prob = 0.975)))
abline(0, 1, col = 'red')
dev.off()

pdf('fig_sigCopy.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(1/samp[, 'tauCopy'], type = 'l', xlab = 'Iterations', ylab = expression(sigma[copy]))
abline(h = sigCopy, col = 'red')
dev.off()


## save output and input
write.csv(samp, file = 'mcmcRes_noPrimerAffin.csv', row.names = FALSE)
save(sim, file = 'mcmcInput_noPrimerAffin.RData')
