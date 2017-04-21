
knitr::opts_chunk$set(echo = TRUE)
library(nimble)


# number of orders
nOrd <- 7

# number of species
nspp <- nOrd * 10

# vector identifying which species are in which orders
ordID <- rep(1:nOrd, each = nspp/nOrd)

# simulate random effect
ordSD <- 0.5
set.seed(1)
ordEffect <- exp(rnorm(nOrd, mean = 0, sd = ordSD))

# generate true abundances for each species
set.seed(1)
x <- sample(1:40, nspp, rep = TRUE)

# simulate number of reads for each species
totReads <- 10^6
set.seed(1)
nreads <- rmulti(1, totReads, x*ordEffect[ordID])


# palette(hsv(h = seq(0, 0.8, length.out = max(ordID)), 
#             s = 1-0.3*seq(-1, 1, length.out = max(ordID))^2, 
#             v = 0.7 + 0.3*seq(-1, 1, length.out = max(ordID))^2, 
#             alpha = 0.7))
# plot(x, nreads, bg = ordID, pch = 21, 
#      xlab = 'Abundance', ylab = 'Number of reads')


modCode <- nimbleCode({
    ## reads model
    nreads[1:S] ~ dmulti(alpha[1:S], totReads)
    
    # priors
    x[1:S] ~ dmulti(p0[1:S], N) # p0 defined in constants 
    tauCopy ~ dgamma(0.01, 0.0001)
    
    # order-level random effect
    for(j in 1:nOrd) {
        logNu[j] ~ dnorm(0, tauCopy)
        nu[j] <- exp(logNu[j])
    }
    
    # define multinom param
    for(i in 1:S) {
        alpha[i] <- x[i] * nu[ordID[i]]
    }
})

# model constants, data and inits
modConstants <- list(S = nspp, N = sum(x), totReads = totReads, ordID = ordID, nOrd = nOrd,
                     p0 = rep(1/nspp, nspp))

modData <- list(nreads = nreads)

modInits <- list( 
    # priors
    tauCopy = 1,
    x = as.numeric(table(cut(1:modConstants$N, nspp))), # uniform partition
    
    # hyperdistribution
    logNu = rep(0, modConstants$nOrd),
    nu = rep(1, modConstants$nOrd),
    
    # deterministic relationships
    alpha = rep(1/nspp, nspp)
)

# build model
mod <- nimbleModel(code = modCode, 
                   constants = modConstants, data = modData, inits = modInits)


Cmod <- compileNimble(mod)
modConf <- configureMCMC(mod)


modConf$getMonitors()
modConf$addMonitors('logNu')


# compile MCMC
modMCMC <- buildMCMC(modConf)
CmodMCMC <- compileNimble(modMCMC)

# set MCMC iterations
mcmcN <- 5e+03
burn <- 4e+02
niter <- (mcmcN + burn) * modConf$thin
CmodMCMC$run(niter)

# the posterior sample    
samp <- as.matrix(CmodMCMC$mvSamples)[-(1:burn), ]
write.csv(samp, file = '~/Dropbox/hawaiiDimensions/mol2ecol/mcmcRes_multinomRand.csv', 
          row.names = FALSE)

# # load mcmc output from background run
# samp <- as.matrix(read.csv('mcmcRes_multinomRand.csv', header = TRUE))
# colnames(samp) <- strsplit(gsub('"', '', 
#                                 readLines('mcmcRes_multinomRand.csv', n = 1)), 
#                            ',')[[1]]


# par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
# plot(x, colMeans(samp[, grep('x', colnames(samp))]), bg = ordID, pch = 21, 
#      xlab = 'Actual Abundance', ylab = 'Estimated')
# abline(0, 1)


# plot(round(seq(1, nrow(samp), length.out = nrow(samp)/10)), 
#      samp[round(seq(1, nrow(samp), length.out = nrow(samp)/10)), 'tauCopy']^-0.5, 
#      type = 'l', ylim = range(ordSD, samp[, 'tauCopy']^-0.5), 
#      xlab = 'Iterations', ylab = expression(sigma[order]))
# abline(h = ordSD, col = 'red')


# plot(ordEffect, colMeans(samp[, grep('Nu', colnames(samp))]), bg = unique(ordID), pch = 21, 
#      xlab = 'Actual Order Effect', ylab = 'Estimated Order Effect', 
#      xlim = range(ordEffect, samp[, grep('Nu', colnames(samp))]), 
#      ylim = range(ordEffect, samp[, grep('Nu', colnames(samp))]))
# abline(0, 1)
print('done')
