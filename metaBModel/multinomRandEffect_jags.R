setwd('~/Dropbox/hawaiiDimensions/mol2ecol/metaBModel')
library(rjags)

# number of orders
nOrd <- 10

# number of species
nspp <- nOrd * 8

# vector identifying which species are in which orders
ordID <- rep(1:nOrd, each = nspp/nOrd)

# simulate random effect
ordSD <- 0.1
set.seed(1)
ordEffect <- exp(rnorm(nOrd, mean = 0, sd = ordSD))

# generate true abundances for each species
set.seed(1)
x <- sample(seq(1, 80, length.out = nspp))

# simulate number of reads for each species
totReads <- 10^6
set.seed(1)
xreads <- rmultinom(1, totReads, x*ordEffect[ordID])[, 1]

# define model
cat('model{
    
    ## reads model
    xreads[1:S] ~ dmulti(alpha[1:S], totReads)
    
    # priors
    sigX ~ dunif(0, 100)
    tauX <- 1/(sigX^2)
    sigCopy ~ dunif(0, 100)
    tauCopy <- 1/(sigCopy^2)
    
    # species-level effect
    for(i in 1:S) {
    logX[i] ~ dnorm(N/S, tauX)
    x[i] <- exp(logX[i])
    }
    # x[1:S] ~ ddirch(p0[1:S]) #dmulti(p0[1:S], N) # p0 defined in constants

    # order-level random effect
    for(j in 1:nOrd) {
    logNu[j] ~ dnorm(0, tauCopy)
    nu[j] <- exp(logNu[j])
    }
    
    # define multinom param
    for(i in 1:S) {
    alpha[i] <- x[i] * nu[ordID[i]]
    }
    
    }',
    file = 'multiRandEffect.jag' )

# model data and constants
modDat <- list(S = nspp, N = sum(x), 
               totReads = totReads, ordID = ordID, nOrd = nOrd,
               # p0 = rep(sum(x)/nspp, nspp), 
               xreads = xreads)

## compile model
mod <- jags.model(file = 'multiRandEffect.jag',
                  data = modDat,
                  n.chains = 1,
                  n.adapt = 100)

# run mcmc
thin <- 10
burn <- 1000
iter <- 5000
sampJAGS <- coda.samples(mod,
                         var = c('sigCopy', 'x', 'logNu'),
                         n.iter = (iter + burn) * thin,
                         thin = thin)
sampJAGS <- as.matrix(sampJAGS[[1]])[-(1:burn), ]

# plot abundance
ncol <- max(ordID) - 3
palette(c(hsv(h = seq(0, 0.8, length.out = ncol), 
            s = 1-0.3*seq(-1, 1, length.out = ncol)^2, 
            v = 0.7 + 0.3*seq(-1, 1, length.out = ncol)^2), 
          'black', 'white', 'gray'))


xest <- sum(x) * sampJAGS[, grep('x', colnames(sampJAGS))] / 
    rowSums(sampJAGS[, grep('x', colnames(sampJAGS))])
xestCI <- apply(xest, 2, quantile, prob = c(0.025, 0.975))
plot(x, colMeans(xest), cex = 0.5,
     panel.first = {
         segments(x0 = x, y0 = xestCI[1, ], y1 = xestCI[2, ])
         points(x, colMeans(xest), col = 'white', pch = 16, cex = 0.5)
     }, 
     ylim = range(xestCI))
abline(0, 1, col = 'red')


plot(sampJAGS[, 'sigCopy'], type = 'l')
abline(h = ordSD, col = 'red')


