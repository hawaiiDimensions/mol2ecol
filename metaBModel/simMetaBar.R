library(ape)
library(TreeSim)
# library(pika)


## ========================================
## function to simulate community phylogeny
## ========================================

## assumes homogeneous birth-death tree for spp (with age of MRCA = 25), then coalescent tree 
## for population samples whoes abundances are given by `abund` and who each get assigned a 
## random Ne according to `rexp(1, 10^6)`
makeCommTree <- function(abund) {
    nspp <- length(abund)
    
    ## simulate tree
    T0 <- 25
    tre <- sim.bd.taxa.age(nspp, numbsim = 1, lambda = 1, mu = 0.9, 
                           age = T0, mrca = TRUE)[[1]]
    
    ## simulate coalescent populations at tips
    for(i in 1:length(abund)) {
        x <- abund[i]
        if(x > 1) {
            ## coalescent
            out <- rcoal(x, tip.label = paste('c', i, ':', 1:x, sep = ''))
            
            ## rescale edges by Ne (note that timescale of phylo is in MY, thus `rexp(1, 1)`)
            out$edge.length <- out$edge.length * rexp(1, 1)
            r0 <- max(node.depth.edgelength(out))
            r0max <- tre$edge.length[tre$edge[, 2] == which(tre$tip.label == 
                                                                sprintf('t%s', i))] * 0.9
            ## if coalescent is too old for that species edge, re-scale it
            if(r0 > r0max) {
                out$edge.length <- (out$edge.length / r0) * r0max
                r0 <- r0max
            }
            
            ## graft population onto spp phylogeny
            tre <- bind.tree(tre, out, where = which(tre$tip.label == sprintf('t%s', i)), 
                             position = r0)
            tre <- drop.tip(tre, which(tre$tip.label == sprintf('t%s', i)))
        }
    }
    
    return(tre)
}


## =========================================
## function to simulate metabarcoding output
## =========================================

## assumes copy number and primar affinity evolve by independent BM (log and logit transformed,
## repectively). returns resulting number of reads for each spp, as well as their actual
## abundances
# simMetaBar <- function(abund, minAffin = 0.8, sigCopy, sigPrimer, nreads = 10^6) {
simMetaBar <- function(abund, sigCopy, nreads = 10^6) {
    nspp <- length(abund)
    
    ## simulate order-level tree
    T0 <- 25
    tre <- sim.bd.taxa.age(ifelse(nspp/10 < 3, 3, round(nspp/10)), 
                           numbsim = 1, lambda = 1, mu = 0.9, 
                           age = T0, mrca = TRUE)[[1]]
    
    ## simulate copy number evolution as log-normal
    cn <- exp(rTraitCont(tre, sigma = sqrt(sigCopy)))
    
    ## assign spp to orders and replicate cn
    ordAbund <- split(abund, cut(abund, length(tre$tip.label)))
    ordID <- rep(1:length(ordAbund), sapply(ordAbund, length))
    cn <- cn[ordID]
    
    ## simulate reads as a multinomial sample given abund, cn
    reads <- rmultinom(1, nreads, abund*cn)[, 1]
    
    ## return starting abundance and resultant reads
    return(list(abund = abund, reads = reads, ordID = ordID, tre = tre, copy = cn))
}
