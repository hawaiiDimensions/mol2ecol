library(ape)
library(phylosim)
library(pika)


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

bla <- makeCommTree(1:10)
plot(bla, show.tip.label = FALSE)



