## simple word count function for short, uncluttered text

wc <- function(x) {
    x <- strsplit(x, ' |\n')[[1]]
    x <- x[x != '']
    length(x)
}

## function to count how many times each citation is used
library(stringr)

ncite <- function(tex) {
    keys <- read.csv('~/Dropbox/hawaiiDimensions/mol2ecol/keyHash.csv', as.is = TRUE)
    keys <- unique(keys$newKey)
    
    x <- readLines(tex)
    
    out <- sapply(keys, function(s) sum(str_count(x, s)))
    
    return(sort(out[out > 0]))
}

write.csv(cbind(ncite('ms/overleaf/main_clean.tex')), 'ncite.csv')
