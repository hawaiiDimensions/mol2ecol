setwd('~/Dropbox/hawaiiDimensions/mol2ecol')

cleanBib <- function(bib) {
    x <- readLines(bib)
    
    dubSpace <- TRUE
    while(dubSpace) {
        x <- gsub('  ', ' ', x)
        dubSpace <- any(grepl('  ', x))
    }
    
    x <- x[x != '']
    x <- x[substring(x, 1, 1) != '%']
    
    entryStart <- grep('@.*\\{', x)
    entryStart <- entryStart[substring(x[entryStart], 1, 1) == '@']
    
    entries <- split(x, cut(1:length(x), c(1, entryStart[-1], length(x)+1) - 1))
    names(entries) <- NULL
    
    cleanEntries <- lapply(entries, function(e) {
        fieldStart <- grep('=', e)
        fieldStart <- fieldStart[unlist(gregexpr('=', e[fieldStart])) <= 16]
        fieldStart <- fieldStart[!grepl('S= cAz', e[fieldStart])]
        
        fields <- split(e, cut(1:length(e), c(1, fieldStart, length(x)+1) -1))
        names(fields) <- NULL
        
        ## clean up bib key
        if(grepl('-', fields[[1]][1])) {
            ss <- c(0, 2) + max(unlist(gregexpr('-', fields[[1]][1])))
            fields[[1]][1] <- paste0(substring(fields[[1]][1], 1, ss[1]-1), 
                                     substring(fields[[1]][1], ss[2]+1))
            fields[[1]][1] <- tolower(fields[[1]][1])
        }
        
        bad <- sapply(fields, function(f) {
            any(grepl('abstract =|affiliation =|month =|language =|keywords =', f, 
                      ignore.case = TRUE))
        })
        fields <- unlist(fields[!bad])
        
        if(fields[length(fields)] != '}') {
            if(substring(fields[length(fields)], nchar(fields[length(fields)])) == ',') 
                fields[length(fields)] <- substring(fields[length(fields)], 
                                                    1, nchar(fields[length(fields)]) - 1)
            fields <- c(fields, '}')
        }
        
        fields <- c(fields, '')
        
        return(fields)
    })
    
    keyHash <- lapply(1:length(entries), function(i) {
        keys <- c(entries[[i]][1], cleanEntries[[i]][1])
        keys <- gsub('.*\\{|,.*', '', keys)
        return(keys)
    })
    
    keyHash <- do.call(rbind, keyHash)
    
    dup <- duplicated(keyHash[, 2])
    dup <- keyHash[, 2] %in% keyHash[dup, 2]
    
    addThis <- sapply(which(dup), function(i) {
        e <- cleanEntries[[i]]
        s <- tolower(substring(gsub('journal = |"| |\\.|,', '', 
                                    e[grep('journal =', e)]), 1, 9))
        
        cleanEntries[[i]][1] <<- gsub(',', paste0(s, ','), cleanEntries[[i]][1])
        return(s)
    })
    
    keyHash[dup, 2] <- paste0(keyHash[dup, 2], addThis)
    
    keyHash <- cbind(keyHash, c('', 'check')[as.numeric(dup) + 1])
    colnames(keyHash) <- c('oldKey', 'newKey', 'check')
    keyHash <- keyHash[order(keyHash[, 2]), ]
    
    writeLines(unlist(cleanEntries), gsub('.bib', '_clean.bib', bib))
    
    return(keyHash)
}

keyHash <- cleanBib('ms/overleaf/references.bib')
write.csv(keyHash, file = 'keyHash.csv', row.names = FALSE)

## GO CHECK THE FILES AND COME BACK!


## fix inline citations

foo <- readLines('ms/overleaf/main.tex')
keyHash <- read.csv('keyHash.csv', as.is = TRUE)

for(i in 1:nrow(keyHash)) {
    foo <- gsub(keyHash$oldKey[i], keyHash$newKey[i], foo)
}

writeLines(foo, 'ms/overleaf/main_clean.tex')
