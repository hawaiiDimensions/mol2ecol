setwd('~/Dropbox/hawaiiDimensions/mol2ecol/ms')

x <- readLines('mol2ecol_dirty.tex')

## clean up citations and references

citeLines <- grep('\\textbackslash{}cite\\{', x, fixed = TRUE)
x[citeLines] <- gsub('\\textbackslash{}cite\\{', '\\cite{', x[citeLines], fixed = TRUE)
x[citeLines] <- gsub('\\}', '}', x[citeLines], fixed = TRUE)

refLines <- grep('\\textbackslash{}ref\\{', x, fixed = TRUE)
x[refLines] <- gsub('\\textbackslash{}ref\\{', '\\ref{', x[refLines], fixed = TRUE)
x[refLines] <- gsub('\\}', '}', x[refLines], fixed = TRUE)


## remove all labels
x <- gsub('\\\\label.*', '', x)

## remove quotes
x <- x[-grep('begin\\{quote\\}|end\\{quote\\}', x)]

## remove 'def'
x <- x[-grep('\\\\def', x)]

## add appropriate header
x <- c(readLines('header.tex'), x[-(1:grep('\\begin{document}', x, fixed = TRUE))])

writeLines(x, 'mol2ecol.tex')
