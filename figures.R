library(igraph)

## figure showing cycles of equilibrium and non-equilibrium

ncol <- 50
mm <- c(-1.5, 1.5)

zcol <- expand.grid(seq(0, 1, length.out = ncol), seq(0, 1, length.out = ncol))
zgrid <- expand.grid(seq(mm[1], mm[2], length.out = ncol), 
                     seq(mm[1], mm[2], length.out = ncol))
colz <- rgb(0.4*(rowSums(zcol)/2), 0.3+0.6*(1 - zcol[, 1]), 0.3+0.6*(1 - zcol[, 2]))

ecoCycle <- make_graph(c(1, 2, 2, 1, 2, 3, 3, 4, 4, 1), directed = TRUE)
evoCycle <- make_graph(c(1, 4, 4, 3, 3, 1), directed = TRUE)

pdf('ms/overleaf/fig_cycles.pdf', width = 4, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1.5, 0, 0), xpd = FALSE)
plot(1, type = 'n', xlim = mm, ylim = mm, 
     xaxs = 'i', yaxs = 'i', axes = FALSE, frame.plot = TRUE,
     xlab = 'Demogrphic of phylogenetic departure\nfrom equilibrium', 
     ylab = 'Deviation from ahistorical\necological theory')

rect(xleft = zgrid[, 1], 
     xright = zgrid[, 1] + diff(par('usr')[1:2])/(ncol-1), 
     ybottom = zgrid[, 2], 
     ytop = zgrid[, 2] + diff(par('usr')[3:4])/(ncol-1), 
     col = colz, border = colz)

abline(h = mean(mm),v = mean(mm), col = 'gray')

plot(ecoCycle, 
     layout = matrix(c(1, 1, 1, 2, 2, 2, 2, 1), ncol = 2, byrow = TRUE), 
     vertex.size = 50, vertex.color = 'transparent', vertex.frame.color = 'transparent', 
     vertex.label = c('I', 'II', 'III', 'IV'), vertex.label.color = 'black', 
     edge.curved = c(0.3, 0.3, rep(0.2, 3)), edge.color = 'black', edge.width = 3,
     add = TRUE)

plot(evoCycle, 
     layout = matrix(c(1, 1, 1, 2, 2, 2, 2, 1), ncol = 2, byrow = TRUE), 
     vertex.size = 50, vertex.color = 'transparent', vertex.frame.color = 'transparent', 
     vertex.label = c('I', 'II', 'III', 'IV'), vertex.label.color = 'black', 
     edge.curved = c(0.2, 0.2, -0.2), edge.color = 'white', edge.width = 3,
     add = TRUE)

dev.off()


## figure showing age abundance relationship

pdf('ms/overleaf/fig_age-abund.pdf', width = 4, height = 4)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1, 0.75, 0), xpd = FALSE)
plot(1:10, type = 'b', col = colz[1], pch = 16, 
     xlab = 'Lineage age', ylab = 'Lineage abundance', 
     axes = FALSE, frame.plot = TRUE, cex = 1.5)
points(rep(3.5, 10), col = rgb(0.4, 0.8, 0.3), type = 'b', pch = 16, cex = 1.5)
points(10:1, col = colz[length(colz)], type = 'b', pch = 16, cex = 1.5)

dev.off()
