library(igraph)
library(socorro)

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

par(mar = c(3.25, 3.25, 0, 0) + 0.5, mgp = c(1.5, 0, 0), xpd = FALSE)
plot(1, type = 'n', xlim = mm, ylim = mm, 
     xaxs = 'i', yaxs = 'i', axes = FALSE, frame.plot = TRUE,
     xlab = '', ylab = '')
mtext('Demogrphic or phylogenetic departure\nfrom equilibrium', side = 1, line = 2.5)
mtext('Deviation from ahistorical\necological theory', side = 2, line = 1.75)
axisArrows(lwd = 2, length = 0.1, offset = 0.125)
mtext('less', side = 1, line = 0, at = par('usr')[1] + 0.05*diff(range(par('usr')[1:2])))
mtext('more', side = 1, line = 0, at = par('usr')[2] - 0.05*diff(range(par('usr')[1:2])))
mtext('less', side = 2, line = 0.25, at = par('usr')[3] + 0.05*diff(range(par('usr')[3:4])))
mtext('more', side = 2, line = 0.25, at = par('usr')[4] - 0.05*diff(range(par('usr')[3:4])))

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

nlineage <- 10
y1 <- seq(1, 10, length.out = nlineage)
y2 <- rep(3.25, nlineage)
y3 <- seq(8, 3, length.out = nlineage)

y1sd <- abs(rnorm(nlineage, qpois(0.35, y1) - y1, sd = 0.1))
y2sd <- abs(rnorm(nlineage, 0))
y3sd <- abs(rnorm(nlineage, qpois(0.35, y3) - y3, sd = 0.1))

x <- seq(1, max(y1 + y1sd), length.out = nlineage)

col1 <- rgb(0.1, 0.9, 0.9)
col2 <- rgb(0, 0.7, 0.3)
col3 <- rgb(0.5, 0.25, 0.25)

par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(1, 0.75, 0), xpd = FALSE)
with(data.frame(x = jitter(x)), 
     plot(x, y1, type = 'b', col = col1, pch = 16, 
          xlab = 'Lineage age', ylab = 'Lineage abundance', 
          axes = FALSE, frame.plot = TRUE, cex = 1.5, asp = 1, 
          ylim = range(y1 - y1sd, y1 + y1sd), 
          panel.first = {
              segments(x0 = x, y0 = y1-y1sd, y1 = y1+y1sd, col = col1)
          }))

with(data.frame(x = jitter(x)), {
    segments(x0 = x, y0 = y2-y2sd, y1 = y2+y2sd, col = col2)
    points(x, y2, col = col2, type = 'b', pch = 16, cex = 1.5)
})

with(data.frame(x = jitter(x)), {
    segments(x0 = x, y0 = y3-y3sd, y1 = y3+y3sd, col = col3)
    points(x, seq(8, 3, length.out = 10), col = col3, 
           type = 'b', pch = 16, cex = 1.5)
})

polygon(x = 8 + c(0, cos(pi*40/180) * strwidth('evo-eco drift'), 
                  cos(pi*40/180) * strwidth('evo-eco drift') + cos(pi*120/180) * 
                      strheight('evo-eco drift'), 
                  cos(pi*120/180) * strheight('evo-eco drift')), 
        y = 8 + c(0, sin(pi*40/180) * strwidth('evo-eco drift'), 
                  sin(pi*40/180) * strwidth('evo-eco drift') + sin(pi*120/180) * 
                      strheight('evo-eco drift'), 
                  sin(pi*120/180) * strheight('evo-eco drift')), 
        col = 'white', border = 'white', lwd = 3)

text('evo-eco drift', x = 8, y = 8, srt = 40, col = rgb(0.1, 0.9, 0.9), adj = c(0, 0))

rect(xleft = 4.5, xright = 4.5 + strwidth('eco disturbance'), 
     ybottom = 2.25, ytop = 2.25 + strheight('eco disturbance'), 
     col = 'white', border = 'white', lwd = 3)
text('eco disturbance', x = 4.5, y = 2.25, adj = c(0, 0), col = rgb(0, 0.7, 0.3))

polygon(x = 1.5 + c(0, cos(-pi*25/180) * strwidth('evo innovation'), 
                    cos(-pi*25/180) * strwidth('evo innovation') + cos(pi*65/180) * 
                        strheight('evo innovation'), 
                    cos(pi*65/180) * strheight('evo innovation')), 
        y = 8.5 + c(0, sin(-pi*25/180) * strwidth('evo innovation'), 
                    sin(-pi*25/180) * strwidth('evo innovation') + sin(pi*65/180) * 
                        strheight('evo innovation'), 
                    sin(pi*65/180) * strheight('evo innovation')), 
        col = 'white', border = 'white', lwd = 3)
text('evo innovation', x = 1.5, y = 8.5, adj = c(0, 0), col = rgb(0.5, 0.25, 0.25), 
     srt = -25)

dev.off()
