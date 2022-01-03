# Visualize sim study results
# Seth Temple, sdtemple@uw.edu

# Inputs ------------------------------------------------------------------

flavor <- 'vanilla'
beta <- 3.2
gamma <- 1
bg <- c(beta, gamma)
N <- 300

file <- paste(flavor, beta, gamma, N, 'vec', sep = '-')
file <- paste(file, '.rds', sep = '')
vec <- readRDS(file)

file <- paste(flavor, beta, gamma, N, 'mat', sep = '-')
file <- paste(file, '.rds', sep = '')
mat <- readRDS(file)

# Graphing Parameters -----------------------------------------------------

idx <- 1
val <- bg[idx]
xlim <- c(val - val/2, val + val/2)

lty <- 1
lwd <- 2
if(idx == 1){
  xlab <- expression(beta)
} else{
  xlab <- expression(gamma)
}
ylim <- NULL
main <- NA

# colors
acol <- rcartocolor::carto_pal(6, "Magenta")
acol <- c(acol[6], acol[1], acol[2], acol[3], acol[4], acol[5])

# legend
x <- 'topright'
title <- '% obs.'

# Graphing ----------------------------------------------------------------

plot(density(mat[,idx,1]), 
     col = acol[1], 
     lwd = lwd, 
     lty = lty, 
     ylim = ylim,
     xlim = xlim,
     xlab = xlab, 
     main = main)
for(i in 2:6){lines(density(mat[,idx,i]), col = acol[i], lwd = 2, lty = 1)}
legend(x, 
       c('100','80','60','40','20','0'), 
       lwd = rep(2,6), 
       lty = rep(1,6), 
       col = acol[c(1,6:3,2)], 
       title = title)
abline(v = val, 
       col = 'black', 
       lwd = lwd, 
       lty = lty)

