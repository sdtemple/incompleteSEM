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

lty <- 1
lwd <- 2
xlab <- expression(R[0])
ylim <- c(0,5)
val <- beta / gamma; xlim <- c(val - val / 3, val + val / 3)
main <- NA

# colors
acol <- rcartocolor::carto_pal(6, "Magenta")
acol <- c(acol[6], acol[1], acol[2], acol[3], acol[4], acol[5])

# legend
x <- 'topright'
title <- '% obs.'

# Graphing ----------------------------------------------------------------

z <- mat[,1,1] / mat[,2,1]
plot(density(z), 
     col = acol[1], 
     lwd = lwd, 
     lty = lty, 
     ylim = ylim,
     xlim = xlim,
     xlab = xlab, 
     main = main)
for(i in 2:6){
  z <- mat[,1,i] / mat[,2,i]
  lines(density(z), 
        col = acol[i], 
        lwd = 2, 
        lty = 1)
}
legend(x, 
       c('100','80','60','40','20','0'), 
       lwd = rep(2,6), 
       lty = rep(1,6), 
       col = acol[c(1,6:3,2)], 
       title = title)
abline(v = beta / gamma, 
       col = 'black', 
       lwd = lwd, 
       lty = lty)

