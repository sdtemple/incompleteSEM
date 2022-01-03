# Evaluating simulation studies
# Seth Temple, sdtemple@uw.edu

# inputs
beta <- 1.5
gamma <- 1
N <- 400

# loading
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-vec.rds', sep = '')
vec <- readRDS(fn)
fn <- paste('../data/ii-', beta, '-', gamma, '-', N, '-mat.rds', sep = '')
mat <- readRDS(fn)

# table
bg <- c(beta, gamma)
avg_abs_dev <- function(mat, bgr, idx){return(apply(apply(mat[,idx,,drop=F], 1, function(x){abs(x - bgr[idx])}), 1, mean))}
table <- avg_abs_dev(mat, bg, 1)
table <- rbind(table, avg_abs_dev(mat, bg, 2))
table <- t(table)
colnames(table) <- c(expression(beta), expression(gamma))
#rownames(table) <- c('complete', 'partial', 'i25', 'i50', 'i75')
rownames(table) <- c('complete', 'partial', 'i20', 'i40', 'i60', 'i80')

# plot(vec, mat[,2,1], pch = 20, ylim = c(0, 3 * gamma))
# points(vec, mat[,2,2], pch = 20, col = 'red')
# points(vec, mat[,2,5], pch = 20, col = 'blue')
# abline(h=gamma)
