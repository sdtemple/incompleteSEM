library(incompleteSEM)

N <- 100
beta <- 4
gamma <- 1
epi <- rgsem(beta, gamma, N)

epi <- caseri(epi)
epi <- sortri(epi)
out <- mle_gsem(epi[,2], epi[,1], N)
print("completely observed")
c(out[1] / out[2], out)
complete_likelihood(epi[,2], epi[,1], beta, gamma, N)
epi2 <- missri(epi, 9/10, 1/2)
epi2 <- sortri(epi2)
r <- epi2[,2]
i <- epi2[,1]

pbla_incomplete(r, i, beta, gamma, N)
pbla_std_gsem(epi[,2], beta, gamma, N)
complete_likelihood(epi[,2], epi[,1], beta, gamma, N)

# two step
# out <- nlm(pbla_two_step, c(1), r, i, N)$estimate
# g <- mle_removal_rate(r, i)
# print('two step')
# unname(c(out / g[1], out, g[1]))

# pbla universal
out <- nlm(pbla_uni_wrap, c(1,1), r=r, i=i, N=N)$estimate
print('universal')
c(out[1] / out[2], out)

# pbla standard
out <- nlm(pbla_gsem, c(1,1), pbla=pbla_std_gsem, r=epi[,2], N=N)$estimate
print('partially observed')
c(out[1] / out[2], out)

#pbla_uni(r, i, 1.5, 1, N)
#pbla_std_gsem(r, 1.5, 1, N)
#epi2

# for(l in 1:10){
#   epi2 <- caseri(epi)
#   epi2 <- missri(epi2, 1, 1)
#   epi2 <- sortri(epi2)
#   r <- epi2[,2]
#   i <- epi2[,1]  
#   print(c(nlm(pbla_two_step, c(1), r=r, i=i, N=100)$estimate, mle_removal_rate(r, i)[1]))
#   print(nlm(pbla_uni_wrap, c(1,1), r=r, i=i, N=100)$estimate)
#   print(nlm(pbla_gsem, c(1,1), pbla=pbla_std_gsem, r=r, N=100)$estimate)
# }

