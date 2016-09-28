
# set seed
set.seed(43782)

# set parameters
n <- 10
n.sims <- 100000
x0 <- 20
b_cons <- 2.5
b_educ <- 0.1
educ <- round(seq(10, 20, length.out = 10)) #round(runif(n, 10, 21))

# do simulation
q.hat <- b0.hat <- b1.hat <- numeric(n.sims)
pb <- txtProgressBar(min = 0, max = n.sims, style = 3)
for (i in 1:n.sims) {
  #print(i)
  e <- rnorm(n)
  inc <- exp(b_cons + b_educ*educ + e)
  m <- lm(log(inc) ~ educ)
  b0.hat[i] <- coef(m)[1]
  b1.hat[i] <- coef(m)[2]
  q.hat[i] <- exp(coef(m)[1] + coef(m)[2]*x0)
  setTxtProgressBar(pb, i)
}

# check that relationship looks reasonable
plot(educ, inc)

# calculate monte carlo statistics
par(mfrow = c(1, 4))
#hist(b0.hat)
#hist(b1.hat)
#hist(q.hat)
true.q <- exp(b_cons + b_educ*x0)
true.q
mean(b0.hat)
mean(b1.hat)
mean(q.hat)

# calculate approximate bias
eXb <- exp(2.5 + 20*0.1)
0.5*(eXb*var(b0.hat) + 400*eXb*var(b1.hat) + 2*20*eXb*cov(b0.hat, b1.hat))
