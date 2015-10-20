
approximate_tau_bias <- function(beta_hat, tau_fn, Sigma, ...) {
	# checks
	if (class(beta_hat) != "numeric") {
		warning("The argument `beta_hat` must be a numeric vector.")
	}
	if (class(tau_fn) != "function") {
		warning("The argument `tau` must be a function.")
	}
	if (class(Sigma) != "matrix") {
		warning("The argument `Sigma` must matrix.")
	}
	if (nrow(Sigma) != ncol(Sigma)) {
		warning("Sigma is not a square matrix.")
	}
	if (nrow(Sigma) != length(beta_hat)) {
		warning("Dimension of Sigma does not match length of beta_hat.")
	}
	H <- hessian(tau_fn, beta_hat)
  approx_tau_bais <- 0.5*sum(H*Sigma)
  return(approx_tau_bais)
}

tau_fn <- function(beta_hat) {
	tau <- exp(beta_hat[1] + beta_hat[2]*20)
}

# bootstrap quantities of interest
bs_fn <- function(data, index) {
	bs_data <- data[index, ]
	bs_m <- lm(log(inc) ~ educ, data = bs_data)
	q_hat <- tau_fn(coef(bs_m))
	return(q_hat)
}

n <- 10
x <- rnorm(n)
b_cons <- 2.5
b_educ <- 0.1
educ <- round(runif(n, 10, 21))
x0 <- 20


n_sims <- 1000
q_hat <- approx_tau_bias <- bs_bias <- numeric(n_sims)
for (i in 1:n_sims) {
	e <- rnorm(n)
	inc <- exp(b_cons + b_educ*educ + e)
	m <- lm(log(inc) ~ educ)
	q_hat[i] <- tau_fn(coef(m))
	df <- data.frame(inc = inc, educ = educ)
	bs <- boot::boot(df, bs_fn, R = 1000, parallel = "multicore", ncpus = 4)
	bs_bias[i] <- apply(bs$t, 2, mean) - bs$t0
	approx_tau_bias[i] <- approximate_tau_bias(coef(m), tau_fn = tau_fn, Sigma = vcov(m))
}
true_tau <- tau_fn(c(b_cons, b_educ))
e_tau_hat <- mean(q_hat)
true_bias <- e_tau_hat - true_tau
true_bias
mean(approx_tau_bias)
mean(bs_bias, na.rm = TRUE)

approximate_tau_bias(coef(m), tau_fn = tau_fn, Sigma = vcov(m))


bs <- boot::boot(df, bs_fn, R = 2000, parallel = "multicore", ncpus = 4)

bs_bias <- apply(bs$t, 2, mean) - bs$t0
true_bias
bs_bias
approximate_tau_bias(coef(m), tau_fn = tau_fn, Sigma = vcov(m))

