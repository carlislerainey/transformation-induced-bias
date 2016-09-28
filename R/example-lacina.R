
# clear workspace
rm(list = ls())

# install packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(grid)
library(gridExtra)

# load and clean data
lacina <- haven::read_dta("data/lacina.dta") %>%
	dplyr::select(lnbdb, lnduration, lnpop, lnmilqual, lngdp, cw, lnmountain, 
				 	democ, ethnicpolar, relpolar) %>%
	na.omit()


# simulations
simulate <- function(n_sims = 100000) {
	f <- lnbdb ~ lnduration + lnpop + lnmilqual + lngdp + cw + lnmountain + 
		democ + ethnicpolar + relpolar
	m <- lm(f, lacina, x = TRUE)
	sigma_hat <- summary(m)$sigma
	lacina_hi <- lacina; lacina_hi$democ <- 1
	m$x_hi <- m$x; m$x_hi[, "democ"] <- 1
	lacina_lo <- lacina; lacina_lo$democ <- 0
	m$x_lo <- m$x; m$x_lo[, "democ"] <- 0
	lacina$ev_log <- predict(m, newdata = lacina)
	lacina$ev <- exp(predict(m, newdata = lacina) + (sigma_hat^2)/2)
	lacina$fd <- exp(predict(m, newdata = lacina_hi) + (sigma_hat^2)/2) - 
		exp(predict(m, newdata = lacina_lo) + (sigma_hat^2)/2)
	coefs <- matrix(NA, nrow = n_sims, ncol = length(coef(m)))
	sigmas <- numeric(n_sims)
	evs <- fds <- matrix(NA, nrow = n_sims, ncol = nrow(lacina))
	# set progress bar
	pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
	# set covariates and probabilty
	for (i in 1:n_sims) {
		# simulate outcome variable
		lacina$sim_y <- rnorm(n = nrow(lacina), 
													mean = lacina$ev_log, 
													sd = sigma_hat)
		ls_fit <- lm(sim_y ~ lnduration + lnpop + lnmilqual + lngdp + cw + lnmountain + 
									 	democ + ethnicpolar + relpolar, data = lacina)
		coefs[i, ] <- coef(ls_fit)
		sigmas[i] <- summary(ls_fit)$sigma
		evs[i, ] <- exp(predict(ls_fit, newdata = lacina) + (sigmas[i]^2)/2)
		fds[i, ] <- exp(predict(ls_fit, newdata = lacina_hi) + (sigmas[i]^2)/2) - 
			exp(predict(ls_fit, newdata = lacina_lo) + (sigmas[i]^2)/2)
		setTxtProgressBar(pb, i)
	}
	e_coef <- apply(coefs, 2, mean)
	e_sigma <- mean(sigmas)
	e_ev <- apply(evs, 2, mean)
	e_fd <- apply(fds, 2, mean)
  # tau-bias in expected value
	tau_bias_ev <- e_ev - lacina$ev
	ci_tau_bias_ev <- exp(m$x%*%e_coef + (e_sigma^2)/2) - lacina$ev
	ti_tau_bias_ev <- e_ev - exp(m$x%*%e_coef + (e_sigma^2)/2)
	# tau-bias in first difference
	tau_bias_fd <- e_fd - lacina$fd
	dif <- (exp(m$x_hi%*%e_coef + (e_sigma^2)/2) - 
						exp(m$x_lo%*%e_coef + (e_sigma^2)/2))
	ci_tau_bias_fd <- dif - lacina$fd
	ti_tau_bias_fd <- e_fd - dif

	res <- list(coefs = coefs,
							evs = evs,
							fds = fds,
							e_coef = e_coef,
							e_sigma = e_sigma, 
							e_ev = e_ev,
							e_fd = e_fd,
							true_coef = coef(m),
							true_ev = lacina$ev,
							true_fd = lacina$fd,
							tau_bias_ev = tau_bias_ev,
							ci_tau_bias_ev = ci_tau_bias_ev,
							ti_tau_bias_ev = ti_tau_bias_ev,
							tau_bias_fd = tau_bias_fd,
							ci_tau_bias_fd = ci_tau_bias_fd,
							ti_tau_bias_fd = ti_tau_bias_fd)
	return(res)
}


sim <- simulate()


df1 <- data.frame(true = sim$true_ev, 
                  qi = "Expected~Value",
									bias = sim$tau_bias_ev, 
									type = "'Total'~tau*'-Bias'")
df2 <- data.frame(true = sim$true_ev,  
                  qi = "Expected~Value",
									bias = sim$ci_tau_bias_ev, 
									type = "'Coefficient-Induced'~tau*'-Bias'")
df3 <- data.frame(true = sim$true_ev,  
                  qi = "Expected~Value",
									bias = sim$ti_tau_bias_ev, 
									type = "'Transformation-Induced'~tau*'-Bias'")
df4 <- data.frame(true = sim$true_fd, 
									qi = "First~Difference",
									bias = sim$tau_bias_fd, 
									type = "'Total'~tau*'-Bias'")
df5 <- data.frame(true = sim$true_fd,  
									qi = "First~Difference",
									bias = sim$ci_tau_bias_fd, 
									type = "'Coefficient-Induced'~tau*'-Bias'")
df6 <- data.frame(true = sim$true_fd,  
									qi = "First~Difference",
									bias = sim$ti_tau_bias_fd, 
									type = "'Transformation-Induced'~tau*'-Bias'")

# df <- rbind(df1, df2, df3,
# 						df4, df5, df6)

df <- rbind(df3, df6)

# write data
saveRDS(df, "data/mc-sims/lacina-sims.RData")
