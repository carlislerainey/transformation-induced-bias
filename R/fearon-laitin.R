
# install packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)

# load data
fl <- foreign::read.dta("data/fearon-laitin.dta")

# something weird is going on
table(fl$onset)  # wtf?
fl$onset[fl$onset == 4] <- 1 # recode weird case

# clean data
fl <- select(fl, country, year, onset, warl, gdpenl, lpopl1, 
                  lmtnest, ncontig, Oil, nwstate, 
                  instab, polity2l, ethfrac, relfrac) %>%
  na.omit() %>%
  mutate(country_year = paste(country, year, sep = " "))



# simulations
simulate <- function(n_sims = 100) {
  f <- onset ~ warl + gdpenl + lpopl1 + 
    lmtnest + ncontig + Oil + nwstate + 
    instab + polity2l + ethfrac + relfrac
	m <- glm(f, fl, family = binomial(link = probit), x = TRUE)
	m$x_hi <- m$x
	m$x_hi[, "gdpenl"] <- m$x[, "gdpenl"] + sd(m$x[, "gdpenl"])
	fl$prob <- predict(m, newdata = fl, type = "response")
	fl$me <- dnorm(m$x%*%coef(m))*coef(m)[3]
	fl$rr <- pnorm(m$x_hi%*%coef(m))/pnorm(m$x%*%coef(m))
	coefs <- matrix(NA, nrow = n_sims, ncol = length(coef(m)))
	probs <- mes <- rrs <- matrix(NA, nrow = n_sims, ncol = nrow(fl))
	# set progress bar
	pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
	# set covariates and probabilty
	for (i in 1:n_sims) {
		# simulate outcome variable
		fl$sim_y <- rbinom(nrow(fl), 1, prob = fl$prob)
		mle_fit <- glm(sim_y ~ warl + gdpenl + lpopl1 + 
		                 lmtnest + ncontig + Oil + nwstate + 
		                 instab + polity2l + ethfrac + relfrac, 
									 family = "binomial", data = fl)
		# mle
		# mle_fit <- glm(y ~ X - 1, family = "binomial") already done above!
		coefs[i, ] <- coef(mle_fit)
		probs[i, ] <- predict(mle_fit, newdata = fl, type = "response")
		mes[i, ] <- dnorm(m$x%*%coef(mle_fit))*coef(mle_fit)[6]
		rrs[i, ] <- pnorm(m$x_hi%*%coef(mle_fit))/pnorm(m$x%*%coef(mle_fit))
		setTxtProgressBar(pb, i)
	}
	e_coef <- apply(coefs, 2, mean)
	e_prob <- apply(probs, 2, mean)
	e_me <- apply(mes, 2, mean)
	e_rr <- apply(rrs, 2, mean)
	# tau-bias in predicted probability
	tau_bias_prob <- e_prob - fl$prob
	ci_tau_bias_prob <- pnorm(m$x%*%e_coef) - fl$prob
	ti_tau_bias_prob <- e_prob - pnorm(m$x%*%e_coef)
	# tau-bias in marginal effects
	tau_bias_me <- e_me - fl$me
	ci_tau_bias_me <- dnorm(m$x%*%e_coef)*e_coef[3] - fl$me
	ti_tau_bias_me <- e_me -  dnorm(m$x%*%e_coef)*e_coef[3]
	# tau-bias in risk ratios
	tau_bias_rr <- e_rr - fl$rr
	ci_tau_bias_rr <- pnorm(m$x_hi%*%e_coef)/pnorm(m$x%*%e_coef) - fl$rr
	ti_tau_bias_rr <- e_rr - pnorm(m$x_hi%*%e_coef)/pnorm(m$x%*%e_coef)

	res <- list(coefs = coefs,
							probs = probs,
							mes = mes,
							e_coef = e_coef,
							e_prob = e_prob, 
							e_me = e_me,
							true_coef = coef(m),
							true_prob = fl$prob,
							true_me = fl$me,
							true_rr = fl$rr,
							tau_bias_prob = tau_bias_prob,
							ci_tau_bias_prob = ci_tau_bias_prob,
							ti_tau_bias_prob = ti_tau_bias_prob,
							tau_bias_me = tau_bias_me,
							ci_tau_bias_me = ci_tau_bias_me,
							ti_tau_bias_me = ti_tau_bias_me,
							tau_bias_rr = tau_bias_rr,
							ci_tau_bias_rr = ci_tau_bias_rr,
							ti_tau_bias_rr = ti_tau_bias_rr)
	return(res)
}


sim <- simulate()


df1 <- data.frame(true = sim$true_prob, 
                  qi = "Predicted~Probability",
									bias = sim$tau_bias_prob, 
									type = "'Total'~tau*'-Bias'")
df2 <- data.frame(true = sim$true_prob,  
                  qi = "Predicted~Probability",
									bias = sim$ci_tau_bias_prob, 
									type = "'Coefficient-Induced'~tau*'-Bias'")
df3 <- data.frame(true = sim$true_prob,  
                  qi = "Predicted~Probability",
									bias = sim$ti_tau_bias_prob, 
									type = "'Transformation-Induced'~tau*'-Bias'")
df4 <- data.frame(true = sim$true_me, 
                  qi = "Marginal~Effect",
                  bias = sim$tau_bias_me, 
                  type = "'Total'~tau*'-Bias'")
df5 <- data.frame(true = sim$true_me,  
                  qi = "Marginal~Effect",
                  bias = sim$ci_tau_bias_me, 
                  type = "'Coefficient-Induced'~tau*'-Bias'")
df6 <- data.frame(true = sim$true_me,  
                  qi = "Marginal~Effect",
                  bias = sim$ti_tau_bias_me, 
                  type = "'Transformation-Induced'~tau*'-Bias'")
df7 <- data.frame(true = sim$true_rr, 
                  qi = "Risk~Ratio",
                  bias = sim$tau_bias_rr, 
                  type = "'Total'~tau*'-Bias'")
df8 <- data.frame(true = sim$true_rr,  
                  qi = "Risk~Ratio",
                  bias = sim$ci_tau_bias_rr, 
                  type = "'Coefficient-Induced'~tau*'-Bias'")
df9 <- data.frame(true = sim$true_rr,  
                  qi = "Risk~Ratio",
                  bias = sim$ti_tau_bias_rr, 
                  type = "'Transformation-Induced'~tau*'-Bias'")
df <- rbind(cbind(fl, df1), cbind(fl, df2), cbind(fl, df3))#, 
            #df4, df5, df6,
            #df7, df8, df9)

ggplot(df, aes(x = true, y = bias)) + 
	geom_point() + 
	facet_grid(qi ~ type, labeller = "label_parsed") + 
  geom_text_repel(data = subset(df, bias > .02 | bias < 0.02 | true > 0.18),
            aes(x = true, y = bias, label = country_year))




