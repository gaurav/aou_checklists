library(rstan)
library(tidyverse)
library(reshape2)

# Simulating data and evaluating parameter recovery --------------------------

# Generate indices for taxonomic groups
nobs <- 834
norder <- 25
nfamily <- 80
ngenus <- 338

order <- sample(norder, nobs, replace = TRUE)
family <- sample(nfamily, nobs, replace = TRUE)
genus <- sample(ngenus, nobs, replace = TRUE)

# generate offset term
years_in_list <- sample(1:50, nobs, replace = TRUE)

# simulate true parameter values from their prior distributions
lambda_0 <- rnorm(1)

sigma_i <- abs(rnorm(1))
sigma_j <- abs(rnorm(1))
sigma_k <- abs(rnorm(1))

pi_i <- rnorm(norder, 0, sigma_i)
tau_j <- rnorm(nfamily, 0, sigma_j)
rho_k <- rnorm(ngenus, 0, sigma_k)

# store these for later
true_ranefs <- list(pi_i = pi_i, tau_j = tau_j, rho_k = rho_k)


log_lambda <- lambda_0 + 
  pi_i[order] + 
  tau_j[family] + 
  rho_k[genus] + 
  log(years_in_list)

# simulate data from prior predictive distribution
stopifnot(length(log_lambda) == nobs)

y <- rpois(nobs, exp(log_lambda))

# evaluate prior predictive distribution for y
# note that this is one draw from the PPD
hist(y)

stan_d <- list(nobs = nobs, 
               norder = norder, 
               nfamily = nfamily, 
               ngenus = ngenus, 
               y = y, 
               order = order, 
               family = family, 
               genus = genus, 
               offset = log(years_in_list))

# Fit model and evaluate parameter recovery -------------------------------
m_fit <- stan("counts_per_name_model.stan", data = stan_d, 
              control = list(max_treedepth = 12))

post <- rstan::extract(m_fit)

compare_ranefs <- function(which_ranef = "pi_i") {
  post[[which_ranef]] %>%
    melt(varnames = c("iter", "idx")) %>%
    tbl_df %>%
    group_by(idx) %>%
    summarize(lo = quantile(value, .025), 
              med = median(value), 
              hi = quantile(value, .975)) %>%
    mutate(true = true_ranefs[[which_ranef]][idx]) %>%
    ggplot(aes(true, med)) + 
    geom_point() + 
    geom_segment(aes(xend = true, y = lo, yend = hi)) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    xlab("True value") + 
    ylab("Estimated value") + 
    ggtitle(paste("Recovery of", which_ranef))
}

compare_ranefs("pi_i")
compare_ranefs("tau_j")
compare_ranefs("rho_k") # some genera are never observed (super wide CIs)

# evaluate recovery of remaining params
traceplot(m_fit, "lambda_0") + 
  geom_hline(aes(yintercept = lambda_0))

traceplot(m_fit, "sigma_i") + 
  geom_hline(aes(yintercept = sigma_i))

traceplot(m_fit, "sigma_j") + 
  geom_hline(aes(yintercept = sigma_j))

traceplot(m_fit, "sigma_k") + 
  geom_hline(aes(yintercept = sigma_k))

