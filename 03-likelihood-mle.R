source("_common.R")

data(aircondit, package = "boot")
y_ac <- aircondit$hours
n_ac <- length(y_ac)
y_ac
c(n = n_ac, mean = mean(y_ac), median = median(y_ac), sd = sd(y_ac))

y_di <- as.numeric(discoveries)
n_di <- length(y_di)
obs <- table(factor(y_di, levels = 0:12)) / n_di
df <- data.frame(k = 0:12, obs = as.numeric(obs), fit = dpois(0:12, mean(y_di)))
ggplot(df, aes(k)) +
  geom_col(aes(y = obs), fill = ie_blue, alpha = .6) +
  geom_point(aes(y = fit), colour = ie_coral, size = 2) +
  geom_line(aes(y = fit), colour = ie_coral) +
  labs(x = "discoveries in a year", y = "relative frequency")
c(mean = mean(y_di), variance = var(y_di))

n <- 10; y <- 8
theta <- seq(0.01, 0.99, by = 0.005)
rel <- dbinom(y, n, theta) / dbinom(y, n, y / n)
ggplot(data.frame(theta, rel), aes(theta, rel)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_hline(yintercept = exp(-qchisq(0.95, 1) / 2), linetype = 3) +
  geom_vline(xintercept = y / n, linetype = 2) +
  labs(x = expression(theta), y = "relative likelihood")
integrate(function(t) dbinom(y, n, t), 0, 1)$value   # 1/(n+1), not 1
dbinom(y, n, 0.5) / dbinom(y, n, 0.8)                # theta=0.5 vs theta=0.8

mu_hat <- mean(y_ac)
se_mu  <- mu_hat / sqrt(n_ac)
loglik_exp <- function(mu) -n_ac * log(mu) - sum(y_ac) / mu
c(mu_hat = mu_hat, se = se_mu,
  numerical = optimize(loglik_exp, c(1, 1000), maximum = TRUE)$maximum)

y_nh <- as.numeric(nhtemp); n_nh <- length(y_nh)   # New Haven mean annual temperature
mu_nh <- mean(y_nh); v_nh <- mean((y_nh - mu_nh)^2)
c(n = n_nh, mu_hat = mu_nh, v_hat = v_nh, s2 = var(y_nh))
c(formula = -n_nh / 2 * log(2 * pi * v_nh) - n_nh / 2,
  dnorm_sum = sum(dnorm(y_nh, mu_nh, sqrt(v_nh), log = TRUE)))

obs_info <- function(t) y / t^2 + (n - y) / (1 - t)^2
exp_info <- function(t) n / (t * (1 - t))
tt <- seq(0.3, 0.95, by = 0.005)
dd <- data.frame(theta = rep(tt, 2), info = c(obs_info(tt), exp_info(tt)),
                 type = rep(c("observed", "expected"), each = length(tt)))
ggplot(dd, aes(theta, info, linetype = type, colour = type)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = c(ie_coral, ie_navy)) + scale_y_log10() +
  geom_vline(xintercept = 0.8, linetype = 3) +
  labs(x = expression(theta), y = "information (log scale)", colour = NULL, linetype = NULL)
c(obs = obs_info(0.6), exp = exp_info(0.6), obs_at_mle = obs_info(0.8), se = 1 / sqrt(obs_info(0.8)))

lam_hat <- 1 / mu_hat
num <- optimize(function(l) loglik_exp(1 / l), c(1e-4, 0.1), maximum = TRUE, tol = 1e-10)$maximum
c(rate_hat = lam_hat, numerical_in_rate = num)
# MLE of P(Y > 100) = exp(-100 / mu) by invariance
exp(-100 / mu_hat)

set.seed(2027)
lam <- 3.1; nn <- 25
est <- replicate(5000, mean(rpois(nn, lam)))
c(mean_of_mle = mean(est), sd_of_mle = sd(est), asymptotic_sd = sqrt(lam / nn))
ggplot(data.frame(est), aes(est)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = ie_blue, alpha = .6) +
  stat_function(fun = dnorm, args = list(mean = lam, sd = sqrt(lam / nn)), colour = ie_coral, linewidth = 1) +
  labs(x = expression(hat(lambda)), y = "density")

I_nh <- diag(c(n_nh / v_nh, n_nh / (2 * v_nh^2)))
sqrt(diag(solve(I_nh)))

library(survival)
t_ov <- ovarian$futime;  delta_ov <- ovarian$fustat
d_ov <- sum(delta_ov);  T_ov <- sum(t_ov)
ll_cens <- function(l) d_ov * log(l) - l * T_ov
mle_ov <- d_ov / T_ov
c(d = d_ov, total_time = T_ov, mle = mle_ov, se = sqrt(d_ov) / T_ov,
  naive_1_over_mean = 1 / mean(t_ov))
# Wald 95% interval and LR test of H0: lambda = 1/1000 (mean survival 1000 days)
mle_ov + c(-1, 1) * 1.96 * sqrt(d_ov) / T_ov
LR_ov <- 2 * (ll_cens(mle_ov) - ll_cens(1 / 1000));  c(LR = LR_ov, p = pchisq(LR_ov, 1, lower.tail = FALSE))
# check against survreg (exponential model): same maximised log-likelihood
fit_ov <- survreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "exponential")
c(rate = unname(exp(-coef(fit_ov))), loglik = as.numeric(logLik(fit_ov)), ours = ll_cens(mle_ov))

yc <- c(-2, -1, 0, 7, 8)
llc <- function(t) -sapply(t, function(s) sum(log(1 + (yc - s)^2)))
scc <- function(t) sapply(t, function(s) sum(2 * (yc - s) / (1 + (yc - s)^2)))
tt <- seq(-5, 12, by = 0.02)
ggplot(data.frame(tt, ll = llc(tt)), aes(tt, ll)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_rug(data = data.frame(yc), aes(x = yc, y = NULL), sides = "b") +
  labs(x = expression(theta), y = "log-likelihood (up to a constant)")
# all roots of the score: sign changes on a fine grid
sg <- scc(tt); roots <- sapply(which(diff(sign(sg)) != 0), function(k) uniroot(scc, tt[k + 0:1])$root)
data.frame(root = roots, loglik = llc(roots), second_deriv = (llc(roots + 1e-3) - 2 * llc(roots) + llc(roots - 1e-3)) / 1e-6)

loglik_pois <- function(lambda) sum(dpois(y_di, lambda, log = TRUE))
opt <- optimize(loglik_pois, c(0.1, 20), maximum = TRUE)
c(numerical = opt$maximum, closed_form = mean(y_di))

# Two parameters: normal, parametrised on the log-scale of the variance
nll <- function(p) -sum(dnorm(y_nh, p[1], exp(p[2] / 2), log = TRUE))   # p = (mu, log v)
fit <- optim(c(50, 0), nll, method = "BFGS", hessian = TRUE)
c(mu = fit$par[1], v = exp(fit$par[2]), loglik = -fit$value)
sqrt(diag(solve(fit$hessian)))     # se(mu_hat), se(log v_hat)

nll_nb <- function(p) -sum(dnbinom(y_di, mu = exp(p[1]), size = exp(p[2]), log = TRUE))
k0 <- mean(y_di)^2 / (var(y_di) - mean(y_di))   # moment estimate of the size as starting value
fnb <- optim(c(log(mean(y_di)), log(k0)), nll_nb, method = "BFGS", hessian = TRUE)
est <- c(mu = exp(fnb$par[1]), size = exp(fnb$par[2]))
se_log <- sqrt(diag(solve(fnb$hessian)))
rbind(estimate = est, se_logscale = se_log)
c(loglik_nb = -fnb$value, loglik_pois = loglik_pois(mean(y_di)))
