source("_common.R")
library(stats4)

data(aircondit, package = "boot")
y_ac <- aircondit$hours; n_ac <- length(y_ac)
y_di <- as.numeric(discoveries); n_di <- length(y_di)
w_gy <- MASS::geyser$waiting

n <- 10; y <- 8; th0 <- 0.5; thh <- y / n
ll <- function(t) y * log(t) + (n - y) * log(1 - t)
sc <- function(t) y / t - (n - y) / (1 - t)
tt <- seq(0.25, 0.97, by = 0.005)
quad <- ll(thh) - 0.5 * n / (thh * (1 - thh)) * (tt - thh)^2
ggplot(data.frame(tt, l = ll(tt), q = quad, tg = ll(th0) + sc(th0) * (tt - th0)), aes(tt)) +
  geom_line(aes(y = l), colour = ie_navy, linewidth = 1) +
  geom_line(aes(y = q), colour = ie_green, linetype = 2) +
  geom_line(aes(y = tg), colour = ie_coral, linetype = 2) +
  geom_vline(xintercept = c(th0, thh), linetype = 3) +
  coord_cartesian(ylim = c(ll(thh) - 4, ll(thh) + 0.5)) +
  annotate("text", x = 0.9, y = ll(thh) - 0.6, label = "quadratic (Wald)", colour = ie_green, hjust = 0.2) +
  annotate("text", x = 0.3, y = ll(th0) + 1.9, label = "tangent at the null", colour = ie_coral, hjust = 0) +
  labs(x = expression(theta), y = "log-likelihood")

W  <- (thh - th0)^2 * n / (thh * (1 - thh))
Sc <- sc(th0)^2 / (n / (th0 * (1 - th0)))
LR <- 2 * (ll(thh) - ll(th0))
stats <- c(Wald = W, Score = Sc, LR = LR)
round(rbind(statistic = stats, p_value = pchisq(stats, df = 1, lower.tail = FALSE)), 4)

z <- qnorm(0.975); se <- sqrt(thh * (1 - thh) / n)
wald <- thh + c(-1, 1) * z * se
score_ci <- as.numeric(prop.test(y, n, correct = FALSE)$conf.int)   # Wilson = inverted score test
cut <- qchisq(0.95, df = 1) / 2                                    # 1.92
f <- function(t) ll(t) - ll(thh) + cut
lr_ci <- c(uniroot(f, c(1e-8, thh), tol = 1e-10)$root, uniroot(f, c(thh, 1 - 1e-8), tol = 1e-10)$root)
round(rbind(Wald = wald, Score = score_ci, LR = lr_ci), 3)

mu_hat <- mean(y_ac)
llmu <- function(mu) -n_ac * log(mu) - sum(y_ac) / mu
wald_mu    <- mu_hat + c(-1, 1) * z * mu_hat / sqrt(n_ac)
wald_logmu <- exp(log(mu_hat) + c(-1, 1) * z / sqrt(n_ac))
g <- function(mu) llmu(mu) - llmu(mu_hat) + cut
lr_mu <- c(uniroot(g, c(1, mu_hat), tol = 1e-10)$root, uniroot(g, c(mu_hat, 5000), tol = 1e-10)$root)
exact <- 2 * sum(y_ac) / qchisq(c(0.975, 0.025), df = 2 * n_ac)
round(rbind(Wald_mu = wald_mu, Wald_log_mu = wald_logmu, LR = lr_mu, Exact = exact), 1)
# Invariance: the LR interval for the rate 1/mu is the reciprocal interval
rbind(rate_from_LR_mu = rev(1 / lr_mu))

mm <- seq(40, 330, by = 1)
rl <- exp(sapply(mm, llmu) - llmu(mu_hat))
ggplot(data.frame(mm, rl), aes(mm, rl)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_hline(yintercept = exp(-cut), linetype = 3) +
  annotate("segment", x = lr_mu[1], xend = lr_mu[2], y = 0.04, yend = 0.04, colour = ie_coral, linewidth = 2) +
  annotate("segment", x = wald_mu[1], xend = wald_mu[2], y = 0.01, yend = 0.01, colour = ie_green, linewidth = 2) +
  annotate("text", x = 300, y = 0.08, label = "LR", colour = ie_coral) +
  annotate("text", x = 300, y = 0.0, label = "Wald", colour = ie_green) +
  labs(x = expression(mu), y = "relative likelihood")

# Full gamma fit, parametrised on the log scale
nll_g <- function(p) -sum(dgamma(y_ac, shape = exp(p[1]), rate = exp(p[2]), log = TRUE))
fit_g <- optim(c(0, log(1 / mean(y_ac))), nll_g, method = "BFGS", hessian = TRUE)
c(alpha_hat = exp(fit_g$par[1]), beta_hat = exp(fit_g$par[2]), loglik = -fit_g$value, conv = fit_g$convergence)
ll_exp   <- sum(dexp(y_ac, rate = 1 / mean(y_ac), log = TRUE))
LR_g <- 2 * (-fit_g$value - ll_exp)
c(ll_exp = ll_exp, LR = LR_g, p_value = pchisq(LR_g, df = 1, lower.tail = FALSE))
# Wald test of alpha = 1 via the delta method: se(alpha) = alpha * se(log alpha)
se_log <- sqrt(diag(solve(fit_g$hessian)))
a_hat <- exp(fit_g$par[1]); se_a <- a_hat * se_log[1]
z_w <- (a_hat - 1) / se_a
c(se_alpha = se_a, z = z_w, p_value = 2 * pnorm(-abs(z_w)))

ll_pois <- sum(dpois(y_di, mean(y_di), log = TRUE))
nll_nb  <- function(p) -sum(dnbinom(y_di, mu = exp(p[1]), size = exp(p[2]), log = TRUE))
k0 <- mean(y_di)^2 / (var(y_di) - mean(y_di))          # moment estimate of the size, used as starting value
fit_nb  <- optim(c(log(mean(y_di)), log(k0)), nll_nb, method = "BFGS", hessian = TRUE)
LR_nb   <- 2 * (-fit_nb$value - ll_pois)
c(mu_hat = exp(fit_nb$par[1]), k_hat = exp(fit_nb$par[2]), loglik_nb = -fit_nb$value,
  loglik_pois = ll_pois, LR = LR_nb)
c(p_naive = pchisq(LR_nb, 1, lower.tail = FALSE), p_boundary = 0.5 * pchisq(LR_nb, 1, lower.tail = FALSE))
# the same model from the naive start log(k) = 0 (k = 1): BFGS runs off towards k = Inf, the Poisson limit
bad <- optim(c(log(mean(y_di)), 0), nll_nb, method = "BFGS")
c(k_hat = exp(bad$par[2]), loglik = -bad$value, convergence = bad$convergence)

y_nh <- as.numeric(nhtemp); n_nh <- length(y_nh)
lp <- function(m) -n_nh / 2 * log(mean((y_nh - m)^2))
m_hat <- mean(y_nh)
k <- function(m) lp(m) - lp(m_hat) + cut
prof_ci <- c(uniroot(k, c(m_hat - 3, m_hat), tol = 1e-10)$root, uniroot(k, c(m_hat, m_hat + 3), tol = 1e-10)$root)
s2_hat <- mean((y_nh - m_hat)^2)
round(rbind(profile_LR = prof_ci, Wald = m_hat + c(-1, 1) * z * sqrt(s2_hat / n_nh),
            t_interval = as.numeric(t.test(y_nh)$conf.int)), 3)
mm <- seq(50.6, 51.7, length.out = 200)
ggplot(data.frame(mm, p = sapply(mm, lp) - lp(m_hat)), aes(mm, p)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_hline(yintercept = -cut, linetype = 2) + geom_vline(xintercept = prof_ci, linetype = 3) +
  labs(x = expression(mu), y = "profile log-likelihood")

lp_a <- function(a) sum(dgamma(y_ac, shape = a, rate = a / mean(y_ac), log = TRUE))
gg <- function(a) lp_a(a) - lp_a(a_hat) + cut
prof_a <- c(uniroot(gg, c(0.05, a_hat), tol = 1e-10)$root, uniroot(gg, c(a_hat, 10), tol = 1e-10)$root)
I_a <- n_ac * (trigamma(a_hat) - 1 / a_hat)          # profile information
wald_a <- exp(log(a_hat) + c(-1, 1) * z / (a_hat * sqrt(I_a)))   # Wald on log(alpha)
round(rbind(profile_LR = prof_a, Wald_log_scale = wald_a), 3)
c(se_alpha_profile = 1 / sqrt(I_a), se_alpha_hessian = se_a)
aa <- seq(0.15, 2.2, length.out = 200)
ggplot(data.frame(aa, p = sapply(aa, lp_a) - lp_a(a_hat)), aes(aa, p)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_hline(yintercept = -cut, linetype = 2) + geom_vline(xintercept = c(a_hat, 1), linetype = c(3, 1), colour = c("grey40", ie_coral)) +
  coord_cartesian(ylim = c(-4, 0.2)) + labs(x = expression(alpha), y = "profile log-likelihood")

S_a <- function(a) n_ac * (log(a) - digamma(a) - log(mean(y_ac)) + mean(log(y_ac)))
I_a_fun <- function(a) n_ac * (trigamma(a) - 1 / a)
a <- 1; iter <- data.frame()
for (k in 0:5) {
  iter <- rbind(iter, data.frame(k = k, alpha = a, score = S_a(a), info = I_a_fun(a)))
  a <- a + S_a(a) / I_a_fun(a)
}
print(iter, digits = 6, row.names = FALSE)

nr_rate <- function(l0, K = 6) { l <- numeric(K + 1); l[1] <- l0; for (k in 1:K) l[k + 1] <- l[k] * (2 - mean(y_ac) * l[k]); l }
rbind(start_0.005 = nr_rate(0.005), start_0.02 = nr_rate(0.02)) |> signif(4)
1 / mean(y_ac)

yc <- c(-2, -1, 0, 7, 8); nc <- length(yc)
Sc_ <- function(t) sum(2 * (yc - t) / (1 + (yc - t)^2))
Ic_ <- function(t) sum(2 * (1 - (yc - t)^2) / (1 + (yc - t)^2)^2)
run <- function(t0, info, K = 25) { t <- t0; for (k in 1:K) t <- t + Sc_(t) / info(t); t }
starts <- c(-3, 0, 3, 4, 6, 10)
data.frame(start = starts,
           newton   = sapply(starts, function(s) run(s, Ic_)),
           fisher   = sapply(starts, function(s) run(s, function(t) nc / 2)))

fit_mle <- function(nll, start, ...) {
  fit <- optim(start, nll, method = "BFGS", hessian = TRUE, ...)
  if (fit$convergence != 0) warning("optim() did not report convergence")
  ev <- eigen(fit$hessian, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) warning("Hessian not positive definite: not a proper maximum")
  list(par = fit$par, se = sqrt(diag(solve(fit$hessian))), loglik = -fit$value, eigen = ev)
}
fit <- fit_mle(nll_g, c(0, log(1 / mean(y_ac))))
rbind(estimate_log = fit$par, se_log = fit$se)
fit$eigen

nll_exp_log <- function(logmu) -sum(dexp(y_ac, rate = exp(-logmu), log = TRUE))
m1 <- mle(nll_exp_log, start = list(logmu = 4), method = "BFGS")
coef(summary(m1))
exp(confint(m1))     # profile-likelihood interval for mu: it matches the LR interval above

n0 <- 20
c(LR_upper = 1 - exp(-cut / n0),
  exact_Clopper_Pearson = binom.test(0, n0)$conf.int[2],
  Wilson_upper = prop.test(0, n0, correct = FALSE)$conf.int[2],
  rule_of_three = 3 / n0)

ggplot(data.frame(w = w_gy), aes(w)) + geom_histogram(aes(y = after_stat(density)), bins = 30, fill = ie_blue, alpha = .6) +
  labs(x = "waiting time (min)", y = "density")

nll_mix <- function(par, y) {
  p <- plogis(par[1])
  -sum(log(p * dnorm(y, par[2], exp(par[4])) + (1 - p) * dnorm(y, par[3], exp(par[5]))))
}
show_fit <- function(fit) round(c(p = plogis(fit$par[1]), mu1 = fit$par[2], mu2 = fit$par[3],
                                  sigma1 = exp(fit$par[4]), sigma2 = exp(fit$par[5]),
                                  loglik = -fit$value, conv = fit$convergence), 3)
good <- optim(c(0, 50, 80, log(5), log(5)), nll_mix, y = w_gy, method = "BFGS", hessian = TRUE)
sym  <- optim(c(0, mean(w_gy), mean(w_gy), log(sd(w_gy)), log(sd(w_gy))), nll_mix, y = w_gy, method = "BFGS", hessian = TRUE)
rbind(sensible_start = show_fit(good), symmetric_start = show_fit(sym))
rbind(eigen_good = eigen(good$hessian)$values, eigen_sym = eigen(sym$hessian)$values)

spike <- function(s) -nll_mix(c(0, w_gy[1], 80, log(s), log(7.5)), w_gy)
sapply(c(1e-2, 1e-4, 1e-6, 1e-8), spike)

set.seed(2027)
yn <- rnorm(30, mean = 3, sd = 1)
nll_ab <- function(p) -sum(dnorm(yn, mean = p[1] + p[2], sd = 1, log = TRUE))
o1 <- optim(c(0, 0), nll_ab, method = "BFGS", hessian = TRUE)
o2 <- optim(c(5, -5), nll_ab, method = "BFGS", hessian = TRUE)
rbind(start_0_0 = c(o1$par, sum = sum(o1$par)), start_5_m5 = c(o2$par, sum = sum(o2$par)))
eigen(o1$hessian)$values

em_mix <- function(y, p, m, s, maxit = 500, tol = 1e-10) {
  ll <- numeric(0)
  for (k in 1:maxit) {
    d1 <- p * dnorm(y, m[1], s[1]); d2 <- (1 - p) * dnorm(y, m[2], s[2])
    ll[k] <- sum(log(d1 + d2))                       # observed-data log-likelihood BEFORE the update
    g <- d1 / (d1 + d2)                              # E-step: responsibilities
    p <- mean(g)                                     # M-step
    m <- c(sum(g * y) / sum(g), sum((1 - g) * y) / sum(1 - g))
    s <- sqrt(c(sum(g * (y - m[1])^2) / sum(g), sum((1 - g) * (y - m[2])^2) / sum(1 - g)))
    if (k > 1 && abs(ll[k] - ll[k - 1]) < tol) break
  }
  list(p = p, mu = m, sigma = s, ll = ll)
}
fit_em <- em_mix(w_gy, 0.5, c(50, 80), c(5, 5))
c(p = fit_em$p, mu = fit_em$mu, sigma = fit_em$sigma, loglik = tail(fit_em$ll, 1), iterations = length(fit_em$ll))
all(diff(fit_em$ll) >= -1e-9)                        # monotone increase
ggplot(data.frame(k = seq_along(fit_em$ll), ll = fit_em$ll), aes(k, ll)) + geom_line(colour = ie_navy) + geom_point(size = 1) +
  labs(x = "iteration", y = "log-likelihood")

lr_interval <- function(loglik, mle, lower, upper, level = 0.95) {
  cut <- qchisq(level, df = 1) / 2
  f <- function(t) loglik(t) - loglik(mle) + cut
  c(uniroot(f, c(lower, mle), tol = 1e-10)$root, uniroot(f, c(mle, upper), tol = 1e-10)$root)
}
lr_interval(function(l) sum(dpois(y_di, l, log = TRUE)), mean(y_di), 1, 10)   # Poisson, discoveries

s <- 50; n <- 25; lam <- s / n
ll <- function(l) s * log(l) - n * l
wald <- lam + c(-1, 1) * qnorm(0.975) * sqrt(lam / n)
lr <- lr_interval(ll, lam, 0.5, 6)
sc <- uniroot(function(l) (s / l - n)^2 / (n / l) - qchisq(.95, 1), c(0.5, lam))$root
sc2 <- uniroot(function(l) (s / l - n)^2 / (n / l) - qchisq(.95, 1), c(lam, 6))$root
round(rbind(Wald = wald, Score = c(sc, sc2), LR = lr), 3)

s <- sum(y_ac); lam <- 1 / mean(y_ac)
exact <- qchisq(c(0.025, 0.975), 2 * n_ac) / (2 * s)
wald_l <- lam + c(-1, 1) * z * lam / sqrt(n_ac)             # se(lambda_hat) = lambda_hat / sqrt(n)
wald_ll <- exp(log(lam) + c(-1, 1) * z / sqrt(n_ac))
lr_l <- rev(1 / lr_mu)                                       # invariance
signif(rbind(Exact = exact, Wald = wald_l, Wald_log = wald_ll, LR = lr_l), 3)

vh <- mean((y_nh - mean(y_nh))^2)
LRv <- function(v0) n_nh * (vh / v0 - 1 - log(vh / v0))
c(LR = LRv(1), p = pchisq(LRv(1), 1, lower.tail = FALSE))
c(uniroot(function(v) LRv(v) - qchisq(.95, 1), c(0.5, vh))$root, uniroot(function(v) LRv(v) - qchisq(.95, 1), c(vh, 5))$root)

eta <- numeric(8); eta[1] <- log(1.2)
for (k in 1:7) eta[k + 1] <- eta[k] + exp(-eta[k]) / 2 - 1
rbind(lambda = exp(eta))

nll_c <- function(t) sum(log(1 + (yc - t)^2))
res <- t(sapply(c(-3, 0, 4, 6, 10), function(s) { o <- optim(s, nll_c, method = "BFGS"); c(start = s, theta_hat = o$par, loglik = -o$value, conv = o$convergence) }))
round(res, 3)

set.seed(2027)
boot_LR <- replicate(300, {
  yb <- rexp(n_ac, rate = 1 / mean(y_ac))
  nl <- function(p) -sum(dgamma(yb, shape = exp(p[1]), rate = exp(p[2]), log = TRUE))
  o <- optim(c(0, log(1 / mean(yb))), nl, method = "BFGS")
  2 * (-o$value - sum(dexp(yb, 1 / mean(yb), log = TRUE)))
})
c(boot_p_value = mean(boot_LR >= LR_g), chisq_p = pchisq(LR_g, 1, lower.tail = FALSE))

ll1 <- sum(dnorm(w_gy, mean(w_gy), sqrt(mean((w_gy - mean(w_gy))^2)), log = TRUE))
ll2 <- tail(fit_em$ll, 1)
n <- length(w_gy)
data.frame(model = c("1 normal", "2-component mixture"), npar = c(2, 5), loglik = c(ll1, ll2),
           AIC = -2 * c(ll1, ll2) + 2 * c(2, 5), BIC = -2 * c(ll1, ll2) + log(n) * c(2, 5))

set.seed(2027)
LR_obs <- 2 * (ll2 - ll1)
LR_boot <- replicate(100, {
  yb <- rnorm(n, mean(w_gy), sd(w_gy))
  l1 <- sum(dnorm(yb, mean(yb), sqrt(mean((yb - mean(yb))^2)), log = TRUE))
  f <- em_mix(yb, 0.5, mean(yb) + c(-1, 1) * sd(yb), rep(sd(yb) / 2, 2), maxit = 300, tol = 1e-6)
  2 * (tail(f$ll, 1) - l1)
})
c(LR_obs = LR_obs, boot_max = max(LR_boot), boot_95 = unname(quantile(LR_boot, 0.95)), p_value = (1 + sum(LR_boot >= LR_obs)) / 101)
