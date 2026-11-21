source("_common.R")
suppressPackageStartupMessages({library(lme4); library(MASS)})
data(cbpp, package = "lme4")

cb <- cbpp; cb$prop <- cb$incidence / cb$size; cb$period_n <- as.integer(cb$period)
pool <- aggregate(cbind(incidence, size) ~ period, cbpp, sum); pool$prop <- pool$incidence / pool$size; pool$period_n <- 1:4
ggplot(cb, aes(period_n, prop, group = herd)) + geom_line(colour = ie_navy, alpha = .4) +
  geom_point(aes(size = size), colour = ie_navy, alpha = .6) +
  geom_line(data = pool, aes(group = 1), colour = ie_coral, linewidth = 1.3) +
  scale_size_area(max_size = 4, name = "herd size") + labs(x = "period", y = "proportion of new cases")

fit_glm <- glm(cbind(incidence, size - incidence) ~ period, family = binomial, data = cbpp)
sum(residuals(fit_glm, type = "pearson")^2) / df.residual(fit_glm)

set.seed(1); B <- 1e6; sb2 <- .5
bb <- rnorm(B, 0, sqrt(sb2)); mu <- 4 * exp(bb); y1 <- rpois(B, mu); y2 <- rpois(B, mu)
m_th <- 4 * exp(sb2 / 2)
rbind(theory = c(mean = m_th, var = m_th + m_th^2 * (exp(sb2) - 1), cov = m_th^2 * (exp(sb2) - 1)),
      simulation = c(mean(y1), var(y1), cov(y1, y2)))

y0 <- 3
f0 <- function(b) dpois(y0, exp(b)) * dnorm(b)
exact <- integrate(f0, -Inf, Inf)$value
bhat <- uniroot(function(b) y0 - exp(b) - b, c(0, 3), tol = 1e-12)$root
H <- exp(bhat) + 1
lap <- f0(bhat) * sqrt(2 * pi / H)
c(mode = bhat, curvature = H, Laplace = lap, exact = exact, relative_error = lap / exact - 1)
bg <- seq(-2.5, 3.5, length.out = 300)
ggplot(data.frame(b = rep(bg, 2), f = c(f0(bg), f0(bhat) * exp(-.5 * H * (bg - bhat)^2)),
                  what = rep(c("exact integrand", "Gaussian (Laplace)"), each = length(bg))),
       aes(b, f, colour = what, linetype = what)) + geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(ie_coral, ie_navy), name = NULL) + scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  labs(y = "integrand")

fit_lap <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp)
beta <- fixef(fit_lap); sb <- attr(VarCorr(fit_lap)$herd, "stddev")
h1 <- subset(cbpp, herd == "1"); eta1 <- drop(model.matrix(~ period, cbpp)[cbpp$herd == "1", ] %*% beta)
logh <- function(b) sapply(b, function(bb) sum(dbinom(h1$incidence, h1$size, plogis(eta1 + bb), log = TRUE)) + dnorm(bb, 0, sb, log = TRUE))
exact1 <- integrate(function(b) exp(logh(b)), -Inf, Inf)$value
mode1 <- optimize(logh, c(-4, 4), maximum = TRUE)$maximum
p1 <- plogis(eta1 + mode1); Hh <- sum(h1$size * p1 * (1 - p1)) + 1 / sb^2     # Z'WZ + G^{-1}
lap1 <- unname(exp(logh(mode1)) * sqrt(2 * pi / Hh))
c(mode = mode1, exact = exact1, Laplace = lap1, ratio = lap1 / exact1)

gh_rule <- function(K) {                                     # Golub-Welsch: nodes and weights for exp(-x^2)
  J <- matrix(0, K, K); if (K > 1) J[cbind(1:(K - 1), 2:K)] <- J[cbind(2:K, 1:(K - 1))] <- sqrt((1:(K - 1)) / 2)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
quad <- function(f, bhat, H, K, adaptive = TRUE) {
  r <- gh_rule(K)
  if (adaptive) sqrt(2 / H) * sum(r$w * exp(r$x^2) * f(bhat + sqrt(2 / H) * r$x))
  else          sum(r$w * exp(r$x^2) * f(r$x))                # nodes centred at 0, scale 1
}
tab <- do.call(rbind, lapply(c(3, 30), function(yy) {
  f <- function(b) dpois(yy, exp(b)) * dnorm(b)
  bh <- uniroot(function(b) yy - exp(b) - b, c(0, 6), tol = 1e-12)$root; Hh <- exp(bh) + 1
  Ks <- c(1, 3, 5, 9, 15)
  data.frame(y = yy, K = Ks, adaptive = sapply(Ks, function(K) quad(f, bh, Hh, K)),
             plain = sapply(Ks, function(K) quad(f, bh, Hh, K, adaptive = FALSE)),
             exact = integrate(f, -Inf, Inf, rel.tol = 1e-12)$value)
}))
print(tab, digits = 4)

summary(fit_lap)$coefficients
VarCorr(fit_lap)
fit_agq <- update(fit_lap, nAGQ = 25)
pql <- MASS::glmmPQL(cbind(incidence, size - incidence) ~ period, random = ~ 1 | herd, family = binomial, data = cbpp, verbose = FALSE)
sd_of <- function(f) if (inherits(f, "lme")) as.numeric(VarCorr(f)[1, 2]) else attr(VarCorr(f)$herd, "stddev")
round(rbind(`Laplace (nAGQ = 1)` = c(fixef(fit_lap), sd = sd_of(fit_lap)),
            `AGQ (nAGQ = 25)`    = c(fixef(fit_agq), sd = sd_of(fit_agq)),
            PQL                  = c(fixef(pql), sd = sd_of(pql))), 3)

fit0 <- update(fit_lap, . ~ . - period)
anova(fit0, fit_lap)                                       # LRT for period, 3 df
fit_nore <- glm(cbind(incidence, size - incidence) ~ period, family = binomial, data = cbpp)
T_h <- as.numeric(2 * (logLik(fit_lap) - logLik(fit_nore)))
c(T = T_h, p_naive = pchisq(T_h, 1, lower.tail = FALSE), p_mixture = .5 * pchisq(T_h, 1, lower.tail = FALSE))
confint(update(fit_lap), method = "profile", parm = "beta_")        # profile intervals; done on a copy because profiling changes the internal state of a merMod object

re <- ranef(fit_lap, condVar = TRUE)$herd
round(head(cbind(mode = re[, 1], se = sqrt(attr(ranef(fit_lap, condVar = TRUE)$herd, "postVar")[1, 1, ])), 5), 3)

set.seed(2027)
m <- 2000; n_i <- 10; id <- rep(1:m, each = n_i); x <- rnorm(m * n_i); b <- rnorm(m, 0, 2)[id]
d <- data.frame(y = rbinom(m * n_i, 1, plogis(-1 + 1 * x + b)), x = x, id = factor(id))
fit_cond <- glmer(y ~ x + (1 | id), family = binomial, data = d)
fit_pop  <- glm(y ~ x, family = binomial, data = d)
round(rbind(true_conditional = c(-1, 1), glmer_conditional = fixef(fit_cond), glm_marginal = coef(fit_pop),
            approximation = c(-1, 1) / sqrt(1 + 0.346 * 4)), 3)
xs <- seq(-3, 3, by = .05)
p_marg <- sapply(xs, function(x0) integrate(function(bb) plogis(-1 + x0 + bb) * dnorm(bb, 0, 2), -Inf, Inf)$value)
coef(lm(qlogis(p_marg) ~ xs))                  # the exact marginal curve is close to a logistic one

xs <- seq(-4, 4, by = .05); bq <- qnorm(c(.1, .25, .75, .9), 0, 2)
cond <- expand.grid(x = xs, b = bq); cond$p <- plogis(-1 + cond$x + cond$b)
marg <- data.frame(x = xs, p = sapply(xs, function(x0) integrate(function(bb) plogis(-1 + x0 + bb) * dnorm(bb, 0, 2), -Inf, Inf)$value))
ggplot() + geom_line(data = cond, aes(x, p, group = b), colour = "grey60") +
  geom_line(data = data.frame(x = xs, p = plogis(-1 + xs)), aes(x, p, colour = "typical cluster"), linewidth = 1.2) +
  geom_line(data = marg, aes(x, p, colour = "population average"), linewidth = 1.2) +
  geom_line(data = data.frame(x = xs, p = plogis((-1 + xs) / sqrt(1 + .346 * 4))), aes(x, p, colour = "approximation"), linetype = 2, linewidth = 1) +
  scale_colour_manual(values = c("typical cluster" = ie_green, "population average" = ie_coral, "approximation" = ie_navy), name = NULL) +
  labs(y = "P(y = 1 | x)")

eta <- fixef(fit_lap)[1] + c(0, fixef(fit_lap)[-1])
rbind(typical_herd = plogis(eta),
      population_average = sapply(eta, function(e) integrate(function(b) plogis(e + b) * dnorm(b, 0, sb), -Inf, Inf)$value),
      observed = tapply(cbpp$incidence, cbpp$period, sum) / tapply(cbpp$size, cbpp$period, sum))

data(bacteria, package = "MASS")
fit_b <- glmer(y ~ trt + I(week > 2) + (1 | ID), family = binomial, data = bacteria)
sb_b <- attr(VarCorr(fit_b)$ID, "stddev")
cbind(GLMM_conditional = fixef(fit_b), approx_marginal = fixef(fit_b) / sqrt(1 + .346 * sb_b^2),
      glm_marginal = coef(glm(y ~ trt + I(week > 2), family = binomial, data = bacteria)))

gt <- grouseticks; gt$HEIGHT_s <- as.numeric(scale(gt$HEIGHT))          # standardised altitude
c(mean = mean(gt$TICKS), variance = var(gt$TICKS))
pdisp <- function(fit) sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
fit_g0 <- glm(TICKS ~ YEAR + HEIGHT_s, family = poisson, data = gt)
fit_g1 <- glmer(TICKS ~ YEAR + HEIGHT_s + (1 | LOCATION) + (1 | BROOD), family = poisson, data = gt)
fit_g2 <- update(fit_g1, . ~ . + (1 | INDEX))                              # INDEX = row number: the OLRE
tab_g <- rbind(`Poisson GLM` = c(AIC(fit_g0), pdisp(fit_g0), coef(fit_g0)["HEIGHT_s"], sqrt(vcov(fit_g0)["HEIGHT_s", "HEIGHT_s"])),
               `+ brood, location` = c(AIC(fit_g1), pdisp(fit_g1), fixef(fit_g1)["HEIGHT_s"], sqrt(vcov(fit_g1)["HEIGHT_s", "HEIGHT_s"])),
               `+ OLRE` = c(AIC(fit_g2), pdisp(fit_g2), fixef(fit_g2)["HEIGHT_s"], sqrt(vcov(fit_g2)["HEIGHT_s", "HEIGHT_s"])))
colnames(tab_g) <- c("AIC", "Pearson/df", "beta_height", "se")
round(tab_g, 3)
VarCorr(fit_g2)

set.seed(3)
p_hat <- fitted(fit_lap)                                       # conditional probabilities
Fhi <- pbinom(cbpp$incidence, cbpp$size, p_hat); Flo <- pbinom(cbpp$incidence - 1, cbpp$size, p_hat)
u <- Flo + runif(nrow(cbpp)) * (Fhi - Flo); r_q <- qnorm(u)
ggplot(data.frame(r_q), aes(sample = r_q)) + stat_qq(colour = ie_navy) + stat_qq_line(colour = ie_coral) +
  labs(x = "theoretical quantiles", y = "randomised quantile residuals")
c(mean = mean(r_q), sd = sd(r_q))

# fit  <- glmer(cbind(cases, size - cases) ~ x + (1 | g), family = binomial, data = d)   # Laplace (nAGQ = 1)
# fit25 <- update(fit, nAGQ = 25)                       # adaptive Gauss-Hermite, scalar random effect only
# fitP <- MASS::glmmPQL(y ~ x, random = ~ 1 | g, family = poisson, data = d)             # PQL
# summary(fit); fixef(fit); VarCorr(fit); ranef(fit, condVar = TRUE)                    # estimates, predictions
# anova(fit_small, fit_big)                              # LRT with the same nAGQ; boundary for variances
# confint(fit, method = "profile")                       # also "Wald", "boot"
# predict(fit, newdata = nd, re.form = NA, type = "response")    # typical group (b = 0), NOT the population average
# simulate(fit, nsim = 1000)                             # parametric bootstrap, randomised quantile residuals
# glmer(y ~ x + (1 | g) + (1 | obs), family = poisson, data = d)   # OLRE: obs = 1:nrow(d)
# glmer.nb(y ~ x + (1 | g), data = d)                    # negative binomial GLMM

m <- 6 * exp(.5); c(mean = m, var = m + m^2 * (exp(1) - 1), cov = m^2 * (exp(1) - 1), cor = m^2 * (exp(1) - 1) / (m + m^2 * (exp(1) - 1)))

st <- function(n) sqrt(2 * pi * n) * (n / exp(1))^n
rbind(n = c(5, 10), Stirling = st(c(5, 10)), factorial = factorial(c(5, 10)), rel_error = st(c(5, 10)) / factorial(c(5, 10)) - 1)

tau2 <- .5; fh <- function(b) dpois(1, exp(b)) * dnorm(b, 0, sqrt(tau2))
ex <- integrate(fh, -Inf, Inf, rel.tol = 1e-12)$value
bh <- 0; Hh <- exp(bh) + 1 / tau2
sapply(c(exact = 0, K1 = 1, K3 = 3, K5 = 5), function(K) if (K == 0) ex else quad(fh, bh, Hh, K))

sig <- c(.5, 1, 2, 3); xs <- seq(-3, 3, by = .05)
slope <- sapply(sig, function(s) {
  p <- sapply(xs, function(x0) integrate(function(b) plogis(x0 + b) * dnorm(b, 0, s), -Inf, Inf)$value)
  coef(glm(p ~ xs, family = quasibinomial))[[2]]
})
round(cbind(sigma_b = sig, best_logistic_slope = slope, approx = 1 / sqrt(1 + .346 * sig^2), naive_probit = 1 / sqrt(1 + sig^2)), 3)

data(epil, package = "MASS")
epil$lbase_c <- log(epil$base / 4) - mean(log(epil$base / 4)); epil$lage_c <- log(epil$age) - mean(log(epil$age))
epil$trt <- relevel(factor(epil$trt), "placebo"); epil$obs <- 1:nrow(epil)
ep0 <- glm(y ~ lbase_c * trt + lage_c + V4, family = poisson, data = epil)
ep1 <- glmer(y ~ lbase_c * trt + lage_c + V4 + (1 | subject), family = poisson, data = epil)
ep2 <- update(ep1, . ~ . + (1 | obs))
c(pearson_glm = sum(residuals(ep0, "pearson")^2) / ep0$df.residual)
AIC(ep0, ep1, ep2)
T_ol <- as.numeric(2 * (logLik(ep2) - logLik(ep1))); c(T = T_ol, p_mixture = .5 * pchisq(T_ol, 1, lower.tail = FALSE))
round(rbind(glm = coef(ep0), subject = fixef(ep1), subject_plus_OLRE = fixef(ep2))[, c("trtprogabide", "lbase_c", "lbase_c:trtprogabide")], 3)
round(rbind(glm = sqrt(diag(vcov(ep0))), subject = sqrt(diag(vcov(ep1))), olre = sqrt(diag(vcov(ep2))))[, "trtprogabide"], 3)
VarCorr(ep2)

set.seed(2)
sdl <- function(f) if (inherits(f, "lme")) as.numeric(VarCorr(f)[1, 2]) else attr(VarCorr(f)$id, "stddev")[[1]]
res <- replicate(40, {
  m <- 150; n <- 3; id <- rep(1:m, each = n); x <- rnorm(m * n); b <- rnorm(m, 0, 2)[id]
  d <- data.frame(y = rbinom(m * n, 1, plogis(-1 + x + b)), x = x, id = factor(id))
  f1 <- glmer(y ~ x + (1 | id), family = binomial, data = d); f2 <- update(f1, nAGQ = 15)
  f3 <- MASS::glmmPQL(y ~ x, ~ 1 | id, binomial, d, verbose = FALSE)
  c(slope_Laplace = fixef(f1)[[2]], slope_AGQ = fixef(f2)[[2]], slope_PQL = fixef(f3)[[2]],
    sd_Laplace = sdl(f1), sd_AGQ = sdl(f2), sd_PQL = sdl(f3))
})
round(rbind(mean = rowMeans(res), sd = apply(res, 1, sd)), 3)

nb <- MASS::glm.nb(TICKS ~ YEAR + HEIGHT_s, data = gt)
g_b <- glmer(TICKS ~ YEAR + HEIGHT_s + (1 | BROOD), family = poisson, data = gt)
round(rbind(`NB GLM` = c(AIC(nb), coef(nb)["HEIGHT_s"], sqrt(vcov(nb)["HEIGHT_s", "HEIGHT_s"])),
            `Poisson + brood` = c(AIC(g_b), fixef(g_b)["HEIGHT_s"], sqrt(vcov(g_b)["HEIGHT_s", "HEIGHT_s"])),
            `Poisson + brood + location + OLRE` = c(AIC(fit_g2), fixef(fit_g2)["HEIGHT_s"], sqrt(vcov(fit_g2)["HEIGHT_s", "HEIGHT_s"]))), 3)
