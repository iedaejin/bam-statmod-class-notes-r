source("_common.R")
suppressPackageStartupMessages({library(MASS); library(faraway)})

sh <- subset(ships, service > 0)
sh$year <- factor(sh$year); sh$period <- factor(sh$period)
fs <- glm(incidents ~ type + year + period + offset(log(service)), family = poisson, data = sh)
summary(fs)

bw <- with(birthwt, data.frame(low = low, age = age, lwt = lwt,
  race = factor(race, labels = c("white", "black", "other")),
  smoke = smoke, ptd = as.integer(ptl > 0), ht = ht, ui = ui))
form <- low ~ age + lwt + race + smoke + ptd + ht + ui
fpr <- glm(form, family = binomial(link = "probit"), data = bw)
X <- model.matrix(form, bw); y <- bw$low; mu <- fitted(fpr); eta <- drop(X %*% coef(fpr))
loglik <- function(b) sum(dbinom(y, 1, pnorm(drop(X %*% b)), log = TRUE))
H <- -optimHess(coef(fpr), loglik)                       # numerical observed information
g1 <- 1 / dnorm(eta); g2 <- eta / dnorm(eta)^2           # g'(mu), g''(mu) for g = qnorm
V <- mu * (1 - mu); V1 <- 1 - 2 * mu; w <- 1 / (V * g1^2)
I_exp <- crossprod(X, w * X)
I_obs <- crossprod(X, (w + (y - mu) * (V1 * g1 + V * g2) / (V^2 * g1^3)) * X)
rel <- function(A, B) max(abs(A - B)) / max(abs(B))
c(formula_vs_numerical = rel(I_obs, H), expected_vs_numerical = rel(I_exp, H))

irls <- function(X, y, family = binomial(), offset = rep(0, length(y)), weights = rep(1, length(y)),
                 tol = 1e-10, maxit = 50) {
  mu <- (y + mean(y)) / 2;  eta <- family$linkfun(mu)       # starting values inside the range
  dev <- sum(family$dev.resids(y, mu, weights));  trace <- numeric(0);  converged <- FALSE
  for (iter in 1:maxit) {
    d    <- family$mu.eta(eta)                              # d mu / d eta = 1 / g'(mu)
    z    <- eta - offset + (y - mu) / d                     # working response
    w    <- weights * d^2 / family$variance(mu)             # working weights
    beta <- solve(crossprod(X, w * X), crossprod(X, w * z))
    eta  <- drop(X %*% beta) + offset                       # update with the NEW beta
    mu   <- family$linkinv(eta)
    dev_new <- sum(family$dev.resids(y, mu, weights));  trace <- c(trace, dev_new)
    if (abs(dev_new - dev) / (abs(dev_new) + 0.1) < tol) { converged <- TRUE; dev <- dev_new; break }
    dev <- dev_new
  }
  w <- weights * family$mu.eta(eta)^2 / family$variance(mu)  # weights at the final beta
  list(coefficients = drop(beta), deviance = dev, iter = iter, converged = converged,
       trace = trace, cov_unscaled = solve(crossprod(X, w * X)))
}
my  <- irls(X, y, binomial())
fit <- glm(form, family = binomial, data = bw, control = glm.control(epsilon = 1e-10))
c(iterations = my$iter, converged = my$converged, max_abs_diff = max(abs(my$coefficients - coef(fit))))
round(my$trace, 4)
max(abs(sqrt(diag(my$cov_unscaled)) - sqrt(diag(vcov(fit)))))

Xs <- model.matrix(~ type + year + period, data = sh)
my_s <- irls(Xs, sh$incidents, poisson(), offset = log(sh$service))
c(iterations = my_s$iter, max_abs_diff = max(abs(my_s$coefficients - coef(fs))))

path <- function(link, maxit = 40, show = 9) {
  fam <- binomial(link = link); mu <- (y + mean(y)) / 2; eta <- fam$linkfun(mu)
  B <- matrix(NA, maxit, ncol(X))
  for (t in 1:maxit) {
    d <- fam$mu.eta(eta); z <- eta + (y - mu) / d; w <- d^2 / fam$variance(mu)
    B[t, ] <- solve(crossprod(X, w * X), crossprod(X, w * z)); eta <- drop(X %*% B[t, ]); mu <- fam$linkinv(eta)
  }
  data.frame(link = link, iter = 1:show, err = pmax(apply(B[1:show, ], 1, function(b) max(abs(b - B[maxit, ]))), 1e-16))
}
ggplot(rbind(path("logit"), path("probit")), aes(iter, err, colour = link)) +
  geom_line(linewidth = 1) + geom_point(size = 2.2) + scale_y_log10() +
  scale_colour_manual(values = ie_pal[1:2]) + labs(x = "IRLS iteration", y = "max |beta(t) - beta-hat|", colour = NULL)

fit <- glm(form, family = binomial, data = bw)
est <- coef(summary(fit))
wald <- exp(cbind(OR = est[, 1], lower = est[, 1] - 1.96 * est[, 2], upper = est[, 1] + 1.96 * est[, 2]))
prof <- exp(confint(fit))
round(cbind(wald, prof_lower = prof[, 1], prof_upper = prof[, 2])[c("smoke", "ht", "ptd"), ], 3)
est["smoke", ]

fit0 <- update(fit, . ~ . - race)
anova(fit0, fit, test = "Chisq")
LR <- deviance(fit0) - deviance(fit); df <- df.residual(fit0) - df.residual(fit)
c(LR = LR, df = df, p = pchisq(LR, df, lower.tail = FALSE), p_by_hand = exp(-LR / 2))
c(AIC_without = AIC(fit0), AIC_with = AIC(fit), McFadden = 1 - deviance(fit) / fit$null.deviance)

anova(fs, test = "Chisq")

X2 <- sum(residuals(fs, type = "pearson")^2)
c(X2 = X2, D = deviance(fs), df = df.residual(fs), phi_pearson = X2 / df.residual(fs),
  phi_deviance = deviance(fs) / df.residual(fs))
fq <- update(fs, family = quasipoisson)
round(cbind(est = coef(fs), se_poisson = coef(summary(fs))[, 2], se_quasi = coef(summary(fq))[, 2],
            ratio = coef(summary(fq))[, 2] / coef(summary(fs))[, 2]), 3)[1:4, ]
anova(update(fq, . ~ . - type), fq, test = "F")

rqr <- function(y, F, ...) {                  # F: fitted CDF, its parameters go in "..."
  lo <- F(y - 1, ...); hi <- F(y, ...); u <- runif(length(y))
  up_lo <- F(y, ..., lower.tail = FALSE); up_hi <- F(y - 1, ..., lower.tail = FALSE)
  ifelse(lo < 0.5, qnorm(lo + u * (hi - lo)), qnorm(up_lo + u * (up_hi - up_lo), lower.tail = FALSE))
}
set.seed(2027)
r_ships <- rqr(sh$incidents, ppois, lambda = fitted(fs))
r_bin   <- rqr(bw$low, pbinom, size = 1, prob = fitted(fit))
c(shapiro_ships = shapiro.test(r_ships)$p.value, shapiro_binary = shapiro.test(r_bin)$p.value)

op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(fitted(fit), residuals(fit, "pearson"), col = ifelse(bw$low == 1, ie_coral, ie_navy), pch = 16,
     xlab = "fitted probability", ylab = "Pearson residual")
legend("topright", c("y = 1", "y = 0"), col = c(ie_coral, ie_navy), pch = 16, bty = "n")
qqnorm(r_bin, main = "", pch = 16, col = ie_navy, xlab = "N(0,1) quantile", ylab = "randomized quantile residual"); qqline(r_bin, col = ie_coral)
par(op)

rs <- data.frame(y = sh$incidents, mu = fitted(fs), pearson = residuals(fs, "pearson"),
                 deviance = residuals(fs, "deviance"), rstd = rstandard(fs, type = "pearson"),
                 h = hatvalues(fs), cook = cooks.distance(fs), rq = r_ships)
round(head(rs[order(-abs(rs$pearson)), ], 4), 3)
c(sum_pearson2 = sum(rs$pearson^2), sum_deviance2 = sum(rs$deviance^2), trace_H = sum(rs$h))

fit_probit <- update(fit, family = binomial(link = "probit"))
fit_cll <- update(fit, family = binomial(link = "cloglog"))
sapply(list(logit = fit, probit = fit_probit, cloglog = fit_cll), AIC)

xx <- 0:5; yy <- c(0, 0, 1, 0, 1, 1); Xh <- cbind(1, xx)
zz <- 4 * yy - 2                                                  # working response at beta = 0, w = 1/4
solve(crossprod(Xh, 0.25 * Xh), crossprod(Xh, 0.25 * zz))        # one step from beta = 0
full <- irls(Xh, yy, binomial()); rbind(irls = full$coefficients, glm = coef(glm(yy ~ xx, binomial)))
full$iter

yy <- 59; n <- 189
c(by_hand = -2 * (yy * log(yy / n) + (n - yy) * log(1 - yy / n)), glm = fit$null.deviance)
1 - deviance(fit) / fit$null.deviance

c(F_seq_type = (55.44 / 4) / (42.28 / 25), p = pf((55.44 / 4) / (42.28 / 25), 4, 25, lower.tail = FALSE))
c(F_marginal_type = (23.67 / 4) / (42.28 / 25), p = pf((23.67 / 4) / (42.28 / 25), 4, 25, lower.tail = FALSE))

f <- function(y, mu) c(pearson = (y - mu) / sqrt(mu),
  deviance = sign(y - mu) * sqrt(2 * ifelse(y == 0, mu, y * log(y / mu) - (y - mu))))
round(rbind(`y=0, mu=1.3` = f(0, 1.3), `y=6, mu=1.5` = f(6, 1.5)), 3)
c(qnorm(exp(-1.3)), qnorm(ppois(5, 1.5)), qnorm(ppois(6, 1.5)))

Xt <- model.matrix(~ log(Girth) + log(Height), trees)
mg <- irls(Xt, trees$Volume, Gamma(link = "log"))
gg <- glm(Volume ~ log(Girth) + log(Height), family = Gamma(link = "log"), data = trees)
max(abs(mg$coefficients - coef(gg)))
mu_g <- exp(drop(Xt %*% mg$coefficients))
phi <- sum((trees$Volume - mu_g)^2 / mu_g^2) / (nrow(trees) - 3)
c(phi_irls = phi, phi_summary = summary(gg)$dispersion)
cbind(se_irls = sqrt(phi * diag(mg$cov_unscaled)), se_summary = coef(summary(gg))[, 2])

o <- glm(cbind(damage, 6 - damage) ~ temp, family = binomial, data = orings)
b <- coef(summary(o))["temp", ]
wald <- b[1] + c(-1, 1) * 1.96 * b[2]
prof <- confint(o)["temp", ]
round(rbind(wald = wald, profile = prof, `OR per 10F, Wald` = exp(10 * wald), `OR per 10F, profile` = exp(10 * prof)), 3)
