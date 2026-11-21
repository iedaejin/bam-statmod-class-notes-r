source("_common.R")

library(mgcv)
library(MASS)
data(airquality)
aq <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])

long <- rbind(data.frame(var = "Solar.R (Langleys)", x = aq$Solar.R, y = aq$Ozone),
              data.frame(var = "Wind (mph)", x = aq$Wind, y = aq$Ozone),
              data.frame(var = "Temp (F)", x = aq$Temp, y = aq$Ozone))
ggplot(long, aes(x, y)) + geom_point(colour = ie_navy, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ie_coral) +
  facet_wrap(~var, scales = "free_x") + labs(x = NULL, y = "Ozone (ppb)")

lm_raw <- lm(Ozone ~ Solar.R + Wind + Temp, data = aq)
round(coef(summary(lm_raw)), 4)
c(R2 = summary(lm_raw)$r.squared, sigma = summary(lm_raw)$sigma)

ggplot(data.frame(fit = fitted(lm_raw), res = resid(lm_raw)), aes(fit, res)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_point(colour = ie_navy, alpha = .7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ie_coral) +
  labs(x = "Fitted ozone (ppb)", y = "Residual (ppb)")

lm_log <- lm(log(Ozone) ~ Solar.R + Wind + Temp, data = aq)
rl <- rbind(data.frame(var = "Solar.R", x = aq$Solar.R, r = resid(lm_log)),
            data.frame(var = "Wind", x = aq$Wind, r = resid(lm_log)),
            data.frame(var = "Temp", x = aq$Temp, r = resid(lm_log)))
ggplot(rl, aes(x, r)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_point(colour = ie_navy, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, colour = ie_coral, fill = "grey80") +
  facet_wrap(~var, scales = "free_x") + labs(x = NULL, y = "Residual, log scale")

y  <- log(aq$Ozone); xs <- aq[, c("Solar.R", "Wind", "Temp")]
n  <- length(y); beta0 <- mean(y)
f  <- matrix(0, n, 3, dimnames = list(NULL, names(xs)))
for (it in 1:100) {
  f_old <- f
  for (j in 1:3) {
    r  <- y - beta0 - rowSums(f[, -j, drop = FALSE])      # partial residual
    fj <- predict(smooth.spline(xs[[j]], r, df = 3), xs[[j]])$y
    f[, j] <- fj - mean(fj)                               # centre
  }
  if (max(abs(f - f_old)) < 1e-8) break
}
c(sweeps = it, beta0 = beta0)

g  <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "REML")
X  <- model.matrix(g)
A  <- lapply(1:3, function(j) {
  sm <- g$smooth[[j]]; ii <- sm$first.para:sm$last.para; Xj <- X[, ii]
  Xj %*% solve(crossprod(Xj) + g$sp[j] * sm$S[[1]], t(Xj))     # A_j = B_j (B_j'B_j + lambda_j S_j)^{-1} B_j'
})
f2 <- matrix(0, n, 3); tm <- predict(g, type = "terms")
err <- numeric(25)
for (it in 1:25) {
  for (j in 1:3) f2[, j] <- A[[j]] %*% (y - beta0 - rowSums(f2[, -j, drop = FALSE]))
  err[it] <- max(abs(f2 - tm))
}
signif(err[c(1, 5, 10, 15, 20, 25)], 3)

d <- do.call(rbind, lapply(1:3, function(j) rbind(
  data.frame(var = names(xs)[j], x = xs[[j]], fit = tm[, j], how = "gam, REML"),
  data.frame(var = names(xs)[j], x = xs[[j]], fit = f[, j], how = "backfitting, 3 df"))))
d <- d[order(d$how, d$var, d$x), ]
ggplot(d, aes(x, fit, colour = how, linetype = how)) + geom_line(linewidth = 1) +
  facet_wrap(~var, scales = "free_x") + scale_colour_manual(values = c(ie_navy, ie_coral)) +
  labs(x = NULL, y = "centred effect", colour = NULL, linetype = NULL)

set.seed(16)
n_s <- 200
sweeps <- function(r) {
  x1 <- runif(n_s); x2 <- r * x1 + sqrt(1 - r^2) * runif(n_s)
  yy <- sin(2 * pi * x1) + x2 + rnorm(n_s, sd = .3)
  gg <- gam(yy ~ s(x1) + s(x2), sp = c(.01, .01)); Xg <- model.matrix(gg)
  Ag <- lapply(1:2, function(j) {
    sm <- gg$smooth[[j]]; ii <- sm$first.para:sm$last.para; Xj <- Xg[, ii]
    Xj %*% solve(crossprod(Xj) + gg$full.sp[j] * sm$S[[1]], t(Xj))
  })
  ff <- matrix(0, n_s, 2)
  for (it in 1:5000) {
    old <- ff
    for (j in 1:2) ff[, j] <- Ag[[j]] %*% (yy - mean(yy) - rowSums(ff[, -j, drop = FALSE]))
    if (max(abs(ff - old)) < 1e-8) break
  }
  c(cor = cor(x1, x2), sweeps = it,
    concurvity = concurvity(gg, full = FALSE)$worst[2, 3])
}
res <- t(sapply(c(0, .5, .9, .99), sweeps))
data.frame(r = c(0, .5, .9, .99), round(res, 3))

m_log <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "REML")
summary(m_log)

par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))
plot(m_log, shade = TRUE, residuals = TRUE, pch = 16, cex = .5, seWithMean = TRUE,
     shade.col = "grey85")

nd <- data.frame(Solar.R = median(aq$Solar.R), Wind = median(aq$Wind), Temp = c(70, 90))
pr <- predict(m_log, nd, se.fit = TRUE)
cbind(nd, eta = pr$fit, se = pr$se.fit, ppb = exp(pr$fit))
Xp <- predict(m_log, nd, type = "lpmatrix")
d  <- Xp[2, ] - Xp[1, ]
est <- sum(d * coef(m_log)); se <- sqrt(drop(t(d) %*% vcov(m_log) %*% d))
c(difference = est, se = se, ratio = exp(est),
  lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se))

m_gam <- gam(Ozone ~ s(Solar.R) + s(Wind) + s(Temp), family = Gamma(link = "log"),
             data = aq, method = "REML")
round(rbind(edf = summary(m_gam)$edf, `edf log-response model` = summary(m_log)$edf[1:3]), 2)
round(m_gam$scale, 3)

m_id <- gam(Ozone ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "REML")
c(gaussian_identity   = AIC(m_id),
  gaussian_on_log     = AIC(m_log) + 2 * sum(log(aq$Ozone)),
  gamma_log_link      = AIC(m_gam))

pima <- rbind(Pima.tr, Pima.te)
m_glm <- glm(type ~ glu + bmi + age + ped, family = binomial, data = pima)
m_bin <- gam(type ~ s(glu) + s(bmi) + s(age) + s(ped), family = binomial,
             data = pima, method = "REML")
round(summary(m_bin)$s.table, 4)
AIC(m_glm, m_bin)

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(m_bin, select = 1, shade = TRUE, ylab = "effect on logit", shade.col = "grey85")
plot(m_bin, select = 3, shade = TRUE, ylab = "effect on logit", shade.col = "grey85")

m_lin <- gam(log(Ozone) ~ Solar.R + Wind + Temp, data = aq, method = "REML")
anova(m_lin, m_log, test = "F")
AIC(m_lin, m_log)

m_bin_lin_age <- gam(type ~ s(glu) + s(bmi) + age + s(ped), family = binomial,
                     data = pima, method = "REML")
anova(m_bin_lin_age, m_bin, test = "Chisq")

set.seed(3); n_r <- 100
x1 <- rnorm(n_r); x2 <- .8 * x1 + .6 * rnorm(n_r); x1 <- x1 - mean(x1); x2 <- x2 - mean(x2)
yy <- x1 + x2 + rnorm(n_r)
cf <- coef(lm(yy ~ x1 + x2 - 1))
a <- 0; b <- 0; ea <- numeric(8)
for (t in 1:8) {
  a <- sum(x1 * (yy - b * x2)) / sum(x1^2)
  b <- sum(x2 * (yy - a * x1)) / sum(x2^2)
  ea[t] <- a - cf[1]
}
c(ratio = unname((ea[8] / ea[7])), r2 = cor(x1, x2)^2)

nd <- data.frame(Solar.R = median(aq$Solar.R), Wind = c(5, 15), Temp = median(aq$Temp))
Xp <- predict(m_log, nd, type = "lpmatrix"); d <- Xp[2, ] - Xp[1, ]
est <- sum(d * coef(m_log)); se <- sqrt(drop(t(d) %*% vcov(m_log) %*% d))
round(c(est = est, se = se, ratio = exp(est), lo = exp(est - 1.96 * se),
        hi = exp(est + 1.96 * se)), 3)
round(unname(exp(10 * coef(lm_log)["Wind"])), 3)

c(log_scale = AIC(m_log), jacobian = 2 * sum(log(aq$Ozone)),
  on_ozone_scale = AIC(m_log) + 2 * sum(log(aq$Ozone)), gamma = AIC(m_gam))

dd <- data.frame(year = 1860:1959, count = as.numeric(discoveries))
g_lin <- glm(count ~ year, family = poisson, data = dd)
g_sm  <- gam(count ~ s(year), family = poisson, data = dd, method = "REML")
AIC(g_lin, g_sm)
round(summary(g_sm)$s.table, 3)
par(mar = c(4, 4, 1, 1)); plot(g_sm, shade = TRUE, shade.col = "grey85", ylab = "effect on log-rate")
disp <- sum(residuals(g_sm, type = "pearson")^2) / df.residual(g_sm)
round(disp, 3)

nd <- data.frame(glu = c(100, 150), bmi = 32, age = 30, ped = 0.4)
Xp <- predict(m_bin, nd, type = "lpmatrix"); d <- Xp[2, ] - Xp[1, ]
est <- sum(d * coef(m_bin)); se <- sqrt(drop(t(d) %*% vcov(m_bin) %*% d))
round(exp(c(or = est, lo = est - 1.96 * se, hi = est + 1.96 * se)), 2)
round(unname(exp(50 * coef(m_glm)["glu"])), 2)
round(predict(m_bin, data.frame(glu = 120, bmi = 32, age = c(25, 45), ped = .4),
              type = "response"), 3)

m_wl <- gam(log(Ozone) ~ s(Solar.R) + Wind + s(Temp), data = aq, method = "REML")
anova(m_wl, m_log, test = "F")
c(AIC_linear_wind = AIC(m_wl), AIC_smooth_wind = AIC(m_log))

y2 <- swiss$Fertility; X2 <- swiss[, c("Education", "Catholic")]
gs <- gam(Fertility ~ s(Education) + s(Catholic), data = swiss, method = "REML")
tt <- predict(gs, type = "terms"); e2 <- summary(gs)$edf
backfit <- function(df) {
  b0 <- mean(y2); ff <- matrix(0, nrow(swiss), 2)
  for (it in 1:500) {
    old <- ff
    for (j in 1:2) {
      r <- y2 - b0 - rowSums(ff[, -j, drop = FALSE])
      s <- predict(smooth.spline(X2[[j]], r, df = df[j]), X2[[j]])$y
      ff[, j] <- s - mean(s)
    }
    if (max(abs(ff - old)) < 1e-8) break
  }
  c(sweeps = it, max_diff = max(abs(ff - tt)))
}
rbind(df_4 = backfit(c(4, 4)), df_edf_plus_1 = backfit(e2 + 1))
round(rbind(edf = e2, sd_effect = apply(tt, 2, sd)), 2)
