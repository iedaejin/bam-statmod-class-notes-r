source("_common.R")

library(mgcv)
library(nlme)
library(splines)
ox <- as.data.frame(nlme::Oxboys)       # 26 boys, 9 occasions
ox$Subject <- factor(as.character(ox$Subject))

ggplot(ox, aes(age, height, group = Subject)) +
  geom_line(colour = ie_navy, alpha = .55) + geom_point(colour = ie_navy, size = .9, alpha = .7) +
  labs(x = "Standardised age", y = "Height (cm)")

m_naive <- gam(height ~ s(age), data = ox, method = "REML")
lag1 <- function(r, id) {                # within-subject lag-1 correlation of residuals
  d <- split(r, id)
  cor(unlist(lapply(d, function(v) v[-length(v)])), unlist(lapply(d, function(v) v[-1])))
}
c(edf = sum(m_naive$edf) - 1, sigma = sqrt(m_naive$sig2), lag1_within_boy = lag1(residuals(m_naive), ox$Subject))

mc <- MASS::mcycle; y <- mc$accel; x <- mc$times; n <- length(y)
nseg <- 22; deg <- 3; dx <- diff(range(x)) / nseg
knots <- seq(min(x) - deg * dx, max(x) + deg * dx, by = dx)
B <- splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
K <- ncol(B)
S <- crossprod(diff(diag(K), differences = 2))          # S = D'D, rank K - 2
e <- eigen(S, symmetric = TRUE); pos <- e$values > 1e-8
c(K = K, positive = sum(pos), zero = sum(!pos))

fit_pen <- gam(y ~ B - 1, paraPen = list(B = list(S)), method = "REML")
Z <- B %*% e$vectors[, pos] %*% diag(1 / sqrt(e$values[pos]))
g <- factor(rep(1, n))
fit_lme <- lme(y ~ x, random = list(g = pdIdent(~ Z - 1)), method = "REML")
s2 <- fit_lme$sigma^2; s2b <- as.numeric(VarCorr(fit_lme)[1, "Variance"])
round(c(sigma2 = s2, sigma2_b = s2b, ratio = s2 / s2b, lambda_gam = unname(fit_pen$sp),
        edf = sum(fit_pen$edf)), 4)
max(abs(fitted(fit_lme) - fitted(fit_pen)))

dd <- data.frame(x = x, y = y, mixed = fitted(fit_lme), fixed = fitted(fit_lme, level = 0))
ggplot(dd, aes(x, y)) + geom_point(colour = "grey55") +
  geom_line(aes(y = fixed), colour = ie_coral, linetype = 2, linewidth = 1) +
  geom_line(aes(y = mixed), colour = ie_navy, linewidth = 1.1) +
  labs(x = "Time after impact (ms)", y = "Head acceleration (g)")

m1 <- gam(height ~ s(age) + s(Subject, bs = "re"), data = ox, method = "REML")
summary(m1)$s.table
gam.vcomp(m1)

co <- coef(m1); b_blup <- co[grep("Subject", names(co))]
raw <- tapply(ox$height - predict(m1, exclude = "s(Subject)"), ox$Subject, mean)
c(sp = unname(m1$sp["s(Subject)"]), sigma2_over_sigma2b = m1$sig2 / gam.vcomp(m1)$vc["s(Subject)", "std.dev"]^2,
  shrinkage_9_over_9_plus_lambda = 9 / (9 + unname(m1$sp["s(Subject)"])))
range(b_blup / raw)

m2 <- gam(height ~ s(age) + s(Subject, bs = "re") + s(Subject, age, bs = "re"), data = ox, method = "REML")
m3 <- suppressWarnings(gam(height ~ s(age, k = 5) + s(age, Subject, bs = "fs", k = 5), data = ox, method = "REML"))
gam.vcomp(m2)
cmp <- data.frame(model = c("random intercepts", "+ random slopes", "factor smooths (fs)"),
                  conditional_AIC = c(AIC(m1), AIC(m2), AIC(m3)),
                  lag1_within_boy = c(lag1(residuals(m1), ox$Subject), lag1(residuals(m2), ox$Subject),
                                      lag1(residuals(m3), ox$Subject)))
knitr::kable(cmp, digits = c(0, 1, 2), booktabs = TRUE,
             caption = "Three random-effect structures for the Oxboys growth curves.")

nd <- expand.grid(age = seq(-1, 1, length.out = 40), Subject = levels(ox$Subject))
nd1 <- nd[nd$Subject == levels(ox$Subject)[1], ]
d4 <- rbind(data.frame(nd, fit = predict(m1, nd), model = "random intercepts"),
            data.frame(nd, fit = predict(m3, nd), model = "factor smooths"))
pop <- rbind(data.frame(age = nd1$age, fit = predict(m1, nd1, exclude = "s(Subject)"), model = "random intercepts"),
             data.frame(age = nd1$age, fit = predict(m3, nd1, exclude = "s(age,Subject)"), model = "factor smooths"))
ggplot(d4, aes(age, fit, group = Subject)) + geom_line(colour = ie_blue, alpha = .6, linewidth = .5) +
  geom_line(data = pop, aes(group = NULL), colour = ie_coral, linewidth = 1.3) +
  facet_wrap(~model) + labs(x = "Standardised age", y = "Height (cm)")

g2 <- gamm(height ~ s(age), random = list(Subject = ~ 1 + age), data = ox, method = "REML")
names(g2)
VarCorr(g2$lme)
summary(g2$gam)$s.table

library(gamm4)
g4 <- gamm4(height ~ s(age), random = ~(1 + age | Subject), data = ox, REML = TRUE)
lme4::VarCorr(g4$mer)

data(sleepstudy, package = "lme4"); ss <- sleepstudy; ss$Subject <- factor(as.character(ss$Subject))
lm4 <- lme4::lmer(Reaction ~ Days + (Days | Subject), ss)
ga  <- gam(Reaction ~ s(Days, k = 5) + s(Subject, bs = "re") + s(Subject, Days, bs = "re"), data = ss, method = "REML")
gm  <- gamm(Reaction ~ s(Days, k = 5), random = list(Subject = ~ 1 + Days), data = ss, method = "REML")
g4s <- gamm4(Reaction ~ s(Days, k = 5), random = ~(1 + Days | Subject), data = ss, REML = TRUE)
vcs <- function(sd0, sd1, sdE) c(intercept_sd = sd0, slope_sd = sd1, residual_sd = sdE)
vc4 <- as.data.frame(lme4::VarCorr(lm4)); vcg <- lme4::VarCorr(g4s$mer)
rbind(lmer_linear = vcs(vc4$sdcor[1], vc4$sdcor[2], vc4$sdcor[4]),
      gam_re      = vcs(gam.vcomp(ga)$vc["s(Subject)", "std.dev"], gam.vcomp(ga)$vc["s(Subject,Days)", "std.dev"], sqrt(ga$sig2)),
      gamm        = vcs(as.numeric(VarCorr(gm$lme)["(Intercept)", "StdDev"]), as.numeric(VarCorr(gm$lme)["Days", "StdDev"]), gm$lme$sigma),
      gamm4       = vcs(attr(vcg$Subject, "stddev")[1], attr(vcg$Subject, "stddev")[2], attr(vcg, "sc")))
summary(gm$gam)$s.table

days <- data.frame(day = c(1, 2, 3, 4, 7, 8))
round(corMatrix(Initialize(corAR1(0.5, form = ~ day), data = days)), 3)

aq <- airquality; aq$day <- seq_len(nrow(aq)); aq <- na.omit(aq)
f <- log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp)
m0  <- gamm(f, data = aq, method = "REML")
mar <- gamm(f, data = aq, method = "REML", correlation = corAR1(form = ~ day))
anova(m0$lme, mar$lme)
phi_hat <- coef(mar$lme$modelStruct$corStruct, unconstrained = FALSE)
round(rbind(independent = summary(m0$gam)$edf, AR1 = summary(mar$gam)$edf), 2)

set.seed(19); n <- 200; tt <- 1:n; ftrue <- function(t) 2 * sin(2 * pi * t / 100); phi <- 0.7
sim_one <- function() {
  e <- as.numeric(arima.sim(list(ar = phi), n, sd = sqrt(1 - phi^2)))
  d <- data.frame(y = ftrue(tt) + e, tt = tt)
  fa <- gam(y ~ s(tt, k = 30), data = d, method = "REML")
  fb <- gamm(y ~ s(tt, k = 30), data = d, correlation = corAR1(form = ~ tt), method = "REML")
  fc <- gamm(y ~ s(tt, k = 30), data = d, correlation = corAR1(value = phi, fixed = TRUE, form = ~ tt), method = "REML")
  c(edf_ignore = sum(fa$edf) - 1, edf_AR1 = sum(fb$gam$edf) - 1, edf_AR1_fixed = sum(fc$gam$edf) - 1,
    phi_hat = unname(coef(fb$lme$modelStruct$corStruct, unconstrained = FALSE)),
    mse_ignore = mean((fitted(fa) - ftrue(tt))^2), mse_AR1 = mean((fitted(fb$gam) - ftrue(tt))^2),
    mse_AR1_fixed = mean((fitted(fc$gam) - ftrue(tt))^2))
}
sim <- t(replicate(15, sim_one()))
round(colMeans(sim), 2)

ox$occ <- as.integer(ox$Occasion)
gc1 <- gamm(height ~ s(age), random = list(Subject = ~ 1),
            correlation = corAR1(form = ~ occ | Subject), data = ox, method = "REML")
c(phi = unname(coef(gc1$lme$modelStruct$corStruct, unconstrained = FALSE)),
  sd_intercept = as.numeric(VarCorr(gc1$lme)["(Intercept)", "StdDev"]))

ggplot(data.frame(b = unname(b_blup)), aes(sample = b)) + stat_qq(colour = ie_navy) + stat_qq_line(colour = ie_coral) +
  labs(x = "Theoretical quantiles", y = "Predicted random intercept (cm)")

set.seed(5)
b <- rnorm(2e6, sd = 1)
c(simulated = mean(exp(0.5 + b)) / exp(0.5), formula = exp(1 / 2))          # log link, sigma_b = 1
sl <- sapply(c(1, 2), function(sb) {                                       # logit link: slope of the marginal logit
  m <- function(x) integrate(function(u) plogis(x + u) * dnorm(u, sd = sb), -Inf, Inf)$value
  xs <- seq(-2, 2, length.out = 41); unname(coef(lm(qlogis(sapply(xs, m)) ~ xs))[2]) })
rbind(numerical = sl, approximation = 1 / sqrt(1 + 0.346 * c(1, 2)^2))

aqs <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])
y1 <- log(aqs$Ozone); x1 <- aqs$Temp; n1 <- length(y1)
nseg <- 17; dx1 <- diff(range(x1)) / nseg
kn1 <- seq(min(x1) - 3 * dx1, max(x1) + 3 * dx1, by = dx1)
B1 <- splineDesign(kn1, x1, ord = 4, outer.ok = TRUE); K1 <- ncol(B1)
S1 <- crossprod(diff(diag(K1), differences = 2))
e1 <- eigen(S1, symmetric = TRUE); p1 <- e1$values > 1e-8
fp1 <- gam(y1 ~ B1 - 1, paraPen = list(B1 = list(S1)), method = "REML")
Z1 <- B1 %*% e1$vectors[, p1] %*% diag(1 / sqrt(e1$values[p1])); g1 <- factor(rep(1, n1))
fl1 <- lme(y1 ~ x1, random = list(g1 = pdIdent(~ Z1 - 1)), method = "REML")
r1 <- fl1$sigma^2 / as.numeric(VarCorr(fl1)[1, "Variance"])
fg1 <- gamm(y1 ~ s(x1, bs = "ps", k = 20), method = "REML")
vc1 <- VarCorr(fg1$lme)
rg1 <- fg1$lme$sigma^2 / as.numeric(vc1[nrow(vc1) - 1, "Variance"])
round(c(K = K1, ratio_lme = r1, lambda_gam = unname(fp1$sp), edf = sum(fp1$edf),
        gamm_sp = unname(fg1$gam$sp), gamm_ratio = rg1, gamm_edf = sum(fg1$gam$edf)), 3)

g1 <- gamm(height ~ s(age), random = list(Subject = ~ 1), data = ox, method = "REML")
anova(g1$lme, g2$lme)

set.seed(23); n <- 200
res19 <- t(replicate(15, {
  xx <- as.numeric(arima.sim(list(ar = .9), n)); xx <- xx / sd(xx) * 1.2
  e <- as.numeric(arima.sim(list(ar = .5), n, sd = sqrt(1 - .25)))
  d <- data.frame(y = 2 * sin(xx) + e, xx = xx, tt = 1:n)
  a <- gam(y ~ s(xx, k = 20), data = d, method = "REML")
  b <- gamm(y ~ s(xx, k = 20), data = d, correlation = corAR1(form = ~ tt), method = "REML")
  c(edf_ignore = sum(a$edf) - 1, edf_AR1 = sum(b$gam$edf) - 1,
    phi_hat = unname(coef(b$lme$modelStruct$corStruct, unconstrained = FALSE)))
}))
round(colMeans(res19), 2)
