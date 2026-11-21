source("_common.R")

library(survival)
vet <- survival::veteran
vet$trt <- factor(vet$trt, levels = 1:2, labels = c("standard", "test"))
head(vet, 4)
table(status = vet$status, trt = vet$trt)

aml

a1 <- subset(aml, x == "Maintained")
a1 <- a1[order(a1$time), ]
a1$id <- factor(seq_len(nrow(a1)), levels = rev(seq_len(nrow(a1))))
a1$ev <- factor(ifelse(a1$status == 1, "death", "censored"), levels = c("death", "censored"))
ggplot(a1, aes(y = id)) +
  geom_segment(aes(x = 0, xend = time, yend = id), colour = ie_navy) +
  geom_point(aes(x = time, shape = ev, colour = ev), size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(16, 4)) + scale_colour_manual(values = c(ie_coral, ie_navy)) +
  labs(x = "weeks", y = "patient", shape = NULL, colour = NULL)

d <- sum(vet$status); tot <- sum(vet$time)
c(correct = d / tot,
  drop_censored = d / sum(vet$time[vet$status == 1]),
  censored_as_deaths = nrow(vet) / tot)

tt <- seq(0.01, 3, length.out = 300)
shapes <- c(0.7, 1, 1.6)
dd <- do.call(rbind, lapply(shapes, function(k)
  data.frame(t = tt, value = c(k * tt^(k - 1), exp(-tt^k)),
             what = rep(c("hazard h(t)", "survival S(t)"), each = length(tt)),
             shape = paste("k =", k))))
ggplot(dd, aes(t, value, colour = shape)) + geom_line(linewidth = 1) +
  facet_wrap(~what, scales = "free_y") + coord_cartesian(ylim = c(0, NA)) +
  scale_colour_manual(values = ie_pal[1:3]) + labs(x = "t", y = NULL, colour = NULL)

d <- sum(vet$status); tot <- sum(vet$time)
lam <- d / tot; se_lam <- lam / sqrt(d)
c(lambda = lam, se = se_lam, lo = lam * exp(-1.96 / sqrt(d)), hi = lam * exp(1.96 / sqrt(d)),
  mean_days = 1 / lam, median_days = log(2) / lam)
tapply(vet$status, vet$trt, sum) / tapply(vet$time, vet$trt, sum)   # rate by arm

e0 <- survreg(Surv(time, status) ~ 1, data = vet, dist = "exponential")
c(lambda = unname(exp(-coef(e0))), se_log_lambda = unname(sqrt(vcov(e0)[1, 1])),
  one_over_sqrt_d = 1 / sqrt(d))

p1 <- glm(status ~ trt + offset(log(time)), family = poisson, data = vet)
e1 <- survreg(Surv(time, status) ~ trt, data = vet, dist = "exponential")
rbind(poisson_glm = coef(p1), survreg_negated = -coef(e1))

w0 <- survreg(Surv(time, status) ~ 1, data = vet, dist = "weibull")
anova(e0, w0)
c(scale_sigma = w0$scale, shape_k = 1 / w0$scale)

km_all <- survfit(Surv(time, status) ~ 1, data = vet)
tg <- seq(0, 600, by = 2)
sw <- exp(-(tg / exp(coef(w0)))^(1 / w0$scale))
curves <- data.frame(t = tg, exponential = exp(-lam * tg), Weibull = sw)
kmd <- data.frame(t = c(0, km_all$time), s = c(1, km_all$surv))
ggplot() +
  geom_step(data = kmd, aes(t, s, colour = "Kaplan-Meier"), linewidth = 1) +
  geom_line(data = curves, aes(t, exponential, colour = "exponential"), linetype = 2, linewidth = 1) +
  geom_line(data = curves, aes(t, Weibull, colour = "Weibull"), linetype = 2, linewidth = 1) +
  scale_colour_manual(values = c(ie_navy, ie_coral, ie_green)) +
  labs(x = "days", y = "S(t)", colour = NULL)

toy <- data.frame(t = c(3, 4, 6, 8, 8, 10), s = c(1, 1, 0, 1, 1, 1))
summary(survfit(Surv(t, s) ~ 1, data = toy))

am <- subset(aml, x == "Maintained")
tab <- do.call(rbind, lapply(sort(unique(am$time[am$status == 1])), function(tj)
  data.frame(t = tj, n = sum(am$time >= tj), d = sum(am$time == tj & am$status == 1))))
tab$S <- cumprod(1 - tab$d / tab$n)
tab$sum_gw <- cumsum(tab$d / (tab$n * (tab$n - tab$d)))
tab$se <- tab$S * sqrt(tab$sum_gw)
round(tab, 4)
km_aml <- survfit(Surv(time, status) ~ x, data = aml)
summary(km_aml)$table[, c("records", "events", "median", "0.95LCL", "0.95UCL")]

S31 <- tab$S[tab$t == 31]; s31 <- sqrt(tab$sum_gw[tab$t == 31]) / abs(log(S31))
rbind(log = S31 * exp(c(-1, 1) * 1.96 * sqrt(tab$sum_gw[tab$t == 31])),
      loglog = S31^exp(c(1, -1) * 1.96 * s31))
summary(survfit(Surv(time, status) ~ 1, data = am, conf.type = "log-log"))$lower[5]

kmd <- summary(survfit(Surv(time, status) ~ x, data = aml, conf.type = "log-log"), censored = TRUE)
dd <- data.frame(t = kmd$time, s = kmd$surv, lo = kmd$lower, hi = kmd$upper,
                 cens = kmd$n.censor, g = sub("x=", "", as.character(kmd$strata)))
dd0 <- data.frame(t = 0, s = 1, lo = 1, hi = 1, cens = 0, g = unique(dd$g))
dd <- rbind(dd0, dd)
dd <- dd[order(dd$g, dd$t), ]
ggplot(dd, aes(t, s, colour = g, fill = g)) +
  geom_step(linewidth = 1) +
  geom_step(aes(y = lo), linetype = 3, na.rm = TRUE) + geom_step(aes(y = hi), linetype = 3, na.rm = TRUE) +
  geom_point(data = subset(dd, cens > 0), shape = 3, size = 3, stroke = 1.1) +
  scale_colour_manual(values = ie_pal[1:2]) + scale_fill_manual(values = ie_pal[1:2]) +
  labs(x = "weeks", y = "S(t)", colour = NULL, fill = NULL)

tj <- sort(unique(aml$time[aml$status == 1]))
lr <- do.call(rbind, lapply(tj, function(t) {
  at <- aml$time >= t
  n <- sum(at); n1 <- sum(at & aml$x == "Maintained")
  dj <- sum(aml$time == t & aml$status == 1); d1 <- sum(aml$time == t & aml$status == 1 & aml$x == "Maintained")
  data.frame(t = t, n = n, n1 = n1, d = dj, d1 = d1, E1 = n1 * dj / n,
             V1 = ifelse(n > 1, n1 * (n - n1) * dj * (n - dj) / (n^2 * (n - 1)), 0))
}))
round(lr, 3)
O1 <- sum(lr$d1); E1 <- sum(lr$E1); V1 <- sum(lr$V1)
c(O1 = O1, E1 = E1, V1 = V1, X2 = (O1 - E1)^2 / V1, p = pchisq((O1 - E1)^2 / V1, 1, lower.tail = FALSE))
survdiff(Surv(time, status) ~ x, data = aml)

sapply(c(rho0 = 0, rho1 = 1), function(r) {
  sd <- survdiff(Surv(time, status) ~ trt, data = vet, rho = r)
  c(chisq = sd$chisq, p = pchisq(sd$chisq, 1, lower.tail = FALSE))
})

mk <- function(f, d) {
  s <- summary(survfit(f, data = d), censored = TRUE)
  g <- sub(".*=", "", as.character(s$strata))
  x <- data.frame(t = s$time, s = s$surv, g = g)
  x <- rbind(data.frame(t = 0, s = 1, g = unique(g)), x); x[order(x$g, x$t), ]
}
k1 <- mk(Surv(time, status) ~ trt, vet); k2 <- mk(Surv(time, status) ~ celltype, vet)
k1$panel <- "treatment"; k2$panel <- "cell type"
ggplot(rbind(k1, k2), aes(t, s, colour = g)) + geom_step(linewidth = 0.9) +
  facet_wrap(~panel) + scale_colour_manual(values = ie_pal) +
  labs(x = "days", y = "S(t)", colour = NULL)

toy <- data.frame(y = c(2, 3, 4, 6, 8), d = c(1, 0, 1, 1, 1), x = c(1, 0, 0, 1, 0))
U <- function(b) { s <- 0; for (j in which(toy$d == 1)) { R <- toy$y >= toy$y[j]
    w <- exp(b * toy$x[R]); s <- s + toy$x[j] - sum(w * toy$x[R]) / sum(w) }; s }
Info <- function(b) { s <- 0; for (j in which(toy$d == 1)) { R <- toy$y >= toy$y[j]
    w <- exp(b * toy$x[R]); p <- sum(w * toy$x[R]) / sum(w); s <- s + p * (1 - p) }; s }
b <- 0; for (it in 1:4) { b <- b + U(b) / Info(b); print(round(b, 4)) }
c(U0 = U(0), I0 = Info(0), beta_hat = b, se = 1 / sqrt(Info(b)), HR = exp(b),
  score_test = U(0)^2 / Info(0))
round(summary(coxph(Surv(y, d) ~ x, data = toy))$coefficients, 4)

fit_b <- update(coxph(Surv(time, status) ~ trt + karno + celltype + age, data = vet),
                ties = "breslow")
eta <- drop(model.matrix(fit_b) %*% coef(fit_b))
ev <- sort(unique(vet$time[vet$status == 1]))
H0 <- cumsum(sapply(ev, function(t)
  sum(vet$time == t & vet$status == 1) / sum(exp(eta)[vet$time >= t])))
bh <- basehaz(fit_b, centered = FALSE)
c(max_abs_difference = max(abs(bh$hazard[match(ev, bh$time)] - H0)))
round(head(cbind(time = ev, H0 = H0), 4), 4)

fit <- coxph(Surv(time, status) ~ trt + karno + celltype + age, data = vet)   # Efron
fit_e <- update(fit, ties = "exact")
round(cbind(Efron = coef(fit), Breslow = coef(fit_b), Exact = coef(fit_e)), 4)

round(summary(fit)$coefficients[, c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")], 4)
summary(fit)$logtest

fit0 <- coxph(Surv(time, status) ~ trt + karno + age, data = vet)
anova(fit0, fit)

est <- coef(fit); V <- vcov(fit)
hr <- function(contrast, scale = 1) {   # contrast: weights on coefficients
  l <- sum(contrast * est[names(contrast)]) * scale
  Vc <- V[names(contrast), names(contrast)]
  s <- sqrt(drop(t(contrast) %*% Vc %*% contrast)) * abs(scale)
  round(exp(c(HR = l, lower = l - 1.96 * s, upper = l + 1.96 * s)), 3)
}
rbind(`karno, +10 points`     = hr(c(karno = 1), 10),
      `age, +10 years`        = hr(c(age = 1), 10),
      `adeno vs squamous`     = hr(c(celltypeadeno = 1)),
      `adeno vs smallcell`    = hr(c(celltypeadeno = 1, celltypesmallcell = -1)),
      `test vs standard`      = hr(c(trttest = 1)))

nd <- data.frame(trt = "standard", karno = c(40, 80), celltype = "squamous", age = 60)
sf <- summary(survfit(fit, newdata = nd), censored = FALSE)
dp <- rbind(data.frame(t = 0, s = 1, g = "karno 40"), data.frame(t = 0, s = 1, g = "karno 80"),
            data.frame(t = rep(sf$time, 2), s = c(sf$surv), g = rep(c("karno 40", "karno 80"), each = length(sf$time))))
ggplot(dp, aes(t, s, colour = g)) + geom_step(linewidth = 1) + scale_colour_manual(values = ie_pal[1:2]) +
  labs(x = "days", y = "predicted S(t)", colour = NULL)

vet$kgrp <- cut(vet$karno, c(0, 40, 60, 100), labels = c("10-40", "50-60", "70-99"))
s <- summary(survfit(Surv(time, status) ~ kgrp, data = vet), censored = FALSE)
dl <- data.frame(t = s$time, s = s$surv, g = sub("kgrp=", "", as.character(s$strata)))
dl <- subset(dl, s > 0 & s < 1)
ggplot(dl, aes(log(t), log(-log(s)), colour = g)) + geom_step(linewidth = 1) +
  scale_colour_manual(values = ie_pal[1:3]) + labs(x = "log t", y = expression(log(-log(hat(S)))), colour = "Karnofsky")

zp <- cox.zph(fit)
zp

zp_c <- cox.zph(fit, terms = FALSE)               # one column per coef.
r <- residuals(fit, type = "schoenfeld")         # 128 events x 6 coefficients
d_ev <- sum(vet$status)
rstar <- sweep(d_ev * r %*% vcov(fit), 2, coef(fit), "+")
max(abs(rstar - zp_c$y))                          # the same scaled residuals

dz <- data.frame(g = zp_c$x, b = zp_c$y[, "karno"])
ggplot(dz, aes(g, b)) + geom_point(colour = ie_navy, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, colour = ie_coral, fill = ie_coral, alpha = .2) +
  geom_hline(yintercept = coef(fit)["karno"], linetype = 2) + coord_cartesian(ylim = c(-0.2, 0.15)) +
  labs(x = "time (Kaplan-Meier scale)", y = expression(hat(beta)[karno](t)))

vet$kgrp2 <- cut(vet$karno, c(0, 40, 60, 100))
fit_s <- coxph(Surv(time, status) ~ trt + celltype + age + strata(kgrp2), data = vet)
vs <- survSplit(Surv(time, status) ~ ., data = vet, cut = c(90, 180),
                episode = "tgroup", id = "id")
fit_t <- coxph(Surv(tstart, time, status) ~ trt + celltype + age +
                 karno:strata(tgroup), data = vs)
round(coef(summary(fit_t))[6:8, c("coef", "se(coef)", "Pr(>|z|)")], 4)
round(exp(10 * coef(fit_t)[6:8]), 3)       # HR per 10 points, by period
AIC(fit, fit_t)

cox.zph(fit_t)

round(exp(cbind(HR = coef(fit_t), confint(fit_t)))[1:5, ], 2)

set.seed(2027)
u <- runif(1e5); lam <- 0.1; k <- 2
Tsim <- (-log(u))^(1 / k) / lam
c(mean_sim = mean(Tsim), mean_theory = gamma(1 + 1 / k) / lam,
  median_sim = median(Tsim), median_theory = log(2)^(1 / k) / lam)

tt <- c(4, 5, 7, 9, 12, 12, 15, 20); dl <- c(1, 0, 1, 0, 1, 1, 0, 1)
g <- glm(dl ~ 1 + offset(log(tt)), family = poisson)
c(lambda = unname(exp(coef(g))), se_loglambda = unname(sqrt(vcov(g)[1, 1])),
  one_over_sqrt_d = 1 / sqrt(sum(dl)))

set.seed(2027)
n <- 5000; lam <- 0.1; Ti <- rexp(n, lam)
est <- function(C) {
  y <- pmin(Ti, C); d <- Ti <= C
  c(lambda_hat = sum(d) / sum(y), prop_events = mean(d))
}
rbind(administrative = est(rep(10, n)),
      random_independent = est(rexp(n, 0.05)),
      informative = est(ifelse(runif(n) < 0.5, Ti / 2, Inf)))

t9 <- c(2, 3, 3, 5, 5, 6, 8, 9, 12); s9 <- c(1, 1, 0, 1, 1, 0, 1, 0, 1)
km9 <- survfit(Surv(t9, s9) ~ 1, conf.type = "log-log")
summary(km9)
km9

dd <- data.frame(t = c(2, 4, 4, 7, 3, 6, 6, 9, 11), s = c(1, 1, 0, 1, 1, 1, 1, 0, 1),
                 g = rep(c("A", "B"), c(4, 5)))
survdiff(Surv(t, s) ~ g, data = dd)

y <- c(3, 5, 6, 8, 9); dl <- c(1, 1, 0, 1, 1); x <- c(1, 0, 1, 0, 1)
summary(coxph(Surv(y, dl) ~ x))

toy <- data.frame(y = c(2, 3, 4, 6, 8), d = c(1, 0, 1, 1, 1), x = c(1, 0, 0, 1, 0))
ft <- coxph(Surv(y, d) ~ x, data = toy)
basehaz(ft, centered = FALSE)

dt <- data.frame(y = c(1, 1, 2, 3), x = c(1, 0, 1, 0), d = 1)
sapply(c("breslow", "efron", "exact"),
       function(m) coef(coxph(Surv(y, d) ~ x, data = dt, ties = m)))

lg <- na.omit(lung[, c("time", "status", "age", "sex", "ph.ecog")])
lg$ev <- lg$status - 1                               # 1 = death
lg$sex <- factor(lg$sex, 1:2, c("male", "female"))
km <- survfit(Surv(time, ev) ~ sex, data = lg)
km
survdiff(Surv(time, ev) ~ sex, data = lg)

cf <- coxph(Surv(time, ev) ~ sex + age + ph.ecog, data = lg)
round(cbind(HR = exp(coef(cf)), exp(confint(cf))), 3)

cox.zph(cf)

survfit(Surv(time, 1 - ev) ~ 1, data = lg)
