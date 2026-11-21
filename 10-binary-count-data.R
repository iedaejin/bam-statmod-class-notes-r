source("_common.R")
suppressPackageStartupMessages({library(MASS); library(faraway); library(pscl)})

bw <- with(birthwt, data.frame(low = low, age = age, lwt = lwt,
  race = factor(race, labels = c("white", "black", "other")),
  smoke = smoke, ptd = as.integer(ptl > 0), ht = ht, ui = ui))
fit <- glm(low ~ age + lwt + race + smoke + ptd + ht + ui, family = binomial, data = bw)
round(coef(summary(fit)), 4)
round(exp(cbind(OR = coef(fit), confint(fit))), 3)

nd <- data.frame(age = 25, lwt = 120, race = "white", smoke = 0:1, ptd = 0, ht = 0, ui = 0)
eta <- predict(fit, nd); p <- plogis(eta)
rbind(eta = eta, prob = p, odds = p / (1 - p))
c(odds_ratio = unname((p[2] / (1 - p[2])) / (p[1] / (1 - p[1]))), risk_ratio = unname(p[2] / p[1]),
  exp_beta_smoke = unname(exp(coef(fit)["smoke"])))

grid <- expand.grid(age = 25, lwt = 80:250, race = "white", smoke = 0:1, ptd = 0, ht = 0:1, ui = 0)
pr <- predict(fit, grid, se.fit = TRUE)
grid$p <- plogis(pr$fit); grid$lo <- plogis(pr$fit - 1.96 * pr$se.fit); grid$hi <- plogis(pr$fit + 1.96 * pr$se.fit)
grid$smoke <- factor(grid$smoke, labels = c("non-smoker", "smoker"))
grid$ht <- factor(grid$ht, labels = c("no hypertension", "hypertension"))
ggplot(grid, aes(lwt, p, colour = smoke, fill = smoke)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + facet_wrap(~ht) + scale_colour_manual(values = ie_pal[1:2]) + scale_fill_manual(values = ie_pal[1:2]) +
  labs(x = "mother's weight (lb)", y = "P(low birth weight)", colour = NULL, fill = NULL)

fit_probit <- update(fit, family = binomial(link = "probit"))
round(coef(fit) / coef(fit_probit), 2)
c(max_prob_diff = max(abs(fitted(fit) - fitted(fit_probit))), AIC_logit = AIC(fit), AIC_probit = AIC(fit_probit))

p <- fitted(fit); g <- cut(p, quantile(p, 0:5 / 5), include.lowest = TRUE)
cal <- data.frame(n = as.vector(table(g)), mean_pred = tapply(p, g, mean), observed = tapply(bw$low, g, mean))
round(cal, 3)
g10 <- cut(p, quantile(p, 0:10 / 10), include.lowest = TRUE)
O <- tapply(bw$low, g10, sum); n <- tapply(bw$low, g10, length); pb <- tapply(p, g10, mean)
HL <- sum((O - n * pb)^2 / (n * pb * (1 - pb)))
c(HL = HL, p = pchisq(HL, 8, lower.tail = FALSE))

bf <- glm(cbind(disease, nondisease) ~ sex + food, family = binomial, data = babyfood)
round(cbind(OR = exp(coef(bf)), exp(confint(bf))), 3)
c(deviance = deviance(bf), df = df.residual(bf), p = pchisq(deviance(bf), df.residual(bf), lower.tail = FALSE))

x <- 1:6; y <- c(0, 0, 0, 1, 1, 1)                               # complete separation at x = 3.5
fsep <- suppressWarnings(glm(y ~ x, family = binomial))
c(coef(fsep), se = coef(summary(fsep))[2, 2], deviance = deviance(fsep), iterations = fsep$iter)
x2 <- c(1, 2, 3, 3, 4, 5, 6); y2 <- c(0, 0, 0, 1, 1, 1, 1)       # quasi-complete: both outcomes at x = 3
round(coef(summary(suppressWarnings(glm(y2 ~ x2, family = binomial)))), 4)

fb <- glm(y ~ x, family = binomial, method = brglm2::brglmFit)
round(coef(summary(fb)), 3)

fp <- glm(type ~ npreg + glu + bmi + ped + age, family = binomial, data = Pima.tr)
round(coef(summary(fp)), 4)
p_tr <- fitted(fp); y_tr <- Pima.tr$type == "Yes"
p_te <- predict(fp, Pima.te, type = "response"); y_te <- Pima.te$type == "Yes"

auc <- function(p, y) (sum(rank(p)[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
c(train = auc(p_tr, y_tr), test = auc(p_te, y_te),
  pairs = mean(outer(p_te[y_te], p_te[!y_te], ">")))
rates <- function(p, y, cut) { yhat <- p >= cut; c(sens = mean(yhat[y]), spec = mean(!yhat[!y])) }
conf <- function(p, y, cut) { yhat <- p >= cut; c(TP = sum(yhat & y), FP = sum(yhat & !y), FN = sum(!yhat & y), TN = sum(!yhat & !y)) }
rbind(`c = 0.50` = conf(p_te, y_te, 0.5), `c = 0.27` = conf(p_te, y_te, 0.27))

b <- cut(p_te, c(0, 0.1, 0.25, 0.5, 0.75, 1))
round(data.frame(n = as.vector(table(b)), mean_pred = tapply(p_te, b, mean), observed = tapply(y_te, b, mean)), 3)
cost <- function(p, y, cut, cfp = 1, cfn = 3) { yh <- p >= cut; cfp * sum(yh & !y) + cfn * sum(!yh & y) }
sapply(c(`c = 0.50` = 0.5, `c = 0.27` = 0.27, `c* = 0.25` = 0.25), function(k) cost(p_te, y_te, k))

cuts <- seq(0.05, 0.95, by = 0.01)
J <- function(p, y) sapply(cuts, function(k) sum(rates(p, y, k)) - 1)
c_train <- cuts[which.max(J(p_tr, y_tr))]                          # honest: training data only
c_test  <- cuts[which.max(J(p_te, y_te))]                          # LEAKAGE: uses the test labels
c(c_train = c_train, c_test = c_test)
rbind(`chosen on training` = rates(p_te, y_te, c_train), `chosen on test (leaky)` = rates(p_te, y_te, c_test))

al <- rbind(Pima.tr, Pima.te)
set.seed(2027)
res <- t(replicate(200, {
  idx <- sample(nrow(al), 200); tr <- al[idx, ]; te <- al[-idx, ]
  f <- glm(type ~ npreg + glu + bmi + ped + age, family = binomial, data = tr)
  ptr <- fitted(f); pte <- predict(f, te, type = "response"); ytr <- tr$type == "Yes"; yte <- te$type == "Yes"
  c(honest = J(pte, yte)[which.max(J(ptr, ytr))],
    leaked = max(J(pte, yte)))
}))
round(rbind(mean = colMeans(res), sd = apply(res, 2, sd)), 3)
mean(res[, "leaked"] - res[, "honest"])

set.seed(2027); K <- 5; fold <- sample(rep(1:K, length.out = nrow(Pima.tr)))
pcv <- numeric(nrow(Pima.tr))
for (k in 1:K) {
  f <- glm(type ~ npreg + glu + bmi + ped + age, family = binomial, data = Pima.tr[fold != k, ])
  pcv[fold == k] <- predict(f, Pima.tr[fold == k, ], type = "response")
}
c_cv <- cuts[which.max(J(pcv, y_tr))]
c(c_cv = c_cv); rates(p_te, y_te, c_cv)

sh <- subset(ships, service > 0)
sh$year <- factor(sh$year); sh$period <- factor(sh$period)
fs <- glm(incidents ~ type + year + period + offset(log(service)), family = poisson, data = sh)
round(exp(cbind(RR = coef(fs), confint(fs))), 3)
c(deviance = deviance(fs), df = df.residual(fs), Pearson = sum(residuals(fs, "pearson")^2))

ci <- confint(fs); rr <- data.frame(term = names(coef(fs)), est = exp(coef(fs)), lo = exp(ci[, 1]), hi = exp(ci[, 2]))[-1, ]
rr$term <- factor(rr$term, levels = rev(rr$term))
ggplot(rr, aes(est, term)) + geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, colour = ie_navy, orientation = "y") +
  geom_point(colour = ie_coral, size = 2.5) + scale_x_log10() + labs(x = "rate ratio", y = NULL)

nd <- data.frame(type = "B", year = "65", period = "75", service = c(1000, 2000))
predict(fs, nd, type = "response")

fs_free <- glm(incidents ~ type + year + period + log(service), family = poisson, data = sh)
round(coef(summary(fs_free))["log(service)", ], 3)

ins <- Insurance; ins$Group <- factor(ins$Group, ordered = FALSE); ins$Age <- factor(ins$Age, ordered = FALSE)
fi <- glm(Claims ~ District + Group + Age + offset(log(Holders)), family = poisson, data = ins)
round(exp(cbind(RR = coef(fi), confint(fi))), 3)[c("Age>35", "Group>2l", "District4"), ]
c(phi_hat = sum(residuals(fi, "pearson")^2) / df.residual(fi), deviance = deviance(fi), df = df.residual(fi))

fp <- glm(Days ~ Eth + Sex + Age + Lrn, family = poisson, data = quine)
phi_hat <- sum(residuals(fp, "pearson")^2) / df.residual(fp)
mu <- fitted(fp); yq <- quine$Days
aux <- lm(((yq - mu)^2 - yq) / mu ~ 0 + mu)
c(phi_pearson = phi_hat, phi_deviance = deviance(fp) / df.residual(fp))
ct <- coef(summary(aux)); c(alpha = ct[1, 1], t = ct[1, 3], p_one_sided = pt(ct[1, 3], df.residual(aux), lower.tail = FALSE))

fq  <- update(fp, family = quasipoisson)
fnb <- glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine)
c(phi = summary(fq)$dispersion, k = fnb$theta, se_k = fnb$SE.theta)
tab <- cbind(est_poisson = coef(fp), se_poisson = coef(summary(fp))[, 2], se_quasi = coef(summary(fq))[, 2],
             est_nb = coef(fnb), se_nb = coef(summary(fnb))[, 2])
round(tab, 3)
c(AIC_poisson = AIC(fp), AIC_nb = AIC(fnb))

anova(update(fq, . ~ . - Lrn), fq, test = "F")
anova(update(fnb, . ~ . - Lrn), fnb)

cells <- aggregate(Days ~ Eth + Sex + Age + Lrn, data = quine, FUN = function(d) c(m = mean(d), v = var(d), n = length(d)))
cells <- data.frame(cells[, 1:4], cells$Days); cells <- subset(cells, n >= 3)
mm <- exp(seq(log(min(cells$m) * 0.9), log(max(cells$m) * 1.1), length.out = 100))
curves <- rbind(data.frame(model = "Poisson", m = mm, v = mm),
                data.frame(model = "quasi-Poisson", m = mm, v = summary(fq)$dispersion * mm),
                data.frame(model = "negative binomial", m = mm, v = mm + mm^2 / fnb$theta))
ggplot(cells, aes(m, v)) + geom_line(data = curves, aes(colour = model), linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = ie_pal[c(3, 2, 1)]) + scale_x_log10() + scale_y_log10() +
  labs(x = "cell mean", y = "cell variance", colour = NULL)

rqr <- function(y, F, ...) {
  lo <- F(y - 1, ...); hi <- F(y, ...); u <- runif(length(y))
  up_lo <- F(y, ..., lower.tail = FALSE); up_hi <- F(y - 1, ..., lower.tail = FALSE)
  ifelse(lo < 0.5, qnorm(lo + u * (hi - lo)), qnorm(up_lo + u * (up_hi - up_lo), lower.tail = FALSE))
}
set.seed(2027)
r_p <- rqr(quine$Days, ppois, lambda = fitted(fp))
r_n <- rqr(quine$Days, pnbinom, size = fnb$theta, mu = fitted(fnb))
d <- rbind(data.frame(model = "Poisson", lm = log(fitted(fp)), r = r_p), data.frame(model = "negative binomial", lm = log(fitted(fnb)), r = r_n))
ggplot(d, aes(lm, r)) + geom_hline(yintercept = c(-3, 3), linetype = 2, colour = "grey60") + geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(colour = ie_navy, alpha = 0.7) + facet_wrap(~model) + labs(x = expression(log(hat(mu))), y = "randomized quantile residual")
round(rbind(Poisson = range(r_p), NB = range(r_n)), 2)

c(observed = sum(quine$Days == 0), Poisson = sum(dpois(0, fitted(fp))),
  NB = sum(dnbinom(0, size = fnb$theta, mu = fitted(fnb))))

y <- dvisits$doctorco
nll <- function(par) { pi <- plogis(par[1]); mu <- exp(par[2])
  -sum(ifelse(y == 0, log(pi + (1 - pi) * exp(-mu)), log(1 - pi) + dpois(y, mu, log = TRUE))) }
o <- optim(c(0, 0), nll)
z <- zeroinfl(y ~ 1 | 1, dist = "poisson")
rbind(by_hand = c(pi = plogis(o$par[1]), mu = exp(o$par[2]), loglik = -o$value),
      pscl = c(pi = plogis(coef(z)[2]), mu = exp(coef(z)[1]), loglik = as.numeric(logLik(z))))
c(mean = mean(y), variance = var(y), zeros = mean(y == 0), poisson_zeros = dpois(0, mean(y)))

fm <- doctorco ~ sex + age + income + levyplus + freepoor + illness + actdays + hscore
m1 <- glm(fm, family = poisson, data = dvisits)
m2 <- glm.nb(fm, data = dvisits)
m3 <- zeroinfl(fm, data = dvisits)                        # same covariates in both parts
m4 <- zeroinfl(fm, data = dvisits, dist = "negbin")
m5 <- hurdle(fm, data = dvisits)
m6 <- hurdle(fm, data = dvisits, dist = "negbin")
ez <- c(Poisson = sum(dpois(0, fitted(m1))), NB = sum(dnbinom(0, size = m2$theta, mu = fitted(m2))),
        ZIP = sum(predict(m3, type = "prob")[, 1]), ZINB = sum(predict(m4, type = "prob")[, 1]),
        Hurdle = sum(predict(m5, type = "prob")[, 1]), HurdleNB = sum(predict(m6, type = "prob")[, 1]))
round(data.frame(df = sapply(list(m1, m2, m3, m4, m5, m6), function(m) attr(logLik(m), "df")),
                 AIC = sapply(list(m1, m2, m3, m4, m5, m6), AIC), expected_zeros = ez,
                 row.names = names(ez)), 1)
sum(dvisits$doctorco == 0)

tb <- matrix(c(30, 44, 29, 86), 2, byrow = TRUE)
f1 <- glm(low ~ smoke, family = binomial, data = birthwt)
c(by_hand = log((30 * 86) / (44 * 29)), glm = unname(coef(f1)[2]))
c(se_hand = sqrt(sum(1 / tb)), se_glm = coef(summary(f1))[2, 2])
c(risk_ratio = (30 / 74) / (29 / 115))

(77 / 381) / (47 / 447)                         # bottle vs breast, boys: crude OR
exp(-coef(bf)["foodBreast"])                    # model: bottle vs breast, same sex
anova(update(bf, . ~ . - food), bf, test = "Chisq")
bf2 <- glm(cbind(disease, nondisease) ~ sex * food, family = binomial, data = babyfood)
anova(bf, bf2, test = "Chisq")

ell <- function(t) -sum(log(1 + exp(-t * abs(1:6 - 3.5))))
sapply(c(0, 1, 2, 5, 10), ell)

sc <- c(0.9, 0.8, 0.6, 0.5, 0.3, 0.2); lab <- c(1, 0, 1, 1, 0, 0) == 1
auc(sc, lab)
sapply(c(0.55, 0.4, 0.85), function(k) sum(rates(sc, lab, k)) - 1)

d <- data.frame(events = c(10, 25), py = c(500, 1000), grp = c("A", "B"))
f <- glm(events ~ grp + offset(log(py)), family = poisson, data = d)
c(exp(coef(f)), se = coef(summary(f))[2, 2])
exp(confint.default(f)[2, ])

yy <- c(0, 3, 1, 8); m <- c(1, 2, 2, 5); z <- ((yy - m)^2 - yy) / m
c(z, alpha = sum(m * z) / sum(m^2), phi = sum((yy - m)^2 / m) / 3)

drop1(fi, test = "Chisq")
predict(fi, data.frame(District = "1", Group = "<1l", Age = "<25", Holders = 1000), type = "response")
fi_free <- glm(Claims ~ District + Group + Age + log(Holders), family = poisson, data = ins)   # no offset
b <- coef(summary(fi_free))["log(Holders)", ]
c(exponent = unname(b[1]), se = unname(b[2]), z_vs_1 = unname((b[1] - 1) / b[2]))

c(zeros_obs = sum(dvisits$doctorco == 0), pois = sum(dpois(0, fitted(m1))), nb = sum(dnbinom(0, size = m2$theta, mu = fitted(m2))))
c(phi_pearson = sum(residuals(m1, "pearson")^2) / df.residual(m1))
m2s <- glm.nb(doctorco ~ illness + actdays + hscore, data = dvisits)
c(AIC_full = AIC(m2), AIC_small = AIC(m2s))
