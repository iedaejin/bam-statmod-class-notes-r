source("_common.R")
suppressPackageStartupMessages({library(lme4); library(nlme)})
data(sleepstudy, package = "lme4"); data(Dyestuff, package = "lme4")

m_ri_reml <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
m_ri_ml   <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)
m_dy_reml <- lmer(Yield ~ 1 + (1 | Batch), Dyestuff)
m_dy_ml   <- lmer(Yield ~ 1 + (1 | Batch), Dyestuff, REML = FALSE)
rbind(sleep_REML = as.data.frame(VarCorr(m_ri_reml))$vcov, sleep_ML = as.data.frame(VarCorr(m_ri_ml))$vcov,
      dye_REML = as.data.frame(VarCorr(m_dy_reml))$vcov,   dye_ML = as.data.frame(VarCorr(m_dy_ml))$vcov) |>
  `colnames<-`(c("sigma_b^2", "sigma^2")) |> round(1)

y <- sleepstudy$Reaction; X <- model.matrix(~ Days, sleepstudy)
Z <- model.matrix(~ 0 + Subject, sleepstudy); n <- nrow(X); p <- ncol(X)
V_of  <- function(sb2, s2) sb2 * tcrossprod(Z) + s2 * diag(n)
gls   <- function(V) solve(t(X) %*% solve(V) %*% X, t(X) %*% solve(V) %*% y)
lP <- function(par) {                                   # profile log-likelihood, par = log-variances
  V <- V_of(exp(par[1]), exp(par[2])); r <- y - X %*% gls(V)
  -0.5 * (n * log(2 * pi) + determinant(V)$modulus + t(r) %*% solve(V, r))[1]
}
opt <- optim(c(log(1000), log(1000)), function(p) -lP(p), control = list(reltol = 1e-12))
rbind(by_hand = exp(opt$par), lme4 = as.data.frame(VarCorr(m_ri_ml))$vcov)
Vh <- V_of(exp(opt$par[1]), exp(opt$par[2]))
cbind(by_hand = drop(gls(Vh)), lme4 = fixef(m_ri_ml))
cbind(by_hand = sqrt(diag(solve(t(X) %*% solve(Vh) %*% X))), lme4 = sqrt(diag(vcov(m_ri_ml))))

lR <- function(sb2, s2) {
  V <- V_of(sb2, s2); Vi <- solve(V); XtViX <- t(X) %*% Vi %*% X
  r <- y - X %*% solve(XtViX, t(X) %*% Vi %*% y)
  -0.5 * ((n - p) * log(2 * pi) + determinant(V)$modulus + determinant(XtViX)$modulus + t(r) %*% Vi %*% r)[1]
}
vc <- as.data.frame(VarCorr(m_ri_reml))$vcov
c(by_hand = -2 * lR(vc[1], vc[2]), REMLcrit = REMLcrit(m_ri_reml))
opt_r <- optim(c(log(1000), log(1000)), function(p) -lR(exp(p[1]), exp(p[2])), control = list(reltol = 1e-12))
rbind(REML_by_hand = exp(opt_r$par), REML_lme4 = vc)

aov_t <- summary(aov(Yield ~ Batch, Dyestuff))[[1]]
MSB <- aov_t[1, "Mean Sq"]; MSW <- aov_t[2, "Mean Sq"]; m <- 6; nn <- 5
rbind(REML_formula = c(MSB - MSW, MSW * nn) / c(nn, nn), ML_formula = c((m - 1) / m * MSB - MSW, MSW * nn) / c(nn, nn),
      REML_lme4 = as.data.frame(VarCorr(m_dy_reml))$vcov, ML_lme4 = as.data.frame(VarCorr(m_dy_ml))$vcov)

set.seed(2027)
g <- factor(rep(1:6, each = 4)); xx <- rep(1:4, 6)
est <- t(replicate(500, {
  yy <- 1 + 0.5 * xx + rnorm(6)[g] + rnorm(24); d <- data.frame(yy, xx, g)
  c(ML   = as.data.frame(VarCorr(suppressMessages(lmer(yy ~ xx + (1 | g), d, REML = FALSE))))$vcov[1],
    REML = as.data.frame(VarCorr(suppressMessages(lmer(yy ~ xx + (1 | g), d))))$vcov[1])
}))
round(colMeans(est), 3)
d_est <- data.frame(v = c(est), method = rep(colnames(est), each = nrow(est)))
ggplot(d_est, aes(v, fill = method)) + geom_histogram(bins = 30, alpha = .6, position = "identity") +
  geom_vline(xintercept = 1, linetype = 2) + scale_fill_manual(values = c(ML = ie_coral, REML = ie_navy)) +
  labs(x = expression(hat(sigma)[b]^2), y = "count", fill = NULL)

vc_d <- as.data.frame(VarCorr(m_dy_reml))$vcov
k_dye <- 5 * vc_d[1] / (5 * vc_d[1] + vc_d[2])
raw_dye <- tapply(Dyestuff$Yield, Dyestuff$Batch, mean) - fixef(m_dy_reml)
round(rbind(raw = raw_dye, k_times_raw = k_dye * raw_dye, blup = ranef(m_dy_reml)$Batch[, 1]), 2)
k_dye

sl_unb <- subset(sleepstudy, !(Subject %in% c("308", "309", "310")) | Days <= 1)
m_unb <- lmer(Reaction ~ Days + (1 | Subject), sl_unb)
vcu <- as.data.frame(VarCorr(m_unb))$vcov
ni <- as.vector(table(sl_unb$Subject))
raw_u <- tapply(sl_unb$Reaction - model.matrix(~ Days, sl_unb) %*% fixef(m_unb), sl_unb$Subject, mean)
du <- data.frame(raw = as.numeric(raw_u), blup = ranef(m_unb)$Subject[, 1], n = factor(ni))
ku <- unique(round(ni * vcu[1] / (ni * vcu[1] + vcu[2]), 3)); ku
ggplot(du, aes(raw, blup, colour = n)) + geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_abline(slope = ku[1], intercept = 0, linetype = 2) + geom_abline(slope = ku[2], intercept = 0, linetype = 2) +
  geom_point(size = 2.3) + scale_colour_manual(values = c(ie_coral, ie_navy), name = "observations per subject") +
  labs(x = "raw deviation of the subject mean (ms)", y = "BLUP (ms)")

m_pp <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
G <- unclass(VarCorr(m_pp)$Subject)[1:2, 1:2]; s2 <- sigma(m_pp)^2
d335 <- subset(sleepstudy, Subject == "335"); Zi <- cbind(1, d335$Days)
b_ols <- solve(crossprod(Zi), crossprod(Zi, d335$Reaction - Zi %*% fixef(m_pp)))
K <- solve(crossprod(Zi) + s2 * solve(G), crossprod(Zi))
cbind(K_times_ols = drop(K %*% b_ols), ranef = unlist(ranef(m_pp)$Subject["335", ]))
round(K, 3)

cf <- summary(m_pp)$coefficients; cf
c(t = cf["Days", "t value"], p_normal = 2 * pnorm(-abs(cf["Days", "t value"])))
m1_ml <- update(m_pp, REML = FALSE); m0_ml <- update(m1_ml, . ~ . - Days)
anova(m0_ml, m1_ml)                       # LRT for Days, ML fits
2 * (logLik(m_pp) - logLik(update(m_pp, . ~ . - Days)))   # NOT a valid LRT: two REML fits with different X
confint(update(m_pp), parm = "Days", method = "profile")     # on a copy: profiling changes the internal state of a merMod fit

fit_t <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(fit_t)$coefficients                                   # Satterthwaite degrees of freedom
pbkrtest::KRmodcomp(m_pp, update(m_pp, . ~ . - Days))$stats    # Kenward-Roger F test

confint(update(m_dy_reml), method = "profile")      # Dyestuff: only 6 batches
confint(update(m_pp), method = "profile")           # sleepstudy: sig01 (sd intercept), sig02 (sd slope), sig03 (corr)

m_lm <- lm(Reaction ~ Days, sleepstudy)
T_int <- as.numeric(2 * (logLik(m_ri_ml) - logLik(m_lm)))                   # random intercept vs none (ML)
c(T = T_int, p_naive = pchisq(T_int, 1, lower.tail = FALSE), p_mixture = .5 * pchisq(T_int, 1, lower.tail = FALSE))
a <- anova(m_ri_reml, m_pp, refit = FALSE)                                  # random slope + corr vs intercept only
a
T_sl <- a$Chisq[2]
c(p_naive = pchisq(T_sl, 2, lower.tail = FALSE), p_mixture = .5 * pchisq(T_sl, 1, lower.tail = FALSE) + .5 * pchisq(T_sl, 2, lower.tail = FALSE))

set.seed(2027)
g30 <- factor(rep(1:30, each = 5)); x30 <- rep(1:5, 30)
T_sim <- replicate(500, {
  ys <- 2 + .3 * x30 + rnorm(150); d <- data.frame(ys, x30, g30)
  f1 <- suppressMessages(lmer(ys ~ x30 + (1 | g30), d, REML = FALSE))
  max(0, as.numeric(2 * (logLik(f1) - logLik(lm(ys ~ x30, d)))))
})
c(share_zero = mean(T_sim < 1e-6), q95 = unname(quantile(T_sim, .95)), q95_mixture = qchisq(.90, 1), q95_naive = qchisq(.95, 1))
ggplot(data.frame(T = T_sim), aes(T)) + geom_histogram(bins = 40, fill = ie_blue, alpha = .7) +
  geom_vline(xintercept = c(qchisq(.9, 1), qchisq(.95, 1)), linetype = 2) + labs(x = "LRT statistic", y = "count")

knitr::kable(data.frame(
  Task = c("Report variance components", "Compare models with different fixed effects (LRT, AIC, BIC)",
           "Compare models with the same fixed effects, different random effects"),
  Use = c("REML (default in lmer and lme)", "ML (REML = FALSE); anova() refits with ML automatically", "REML or ML, the same for both models")),
  booktabs = TRUE, caption = "Choosing between ML and REML.") |> kableExtra::column_spec(1:2, width = "7cm")

m_unc <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)       # rho = 0
data.frame(model = c("(1 | Subject)", "(Days || Subject), rho = 0", "(Days | Subject)"),
           npar = c(4, 5, 6), AIC = AIC(m_ri_reml, m_unc, m_pp)$AIC, BIC = BIC(m_ri_reml, m_unc, m_pp)$BIC) |> 
  transform(AIC = round(AIC, 1), BIC = round(BIC, 1))
anova(m_unc, m_pp, refit = FALSE)

res <- resid(m_pp); fit <- fitted(m_pp); re <- ranef(m_pp)$Subject
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(fit, res, xlab = "fitted", ylab = "conditional residual", main = "Residuals vs fitted"); abline(h = 0, lty = 2)
qqnorm(res, main = "Residuals"); qqline(res)
qqnorm(re[, 1], main = "BLUPs: intercepts"); qqline(re[, 1])
qqnorm(re[, 2], main = "BLUPs: slopes"); qqline(re[, 2])
par(op)
head(sleepstudy[order(-abs(resid(m_pp, scaled = TRUE))), ], 3)

Gm <- kronecker(diag(18), G); Zm <- as.matrix(getME(m_pp, "Z")); Xm <- getME(m_pp, "X")
Vm <- Zm %*% Gm %*% t(Zm) + s2 * diag(180); Vi <- solve(Vm)
Pm <- Vi - Vi %*% Xm %*% solve(t(Xm) %*% Vi %*% Xm, t(Xm) %*% Vi)
c(max_abs_diff = max(abs(resid(m_pp) - s2 * drop(Pm %*% sleepstudy$Reaction))))
summary(sqrt(s2^2 * diag(Pm)) / sqrt(s2))      # sd(e_j)/sigma: below 1

# fit_reml <- lmer(y ~ x + (x | g), data = d)               # REML (default)
# fit_ml   <- update(fit_reml, REML = FALSE)                # ML, for LRT / AIC with different fixed effects
# VarCorr(fit_reml); fixef(fit_reml); vcov(fit_reml)        # variance components, beta-hat, Var(beta-hat)
# ranef(fit_reml, condVar = TRUE)                           # BLUPs with conditional variances
# anova(fit_small, fit_big)                                 # LRT; refits with ML automatically
# anova(fit_a, fit_b, refit = FALSE)                        # keep REML: same fixed effects
# confint(fit_reml, method = "profile")                     # profile intervals (also "boot")
# bootMer(fit_reml, FUN = fixef, nsim = 500)                # parametric bootstrap
# resid(fit_reml, scaled = TRUE); plot(fit_reml)            # level-1 residuals
# lme(y ~ x, random = ~ x | g, data = d, method = "REML")  # same model in nlme; intervals(fit)

set.seed(11); n8 <- 8
W <- matrix(rnorm(n8^2), n8); V8 <- crossprod(W) + diag(n8)        # random positive definite V
X8 <- cbind(1, matrix(rnorm(n8 * 2), n8))
A8 <- qr.Q(qr(X8), complete = TRUE)[, -(1:3)]                      # orthonormal basis of the complement of col(X)
max(abs(t(A8) %*% X8))                                              # A'X = 0
Vi8 <- solve(V8); P8 <- Vi8 - Vi8 %*% X8 %*% solve(t(X8) %*% Vi8 %*% X8, t(X8) %*% Vi8)
max(abs(A8 %*% solve(t(A8) %*% V8 %*% A8) %*% t(A8) - P8))          # identity (i)
c(det_left = determinant(t(A8) %*% V8 %*% A8)$modulus,
  det_right = determinant(V8)$modulus + determinant(t(X8) %*% Vi8 %*% X8)$modulus - determinant(crossprod(X8))$modulus)  # (ii) on log scale

es <- as.data.frame(nlme::ergoStool); es$Subject <- factor(as.character(es$Subject))
f1 <- lmer(effort ~ Type + (1 | Subject), es); f0 <- lmer(effort ~ 1 + (1 | Subject), es)
a_ml <- anova(f0, f1)                                        # refits with ML
T_ml <- a_ml$Chisq[2]
T_reml <- 2 * (logLik(f1) - logLik(f0))                       # wrong: REML with different X
c(T_ML = T_ml, p_ML = a_ml$`Pr(>Chisq)`[2], T_REML_wrong = as.numeric(T_reml))
aov_es <- summary(aov(effort ~ Type + Error(Subject), es)); aov_es
F_es <- aov_es[[2]][[1]]['Type', 'F value']

ox <- as.data.frame(nlme::Oxboys); ox$Subject <- factor(as.character(ox$Subject))
a_ox <- anova(lmer(height ~ age + (1 | Subject), ox), lmer(height ~ age + (age | Subject), ox), refit = FALSE)
T_ox <- a_ox$Chisq[2]
c(T = T_ox, p_naive = pchisq(T_ox, 2, lower.tail = FALSE),
  p_mixture = .5 * pchisq(T_ox, 1, lower.tail = FALSE) + .5 * pchisq(T_ox, 2, lower.tail = FALSE))

set.seed(5)
lm0 <- lm(Yield ~ 1, Dyestuff)
T_obs <- as.numeric(2 * (logLik(m_dy_ml) - logLik(lm0)))
B <- 600
T_b <- replicate(B, {
  d <- Dyestuff; d$Yield <- simulate(lm0)[[1]]
  max(0, as.numeric(2 * (logLik(suppressMessages(lmer(Yield ~ 1 + (1 | Batch), d, REML = FALSE))) - logLik(lm(Yield ~ 1, d)))))
})
c(T_obs = T_obs, p_boot = (1 + sum(T_b >= T_obs)) / (B + 1), p_mixture = .5 * pchisq(T_obs, 1, lower.tail = FALSE),
  p_naive = pchisq(T_obs, 1, lower.tail = FALSE), P_T_equals_0 = mean(T_b < 1e-6))

x0 <- c(1, 10); fe <- fixef(m_pp); Vb <- as.matrix(vcov(m_pp))
G <- unclass(VarCorr(m_pp)$Subject)[1:2, 1:2]; s2 <- sigma(m_pp)^2
mu0 <- sum(x0 * fe)
v_mean <- drop(t(x0) %*% Vb %*% x0)                  # uncertainty of the population mean
v_new  <- v_mean + drop(t(x0) %*% G %*% x0) + s2     # formula (c12-pred), z0 = x0
rbind(mean_CI = mu0 + c(-1.96, 1.96) * sqrt(v_mean), new_subject_PI = mu0 + c(-1.96, 1.96) * sqrt(v_new))
c(sd_mean = sqrt(v_mean), sd_new = sqrt(v_new))
