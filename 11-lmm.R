source("_common.R")
suppressPackageStartupMessages({library(lme4); library(nlme)})

data(sleepstudy, package = "lme4")
ggplot(sleepstudy, aes(Days, Reaction, group = Subject)) +
  geom_line(colour = ie_navy, alpha = .45) + geom_point(colour = ie_navy, size = .8, alpha = .6) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = ie_coral, linewidth = 1.2) +
  scale_x_continuous(breaks = 0:9) + labs(x = "days of sleep restriction", y = "reaction time (ms)")

fit_cp <- lm(Reaction ~ Days, data = sleepstudy)
r  <- matrix(resid(fit_cp), nrow = 10)      # column = subject, row = day
cr <- cor(t(r))                             # 10 x 10 day-by-day correlations over the 18 subjects
mean(cr[upper.tri(cr)])
summary(fit_cp)$coefficients

fit_np <- lme4::lmList(Reaction ~ Days | Subject, data = sleepstudy)       # no pooling
fit_pp <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)      # partial pooling
sub <- c("330", "335", "370")
tab <- cbind(coef(fit_np)[sub, ], coef(fit_pp)$Subject[sub, ])
colnames(tab) <- c("int (no pool)", "slope (no pool)", "int (partial)", "slope (partial)")
round(tab, 2)
rbind(complete = coef(fit_cp), population = fixef(fit_pp))

subs <- c("308", "309", "330", "335", "350", "370")
ln <- rbind(
  data.frame(Subject = rownames(coef(fit_np)), a = coef(fit_np)[, 1], b = coef(fit_np)[, 2], fit = "no pooling"),
  data.frame(Subject = rownames(coef(fit_pp)$Subject), a = coef(fit_pp)$Subject[, 1],
             b = coef(fit_pp)$Subject[, 2], fit = "partial pooling"),
  data.frame(Subject = levels(sleepstudy$Subject), a = coef(fit_cp)[[1]], b = coef(fit_cp)[[2]], fit = "complete pooling"))
ln <- ln[ln$Subject %in% subs, ]
ggplot(sleepstudy[sleepstudy$Subject %in% subs, ], aes(Days, Reaction)) +
  geom_point(colour = "grey30", size = 1) +
  geom_abline(data = ln, aes(intercept = a, slope = b, colour = fit, linetype = fit), linewidth = .8) +
  scale_colour_manual(values = c("complete pooling" = "grey45", "no pooling" = ie_coral, "partial pooling" = ie_navy), name = NULL) +
  scale_linetype_manual(values = c("complete pooling" = "dotted", "no pooling" = "dashed", "partial pooling" = "solid"), name = NULL) +
  facet_wrap(~ Subject, nrow = 2) + labs(x = "days", y = "reaction time (ms)")

knitr::kable(data.frame(
  Random = c("levels are a sample from a larger population (subjects, herds, schools)",
             "we want conclusions beyond the observed levels",
             "many levels (say 5 or more)",
             "group effects plausibly unrelated to the covariates"),
  Fixed = c("levels are the whole population of interest (the 3 fertilisers tested)",
            "we care about these specific levels",
            "very few levels",
            "group effects may be related to the covariates (confounding)")),
  booktabs = TRUE, caption = "Guide to fixed versus random group effects.") |>
  kableExtra::column_spec(1:2, width = "7cm")

set.seed(2027)
m <- 60; n <- 5; id <- rep(1:m, each = n)
one_sim <- function() {
  u <- rnorm(m); x <- u[id] + rnorm(m * n)       # covariate partly driven by u
  b <- u + rnorm(m, 0, .5)                       # group effect correlated with the covariate
  y <- 1 * x + b[id] + rnorm(m * n)              # true slope beta = 1
  d <- data.frame(y, x, id = factor(id)); d$xbar <- ave(d$x, d$id)
  c(pooled = coef(lm(y ~ x, d))[[2]],
    random = fixef(lmer(y ~ x + (1 | id), d))[[2]],
    fixed  = coef(lm(y ~ x + id, d))[[2]],
    random_plus_group_mean = fixef(lmer(y ~ x + xbar + (1 | id), d))[[2]])
}
sims <- replicate(200, one_sim())
round(rbind(mean = rowMeans(sims), sd = apply(sims, 1, sd)), 3)

data(Dyestuff, package = "lme4")
tab_aov <- summary(aov(Yield ~ Batch, Dyestuff))[[1]]
tab_aov
MSB <- tab_aov[1, "Mean Sq"]; MSW <- tab_aov[2, "Mean Sq"]; n_s <- 5
c(sigma2 = MSW, sigma_b2 = (MSB - MSW) / n_s)
fit_dy <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff)
VarCorr(fit_dy)
vc_dy <- as.data.frame(VarCorr(fit_dy))$vcov
c(ICC = vc_dy[1] / sum(vc_dy), deff = 1 + (n_s - 1) * vc_dy[1] / sum(vc_dy))
ggplot(Dyestuff, aes(Batch, Yield)) + geom_point(colour = ie_navy, size = 2, alpha = .8) +
  stat_summary(fun = mean, geom = "point", colour = ie_coral, shape = 18, size = 4) +
  geom_hline(yintercept = mean(Dyestuff$Yield), linetype = 2, colour = "grey45")

m_ri <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
summary(m_ri)$coefficients
VarCorr(m_ri)
vc_ri <- as.data.frame(VarCorr(m_ri))$vcov
c(ICC = vc_ri[1] / sum(vc_ri), deff = 1 + 9 * vc_ri[1] / sum(vc_ri))

m_ris <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
print(VarCorr(m_ris), comp = c("Variance", "Std.Dev."))
fixef(m_ris)
G <- unclass(VarCorr(m_ris)$Subject)[1:2, 1:2]
s2 <- sigma(m_ris)^2
Zi <- cbind(1, 0:9)
Vi <- Zi %*% G %*% t(Zi) + s2 * diag(10)         # V_i = Z_i G Z_i' + sigma^2 I
round(diag(Vi))                                     # the variance grows with Days
round(cov2cor(Vi)[c(1, 5, 10), c(1, 5, 10)], 2)    # correlation of Days 0, 4, 9

emp <- aggregate(Reaction ~ Days, sleepstudy, sd)
imp <- data.frame(Days = 0:9, ri = sqrt(sum(vc_ri)), rs = sqrt(diag(Vi)))
ggplot(imp, aes(Days)) +
  geom_line(aes(y = ri, colour = "random intercept"), linewidth = 1) +
  geom_line(aes(y = rs, colour = "intercept + slope"), linewidth = 1) +
  geom_point(data = emp, aes(y = Reaction)) +
  scale_colour_manual(values = c("random intercept" = ie_coral, "intercept + slope" = ie_navy), name = NULL) +
  scale_x_continuous(breaks = 0:9) + labs(y = "SD of reaction time (ms)")

m_c <- lmer(Reaction ~ I(Days - 4.5) + (I(Days - 4.5) | Subject), data = sleepstudy)
print(VarCorr(m_c), comp = c("Variance", "Std.Dev."))
c(logLik_original = as.numeric(logLik(m_ris)), logLik_centred = as.numeric(logLik(m_c)))
-G[1, 2] / G[2, 2]    # centre that makes intercept and slope uncorrelated

X <- getME(m_ris, "X"); Z <- getME(m_ris, "Z")
c(n = nrow(X), p = ncol(X), q = ncol(Z))
as.matrix(Z[1:3, 1:4])           # first three rows: columns = (int_1, slope_1, int_2, slope_2)

m_ml <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
Zm <- as.matrix(getME(m_ml, "Z")); Xm <- getME(m_ml, "X"); ym <- getME(m_ml, "y")
Gm <- kronecker(diag(18), unclass(VarCorr(m_ml)$Subject)[1:2, 1:2])   # b is ordered (int_1, slope_1, int_2, ...)
Vm <- Zm %*% Gm %*% t(Zm) + sigma(m_ml)^2 * diag(180)
res <- ym - Xm %*% fixef(m_ml)
ll_by_hand <- -0.5 * (180 * log(2 * pi) + determinant(Vm)$modulus + t(res) %*% solve(Vm, res))
print(c(by_hand = ll_by_hand[1], lme4 = as.numeric(logLik(m_ml))), digits = 10)

data(Penicillin, package = "lme4")      # diameter of inhibition zone; plate x sample (crossed)
VarCorr(lmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin))
data(Pastes, package = "lme4")          # strength of paste; cask nested in batch
VarCorr(lmer(strength ~ 1 + (1 | batch/cask), Pastes))

g_cs  <- gls(Reaction ~ Days, sleepstudy, correlation = corCompSymm(form = ~ 1 | Subject))
g_ar1 <- gls(Reaction ~ Days, sleepstudy, correlation = corAR1(form = ~ Days | Subject))
rbind(gls_compound_symmetry = c(coef(g_cs), logLik = as.numeric(logLik(g_cs))),
      lmer_random_intercept = c(fixef(m_ri), logLik = as.numeric(logLik(m_ri))))
round(c(rho_cs = coef(g_cs$modelStruct$corStruct, unconstrained = FALSE),
        phi_ar1 = coef(g_ar1$modelStruct$corStruct, unconstrained = FALSE)), 3)
c(gls_cs = AIC(g_cs), gls_ar1 = AIC(g_ar1), lmer_ri = AIC(m_ri), lmer_rs = AIC(m_ris))

tt <- function(x) paste0("\\texttt{", gsub("~", "\\textasciitilde{}", gsub("|", "\\textbar{}", x, fixed = TRUE), fixed = TRUE), "}")
knitr::kable(data.frame(
  `lmer formula` = tt(c("y ~ x + (1 | g)", "y ~ x + (x | g)", "y ~ x + (x || g)", "y ~ x + (1 | g1) + (1 | g2)", "y ~ x + (1 | g1/g2)")),
  Meaning = c("random intercept", "correlated random intercept and slope", "uncorrelated intercept and slope ($\\rho=0$)",
              "crossed random effects", "nested: g2 within g1, i.e. (1 | g1) + (1 | g1:g2)"),
  check.names = FALSE),
  booktabs = TRUE, escape = FALSE, caption = "Random-effects formula syntax in lme4.")

m_lme <- lme(Reaction ~ Days, random = ~ Days | Subject, data = sleepstudy)   # same model in nlme
rbind(lmer = fixef(m_ris), lme = fixef(m_lme))
VarCorr(m_lme)
head(ranef(m_ris)$Subject, 3)            # predicted random effects b~ (BLUPs)
head(coef(m_ris)$Subject, 3)             # beta + b~: subject-specific lines
predict(m_ris, newdata = data.frame(Days = 5, Subject = "308"))          # subject-specific
predict(m_ris, newdata = data.frame(Days = 5, Subject = "new"), allow.new.levels = TRUE)  # population

set.seed(1)
one_trial <- function() {
  m <- 20; n <- 10; g <- rep(1:m, each = n)
  trt <- rep(0:1, each = m / 2)[g]
  y <- rnorm(m, 0, sqrt(.3))[g] + rnorm(m * n, 0, sqrt(.7))
  d <- data.frame(y, trt, g = factor(g))
  c(lm   = summary(lm(y ~ trt, d))$coefficients[2, 4] < .05,
    lmer = abs(summary(lmer(y ~ trt + (1 | g), d))$coefficients[2, 3]) > 1.96)
}
rej <- rowMeans(replicate(500, one_trial()))
rej

es <- as.data.frame(nlme::ergoStool); es$Subject <- factor(as.character(es$Subject))
m_es <- lmer(effort ~ Type + (1 | Subject), es)
summary(m_es)$coefficients
vc_es <- as.data.frame(VarCorr(m_es))$vcov
c(sigma_b2 = vc_es[1], sigma2 = vc_es[2], ICC = vc_es[1] / sum(vc_es))

ox <- as.data.frame(nlme::Oxboys); ox$Subject <- factor(as.character(ox$Subject))
m_ox <- lmer(height ~ age + (age | Subject), ox)
print(VarCorr(m_ox), comp = c("Variance", "Std.Dev."))
fixef(m_ox)
Gx <- unclass(VarCorr(m_ox)$Subject)[1:2, 1:2]; s2x <- sigma(m_ox)^2
Gx[1, 1] / (Gx[1, 1] + s2x)
fixef(m_ox)[2] + c(-1.96, 1.96) * sqrt(Gx[2, 2])

data(Pastes, package = "lme4")
m_pa <- lmer(strength ~ 1 + (1 | batch/cask), Pastes)
vcp <- as.data.frame(VarCorr(m_pa)); vcp[, c("grp", "vcov")]
v <- setNames(vcp$vcov, vcp$grp); round(100 * v / sum(v), 1)

vm <- function(mb, mc, n = 2) v["batch"] / mb + v["cask:batch"] / (mb * mc) + v["Residual"] / (mb * mc * n)
c(`10 batches x 3 casks` = unname(vm(10, 3)), `5 batches x 6 casks` = unname(vm(5, 6)))
