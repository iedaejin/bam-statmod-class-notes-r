source("_common.R")

library(lme4); library(mgcv); library(survival); library(broom)

d <- data.frame(y = c(3, 5, 4, 8, 9, 7), x = 1:6,
                g = factor(c("a", "a", "b", "b", "c", "c")))
model.matrix(y ~ g + x, data = d)

m1 <- lm(mpg ~ wt, data = mtcars)
m2 <- update(m1, . ~ . + hp)          # add a term
formula(m2)
m3 <- update(m2, log(.) ~ .)          # change the response
formula(m3)

m <- lm(mpg ~ wt + hp, data = mtcars)
coef(summary(m))                          # estimate, se, t, p as a matrix
round(confint(m), 3)
# the data actually used (after NA removal)
head(model.frame(m), 3)
predict(m, newdata = data.frame(wt = 3, hp = 120), interval = "prediction")
c(AIC = AIC(m), BIC = BIC(m), logLik = as.numeric(logLik(m)),
  df = df.residual(m), nobs = nobs(m))

tidy(m, conf.int = TRUE)
glance(m)[, c("r.squared", "sigma", "AIC", "BIC", "nobs")]
head(augment(m), 3)[, c(".rownames", "mpg", ".fitted", ".resid", ".hat", ".cooksd")]

# fit <- lm(y ~ x1 + x2 * g, data = d, weights = w, subset = x1 > 0)
# # coefficients, sigma-hat, R^2, F test (all slopes = 0)
# summary(fit)
# # sequential (type I) sums of squares: ORDER MATTERS
# anova(fit)
# anova(fit0, fit)                     # F test of nested models, same rows
# # each term dropped in turn (respects marginality)
# drop1(fit, test = "F")
# car::Anova(fit, type = 2)            # type II tests, order-free
# # residuals vs fitted, QQ, scale-location, Cook's D
# plot(fit, which = 1:4)
# # leverage, influence, studentised
# hatvalues(fit); cooks.distance(fit); rstudent(fit)
# car::vif(fit)                        # variance inflation factors
# # mean;  "prediction": new observation
# predict(fit, newd, interval = "confidence")
# # backward stepwise by BIC (use with care)
# step(fit, k = log(nrow(d)))

# fit <- glm(y ~ x1 + x2, family = binomial(link = "logit"), data = d)
# glm(cbind(succ, fail) ~ x, family = binomial)        # grouped binomial
# glm(count ~ x + offset(log(exposure)), family = poisson)
# # negative binomial (estimates k)
# MASS::glm.nb(count ~ x, data = d)
# # dispersion estimated by Pearson X^2/(n-p)
# glm(y ~ x, family = quasipoisson)
# # Wald z tests, null/residual deviance, AIC, IRLS iterations
# summary(fit)
# # deviance (LR) test of nested models, phi known
# anova(fit0, fit, test = "Chisq")
# # F test when phi is estimated (quasi families)
# anova(fit0, fit, test = "F")
# # odds ratios (or rate ratios) with profile CIs
# exp(cbind(OR = coef(fit), confint(fit)))
# residuals(fit, type = "deviance"); residuals(fit, type = "pearson")
# sum(residuals(fit, "pearson")^2) / df.residual(fit)  # dispersion estimate
# predict(fit, newd, type = "link", se.fit = TRUE)

# REML by default
fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
summary(fm)$coefficients
VarCorr(fm)
head(ranef(fm)$Subject, 3)             # BLUPs of the random effects
# fixed + random: subject-specific intercept and slope
head(coef(fm)$Subject, 3)

# # ML: needed to compare fixed effects
# lmer(y ~ x + (1 | g), data = d, REML = FALSE)
# # refits both with ML before the LR test
# anova(fm0, fm1)
# # keeps REML: valid ONLY if fixed effects are the same
# anova(fm0, fm1, refit = FALSE)
# # profile intervals for fixed effects AND variance components
# confint(fm, method = "profile")
# glmer(cbind(inc, size - inc) ~ period + (1 | herd), family = binomial, data = d)
# # TRUE if a variance is estimated at zero / correlation +-1
# isSingular(fm)
# # alternative with residual correlation structures
# nlme::lme(y ~ x, random = ~ 1 | g, data = d)

aq <- na.omit(airquality[, c("Ozone", "Temp", "Wind", "Solar.R")])
g <- gam(log(Ozone) ~ s(Temp) + s(Wind) + s(Solar.R), data = aq, method = "REML")
# edf, Ref.df, F, p-value for each smooth
summary(g)$s.table
c(AIC = AIC(g), dev_expl = summary(g)$dev.expl, r2 = summary(g)$r.sq)

# gam(y ~ s(x, k = 10, bs = "cr") + z, family = poisson, data = d, method = "REML")
# # tensor product: interaction of smooths on different scales
# gam(y ~ te(x1, x2), data = d)
# # main effects plus pure interaction
# gam(y ~ s(x1) + s(x2) + ti(x1, x2), data = d)
# # a smooth for each level of the factor g
# gam(y ~ s(x, by = g) + g, data = d)
# # random intercept as a smooth (GAMM by gam())
# gam(y ~ s(subject, bs = "re") + s(x), data = d)
# # residual plots, basis-dimension check (k-index, edf, p-value)
# gam.check(g)
# # partial effects with bands
# plot(g, pages = 1, shade = TRUE, residuals = TRUE)
# # non-linear collinearity between smooths
# concurvity(g)
# predict(g, newd, type = "response", se.fit = TRUE)
# # via lme; returns $gam and $lme
# mgcv::gamm(y ~ s(x), random = list(g = ~ 1), data = d)
# # approximate test between nested GAMs
# anova(g0, g1, test = "F")

vet <- survival::veteran
km <- survfit(Surv(time, status) ~ trt, data = vet)
# S(t) at selected times, by arm
summary(km, times = c(90, 180))$surv
cx <- coxph(Surv(time, status) ~ trt + karno + age, data = vet)
round(summary(cx)$coefficients, 4)

# # event: 1 = event, 0 = censored (or TRUE/FALSE)
# Surv(time, event)
# # counting-process format (delayed entry, time-varying x)
# Surv(start, stop, event)
# survfit(Surv(time, event) ~ g, data = d, conf.type = "log-log")
# quantile(km, probs = 0.5)               # median with confidence limits
# # log-rank; rho = 1 Peto-Peto type
# survdiff(Surv(time, event) ~ g, data = d, rho = 0)
# coxph(Surv(time, event) ~ x1 + strata(g) + cluster(id), data = d, ties = "efron")
# # Schoenfeld-residual test of proportional hazards
# cox.zph(fit); plot(cox.zph(fit))
# survSplit(Surv(time, event) ~ ., data = d, cut = c(90, 180), episode = "period")
# coxph(Surv(time, event) ~ x + tt(x), data = d, tt = function(x, t, ...) x * log(t + 20))
# # AFT: coefficients act on log T
# survreg(Surv(time, event) ~ x, data = d, dist = "weibull")
# # predicted survival curves from a Cox model
# survfit(fit, newdata = nd)
# # functional form; "dfbeta" for influence
# residuals(fit, type = "martingale")

# # analysis.R --- question, data, models, results
# library(tidyverse); library(mgcv)           # 1. packages
# set.seed(2027)                              # 2. seed
# # 3. data: raw file is never edited
# dat <- read_csv("data/raw.csv") |>
#   filter(!is.na(y)) |> mutate(g = factor(g))
# stopifnot(nrow(dat) > 100)                  # 4. checks that fail loudly
# fit <- gam(y ~ s(x) + g, data = dat, method = "REML")   # 5. models
# gam.check(fit)                              # 6. diagnostics
# saveRDS(fit, "output/fit.rds")              # 7. results
# sessionInfo()                               # 8. software

msg <- function(expr) {                 # print the message, wrapped to the page width
  m <- tryCatch({ expr; "no error" },
                error = function(e) conditionMessage(e),
                warning = function(w) conditionMessage(w))
  cat(strwrap(m, width = 66), sep = "\n")
}
# 1. different lengths
msg(lm(c(1, 2, 3, 4) ~ c(1, 2, 3)))
# 2. new factor level in predict()
fit_f <- lm(mpg ~ factor(cyl), data = mtcars)
msg(predict(fit_f, newdata = data.frame(cyl = 12)))
# 3. a collinear column
coef(lm(mpg ~ wt + I(2 * wt), data = mtcars))
# 4. models fitted to different rows
a1 <- lm(Ozone ~ Wind, data = airquality)
a2 <- lm(Ozone ~ Wind + Solar.R, data = airquality)
c(nobs(a1), nobs(a2)); msg(anova(a1, a2))
# 5. separation in a logistic model
msg(glm(c(0, 0, 0, 1, 1, 1) ~ c(1, 2, 3, 4, 5, 6), family = binomial))
# 6. basis dimension too large for the data
msg(gam(y ~ s(x, k = 20), data = data.frame(x = rep(1:5, 4), y = {set.seed(2027); rnorm(20)})))
