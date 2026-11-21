source("_common.R")
suppressPackageStartupMessages({library(car); library(ggrepel)})
data(savings, package = "faraway")

fit <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = savings)
n <- nrow(savings); p <- length(coef(fit))
c(n = n, p = p, sigma_hat = sigma(fit))

d <- data.frame(fitted = fitted(fit), res = resid(fit), r = rstandard(fit), country = rownames(savings))
p1 <- ggplot(d, aes(fitted, res)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(colour = ie_navy) + geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ie_coral, linewidth = 0.7) +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs fitted")
p2 <- ggplot(d, aes(sample = r)) + stat_qq_line(colour = ie_coral) + stat_qq(colour = ie_navy) +
  labs(x = "N(0,1) quantiles", y = "Standardized residuals", title = "Normal QQ plot")
p3 <- ggplot(d, aes(fitted, sqrt(abs(r)))) + geom_point(colour = ie_navy) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ie_coral, linewidth = 0.7) +
  labs(x = "Fitted values", y = expression(sqrt(abs(r[i]))), title = "Scale-location")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

h <- hatvalues(fit)
c(sum_h = sum(h), average = mean(h), cutoff_2p_over_n = 2 * p / n)
round(sort(h, decreasing = TRUE)[1:5], 3)
savings["Libya", ]

ti <- rstudent(fit)
cutoff <- qt(1 - 0.05 / (2 * n), df = n - p - 1)
c(max_abs_t = max(abs(ti)), bonferroni_cutoff = cutoff)
round(sort(abs(ti), decreasing = TRUE)[1:3], 3)
outlierTest(fit)
# the identities of Lemma and Theorem, checked directly
i <- which(rownames(savings) == "Zambia")
ri <- rstandard(fit)[i]; hi <- h[i]; ei <- resid(fit)[i]
loo <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = savings[-i, ])
c(t_i = ti[i], from_r = ri * sqrt((n - p - 1) / (n - p - ri^2)),
  loo_prediction_error = savings$sr[i] - predict(loo, savings[i, ]), e_over_1mh = ei / (1 - hi))

dd <- data.frame(h = hatvalues(fit), t = rstudent(fit), D = cooks.distance(fit), country = rownames(savings))
ggplot(dd, aes(h, t)) +
  geom_vline(xintercept = 2 * p / n, linetype = 2, colour = "grey50") +
  geom_hline(yintercept = c(-1, 1) * cutoff, linetype = 3, colour = ie_coral) +
  geom_point(aes(size = D), colour = ie_navy, alpha = 0.7) +
  geom_text_repel(data = subset(dd, h > 2 * p / n | D > 4 / n | abs(t) > 2.5), aes(label = country),
                  colour = ie_coral, size = 3.2, seed = 1) +
  scale_size_area(max_size = 8, name = "Cook's D") +
  labs(x = "Leverage", y = "Studentized residual")
round(sort(cooks.distance(fit), decreasing = TRUE)[1:4], 3)
qf(0.5, p, n - p)

fit_noL <- update(fit, subset = rownames(savings) != "Libya")
tab <- cbind(all = coef(fit), without_Libya = coef(fit_noL),
             diff_in_se = (coef(fit_noL) - coef(fit)) / summary(fit)$coefficients[, "Std. Error"])
round(tab, 4)

data(gala, package = "faraway")
g <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, data = gala)
round(sort(cooks.distance(g), decreasing = TRUE)[1:3], 2)
round(sort(hatvalues(g), decreasing = TRUE)[1:3], 2)

set.seed(606)
x <- runif(80, 0, 10); y <- 2 + 0.4 * x + 0.25 * (x - 5)^2 + rnorm(80, sd = 1)
m_lin <- lm(y ~ x); m_quad <- lm(y ~ x + I(x^2))
q1 <- ggplot(data.frame(x, y), aes(x, y)) + geom_point(colour = ie_navy) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = ie_coral)
q2 <- ggplot(data.frame(f = fitted(m_lin), r = resid(m_lin)), aes(f, r)) + geom_point(colour = ie_navy) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ie_coral) + labs(x = "Fitted", y = "Residual")
gridExtra::grid.arrange(q1, q2, nrow = 1)
anova(m_lin, m_quad)

lake <- data.frame(level = as.numeric(LakeHuron), year = as.numeric(time(LakeHuron)))
lh <- lm(level ~ year, data = lake)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1)); plot(lake$year, resid(lh), type = "b", pch = 16, cex = 0.6,
  col = ie_navy, xlab = "Year", ylab = "Residual", main = "Residuals in time order"); abline(h = 0, lty = 2)
acf(resid(lh), main = "Autocorrelation of residuals"); par(op)
c(lag1_autocorrelation = acf(resid(lh), plot = FALSE)$acf[2])
durbinWatsonTest(lh)

bc <- MASS::boxcox(g, plotit = FALSE)
lambda_hat <- bc$x[which.max(bc$y)]
gs <- lm(sqrt(Species) ~ Area + Elevation + Nearest + Scruz + Adjacent, data = gala)
b1 <- ggplot(data.frame(f = fitted(g), r = resid(g)), aes(f, r)) + geom_point(colour = ie_navy) +
  geom_hline(yintercept = 0, linetype = 2) + labs(x = "Fitted", y = "Residual", title = "Species")
b2 <- ggplot(data.frame(f = fitted(gs), r = resid(gs)), aes(f, r)) + geom_point(colour = ie_navy) +
  geom_hline(yintercept = 0, linetype = 2) + labs(x = "Fitted", y = "Residual", title = "sqrt(Species)")
gridExtra::grid.arrange(b1, b2, nrow = 1)
c(lambda_hat = lambda_hat)
rbind(original = unlist(ncvTest(g)[c("ChiSquare", "p")]), sqrt = unlist(ncvTest(gs)[c("ChiSquare", "p")]))

round(cor(savings[, c("pop15", "pop75", "dpi", "ddpi")]), 2)
v <- vif(fit); round(v, 2)
# VIF by hand for pop75: regress it on the other predictors
R2_75 <- summary(lm(pop75 ~ pop15 + dpi + ddpi, data = savings))$r.squared
c(R2_75 = R2_75, VIF_75 = 1 / (1 - R2_75))

car::confidenceEllipse(fit, which.coef = c("pop15", "pop75"), levels = 0.95, col = ie_navy,
                       xlab = "coefficient of pop15", ylab = "coefficient of pop75", fill = TRUE, fill.alpha = 0.1)
abline(v = confint(fit)["pop15", ], h = confint(fit)["pop75", ], lty = 2, col = "grey40")
points(0, 0, pch = 3, col = ie_coral, lwd = 2)

set.seed(1)
sd_beta1 <- sapply(c(0, 0.9, 0.99), function(rho) {
  sd(replicate(3000, {
    x1 <- rnorm(50); x2 <- rho * x1 + sqrt(1 - rho^2) * rnorm(50)
    coef(lm(1 + x1 + x2 + rnorm(50) ~ x1 + x2))["x1"]
  }))
})
rbind(rho = c(0, 0.9, 0.99), sd = sd_beta1, ratio_to_rho0 = sd_beta1 / sd_beta1[1],
      theory = 1 / sqrt(1 - c(0, 0.9, 0.99)^2))

x0 <- data.frame(pop15 = 35, pop75 = 2, dpi = 1000, ddpi = 3)
rbind(CI = predict(fit, x0, interval = "confidence"), PI = predict(fit, x0, interval = "prediction"))
X <- model.matrix(fit); x0v <- c(1, 35, 2, 1000, 3)
q0 <- drop(t(x0v) %*% solve(crossprod(X)) %*% x0v)
tq <- qt(0.975, n - p)
c(q0 = q0, CI_half = tq * sigma(fit) * sqrt(q0), PI_half = tq * sigma(fit) * sqrt(1 + q0))

f1 <- lm(sr ~ pop15, data = savings)
grid <- data.frame(pop15 = seq(20, 50, length.out = 100))
ci <- predict(f1, grid, interval = "confidence"); pi <- predict(f1, grid, interval = "prediction")
b <- cbind(grid, fit = ci[, 1], cl = ci[, 2], cu = ci[, 3], pl = pi[, 2], pu = pi[, 3])
ggplot(b, aes(pop15)) +
  geom_ribbon(aes(ymin = pl, ymax = pu, fill = "95% prediction interval"), alpha = 0.35) +
  geom_ribbon(aes(ymin = cl, ymax = cu, fill = "95% CI for the mean"), alpha = 0.6) +
  geom_line(aes(y = fit), colour = ie_navy) + geom_point(data = savings, aes(pop15, sr), colour = ie_navy, size = 1.5) +
  scale_fill_manual(values = c("95% CI for the mean" = ie_blue, "95% prediction interval" = "#D2E6F4"), name = NULL) +
  labs(x = "pop15 (% population under 15)", y = "Savings rate sr (%)")

data(teengamb, package = "faraway")
tg <- lm(gamble ~ sex + status + income + verbal, data = teengamb)
ts <- lm(sqrt(gamble) ~ sex + status + income + verbal, data = teengamb)
x0 <- data.frame(sex = 0, status = 40, income = 5, verbal = 6)
pr_raw <- predict(tg, x0, interval = "prediction")
pr_sqrt <- predict(ts, x0, interval = "prediction")
rbind(raw_scale = pr_raw, sqrt_model_back_transformed = pmax(pr_sqrt, 0)^2)
c(ncv_p_raw = ncvTest(tg)$p, ncv_p_sqrt = ncvTest(ts)$p)

loo_mse <- mean((resid(fit) / (1 - hatvalues(fit)))^2)
c(RSS_over_n = mean(resid(fit)^2), sigma2_hat = sigma(fit)^2, LOOCV_MSE = loo_mse)

f4 <- lm(y ~ x, data.frame(x = 1:4, y = c(2, 3, 5, 6)))
round(cbind(h = hatvalues(f4), r = rstandard(f4), t = rstudent(f4), D = cooks.distance(f4)), 3)

set.seed(63)
x <- runif(100, 1, 10)
dA <- data.frame(x, y = 1 + x + rnorm(100))
dB <- data.frame(x, y = 1 + x + 0.3 * x^2 + rnorm(100))
dC <- data.frame(x, y = 1 + x + rnorm(100, sd = 0.2 * x))
plots <- lapply(list(A = dA, B = dB, C = dC), function(d) {
  f <- lm(y ~ x, d)
  ggplot(data.frame(f = fitted(f), r = resid(f)), aes(f, r)) + geom_point(colour = ie_navy, size = 1) +
    geom_hline(yintercept = 0, linetype = 2) + labs(x = "Fitted", y = "Residual")
})
gridExtra::grid.arrange(plots$A + ggtitle("A"), plots$B + ggtitle("B"), plots$C + ggtitle("C"), nrow = 1)

sl <- lm(stack.loss ~ ., data = stackloss)
round(summary(sl)$coefficients, 4); c(sigma = sigma(sl), R2 = summary(sl)$r.squared)
nn <- nrow(stackloss); pp <- length(coef(sl))
infl <- data.frame(h = hatvalues(sl), t = rstudent(sl), D = cooks.distance(sl))
round(infl[order(-infl$D)[1:4], ], 3); c(two_p_over_n = 2 * pp / nn, bonferroni = qt(1 - 0.05 / (2 * nn), nn - pp - 1))
round(vif(sl), 2)
sl2 <- update(sl, subset = -21); round(cbind(all = coef(sl), without_21 = coef(sl2)), 3)
new <- data.frame(Air.Flow = 58, Water.Temp = 20, Acid.Conc. = 85)
rbind(CI = predict(sl, new, interval = "confidence"), PI = predict(sl, new, interval = "prediction"))

gl <- lm(log(Species) ~ Area + Elevation + Nearest + Scruz + Adjacent, data = gala)
round(summary(gl)$coefficients, 4)
c(ncv_p_raw = ncvTest(g)$p, ncv_p_sqrt = ncvTest(gs)$p, ncv_p_log = ncvTest(gl)$p)
c(pct_per_100m = 100 * (exp(100 * coef(gl)["Elevation"]) - 1))

loo <- function(f) mean((resid(f) / (1 - hatvalues(f)))^2)
m_small <- lm(sr ~ pop15 + ddpi, data = savings)
c(full = loo(fit), small = loo(m_small), RSS_n_full = mean(resid(fit)^2), RSS_n_small = mean(resid(m_small)^2))
sapply(c("Libya", "Japan", "Zambia"), function(cn) {
  i <- which(rownames(savings) == cn)
  brute <- savings$sr[i] - predict(update(fit, subset = -i), savings[i, ])
  c(brute = unname(brute), formula = unname(resid(fit)[i] / (1 - hatvalues(fit)[i])))
})
