source("_common.R")
suppressPackageStartupMessages({library(car); library(glmnet); library(ggrepel)})
data(savings, package = "faraway")
data(prostate, package = "faraway")

fit_full <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = savings)
ry <- resid(lm(sr ~ pop75 + dpi + ddpi, data = savings))     # y adjusted for the others
rx <- resid(lm(pop15 ~ pop75 + dpi + ddpi, data = savings))  # pop15 adjusted for the others
c(multiple_regression = coef(fit_full)[["pop15"]], residual_regression = coef(lm(ry ~ rx))[["rx"]])

short <- lm(sr ~ pop75, data = savings)
long  <- lm(sr ~ pop75 + pop15, data = savings)
knitr::kable(round(rbind(`sr ~ pop75` = summary(short)$coefficients["pop75", ],
                         `sr ~ pop75 + pop15` = summary(long)$coefficients["pop75", ]), 3),
             caption = "The coefficient of pop75 with and without pop15 in the model.")

aux <- lm(pop15 ~ pop75, data = savings)
c(marginal = coef(short)[["pop75"]],
  adjusted_plus_bias = coef(long)[["pop75"]] + coef(long)[["pop15"]] * coef(aux)[["pop75"]],
  adjusted = coef(long)[["pop75"]], delta = coef(aux)[["pop75"]], beta15 = coef(long)[["pop15"]])

set.seed(707)
sim <- replicate(2000, {
  x1 <- rnorm(100); x2 <- 0.6 * x1 + 0.8 * rnorm(100)    # Delta = 0.6
  y <- 1 + x1 + 2 * x2 + rnorm(100)
  c(short = coef(lm(y ~ x1))[["x1"]], long = coef(lm(y ~ x1 + x2))[["x1"]])
})
rbind(mean = rowMeans(sim), sd = apply(sim, 1, sd))      # theory: 1 + 0.6 * 2 = 2.2 vs 1

pg_fit <- lm(weight ~ group, data = PlantGrowth)
round(summary(pg_fit)$coefficients, 3)
tapply(PlantGrowth$weight, PlantGrowth$group, mean)
anova(pg_fit)

pg2 <- PlantGrowth; contrasts(pg2$group) <- contr.sum(3)
round(summary(lm(weight ~ group, data = pg2))$coefficients, 3)                    # sum coding
pg3 <- transform(PlantGrowth, group = relevel(group, ref = "trt1"))
round(coef(lm(weight ~ group, data = pg3)), 3)                                    # baseline trt1
c(deviance_treatment = deviance(pg_fit), deviance_sum = deviance(lm(weight ~ group, data = pg2)))

mt <- mtcars; mt$cyl <- factor(mt$cyl)
fit_cyl <- lm(mpg ~ wt + cyl, data = mt)
round(summary(fit_cyl)$coefficients, 3)
anova(lm(mpg ~ wt, data = mt), fit_cyl)                      # 2-df test for the factor

mt$am <- factor(mt$am, labels = c("automatic", "manual"))
fit_am <- lm(mpg ~ wt * am, data = mt)
round(summary(fit_am)$coefficients, 3)
b <- coef(fit_am)
c(slope_auto = b[["wt"]], slope_manual = b[["wt"]] + b[["wt:ammanual"]],
  crossing_weight = -b[["ammanual"]] / b[["wt:ammanual"]],
  gap_at_2 = b[["ammanual"]] + b[["wt:ammanual"]] * 2, gap_at_3 = b[["ammanual"]] + b[["wt:ammanual"]] * 3,
  gap_at_4 = b[["ammanual"]] + b[["wt:ammanual"]] * 4)

grid <- expand.grid(wt = seq(1.5, 5.5, length.out = 60), am = levels(mt$am))
grid$mpg <- predict(fit_am, grid)
ggplot(mt, aes(wt, mpg, colour = am)) + geom_point() + geom_line(data = grid) +
  scale_colour_manual(values = c(automatic = ie_navy, manual = ie_coral), name = NULL) +
  labs(x = "Weight (1000 lb)", y = "Miles per gallon")
range(mtcars$wt[mtcars$am == 0]); range(mtcars$wt[mtcars$am == 1])

fit_wh <- lm(mpg ~ wt * hp, data = mtcars)
round(summary(fit_wh)$coefficients, 4)
b <- coef(fit_wh)
round(setNames(b[["wt"]] + b[["wt:hp"]] * c(100, 150, 250), paste0("hp=", c(100, 150, 250))), 2)

rss <- function(f, d = mtcars) round(deviance(lm(f, d)), 1)
data.frame(
  model = c("I(hp^2)", "hp + I(hp^2)", "wt + wt:hp (no hp)", "wt * hp"),
  RSS = c(rss(mpg ~ I(hp^2)), rss(mpg ~ hp + I(hp^2)), rss(mpg ~ wt + wt:hp), rss(mpg ~ wt * hp)),
  RSS_shifted = c(rss(mpg ~ I((hp - 150)^2)), rss(mpg ~ I(hp - 150) + I((hp - 150)^2)),
                  rss(mpg ~ I(wt - 3) + I(wt - 3):hp), rss(mpg ~ I(wt - 3) * I(hp - 150)))
)

m <- lm(mpg ~ wt + hp + am, data = mt)
anova(m)                                   # sequential
anova(lm(mpg ~ hp + wt + am, data = mt))   # different order, different rows (first two)
Anova(m)                                   # each term after the others (Type II)

full <- lm(lpsa ~ ., data = prostate)
n <- nrow(prostate); p_full <- length(coef(full))
RSS <- deviance(full)
c(AIC = AIC(full), AIC_by_hand = n * log(2 * pi) + n * log(RSS / n) + n + 2 * (p_full + 1),
  extractAIC = extractAIC(full)[2], difference = AIC(full) - extractAIC(full)[2],
  const = n * log(2 * pi) + n + 2)

m_aic <- step(full, trace = 0)
m_bic <- step(full, trace = 0, k = log(n))
sel <- function(m) paste(setdiff(names(coef(m)), "(Intercept)"), collapse = " + ")
data.frame(model = c("full", "step AIC", "step BIC"),
           predictors = c(sel(full), sel(m_aic), sel(m_bic)),
           k = c(10, 7, 5),
           AIC = round(c(AIC(full), AIC(m_aic), AIC(m_bic)), 2),
           BIC = round(c(BIC(full), BIC(m_aic), BIC(m_bic)), 2))

library(leaps)
rs <- summary(regsubsets(lpsa ~ ., data = prostate, nvmax = 8))
data.frame(p = 2:9, RSS = round(rs$rss, 2), Cp = round(rs$cp, 2), adjR2 = round(rs$adjr2, 3), BIC = round(rs$bic, 1))

set.seed(2027)
res <- replicate(100, {
  d <- as.data.frame(matrix(rnorm(100 * 21), 100)); names(d)[1] <- "y"
  m <- step(lm(y ~ ., data = d), trace = 0)
  k <- length(coef(m)) - 1
  sig <- if (k > 0) sum(summary(m)$coefficients[-1, 4] < 0.05) else 0
  c(kept = k, significant = sig)
})
c(mean_kept = mean(res["kept", ]), mean_significant = mean(res["significant", ]),
  share_significant = sum(res["significant", ]) / sum(res["kept", ]))

x <- as.matrix(prostate[, names(prostate) != "lpsa"]); y <- prostate$lpsa
set.seed(606)
foldid <- sample(rep(1:10, length.out = nrow(x)))
cv_lasso <- cv.glmnet(x, y, alpha = 1, foldid = foldid)
cv_ridge <- cv.glmnet(x, y, alpha = 0, foldid = foldid)
cv_enet  <- cv.glmnet(x, y, alpha = 0.5, foldid = foldid)
rbind(lasso = c(lambda.min = cv_lasso$lambda.min, lambda.1se = cv_lasso$lambda.1se),
      ridge = c(cv_ridge$lambda.min, cv_ridge$lambda.1se),
      enet  = c(cv_enet$lambda.min, cv_enet$lambda.1se))
round(cbind(OLS = coef(full), ridge = coef(cv_ridge, s = "lambda.1se")[, 1],
            lasso = coef(cv_lasso, s = "lambda.1se")[, 1], enet = coef(cv_enet, s = "lambda.1se")[, 1]), 3)

op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
plot(cv_lasso$glmnet.fit, xvar = "lambda", label = TRUE)
plot(cv_lasso); par(op)

set.seed(1)
res <- t(sapply(1:50, function(r) {
  id <- sample(nrow(prostate), 30); tr <- prostate[-id, ]; te <- prostate[id, ]
  xtr <- as.matrix(tr[, names(tr) != "lpsa"]); xte <- as.matrix(te[, names(te) != "lpsa"])
  ft <- lm(lpsa ~ ., data = tr)
  fa <- step(ft, trace = 0); fb <- step(ft, trace = 0, k = log(nrow(tr)))
  cvr <- cv.glmnet(xtr, tr$lpsa, alpha = 1, foldid = sample(rep(1:10, length.out = nrow(tr))))
  mse <- function(pred) mean((te$lpsa - pred)^2)
  c(OLS = mse(predict(ft, te)), step_AIC = mse(predict(fa, te)), step_BIC = mse(predict(fb, te)),
    lasso_min = mse(predict(cvr, xte, s = "lambda.min")), lasso_1se = mse(predict(cvr, xte, s = "lambda.1se")))
}))
rbind(mean = colMeans(res), sd = apply(res, 2, sd))

rx <- resid(lm(wt ~ hp, data = mtcars))
c(fwl = sum(rx * mtcars$mpg) / sum(rx^2), lm = coef(lm(mpg ~ wt + hp, data = mtcars))[["wt"]])

m1 <- lm(mpg ~ wt, mtcars); m2 <- lm(mpg ~ wt + hp, mtcars); a <- lm(hp ~ wt, mtcars)
c(delta = coef(a)[["wt"]], marginal = coef(m1)[["wt"]], reconstructed = coef(m2)[["wt"]] + coef(m2)[["hp"]] * coef(a)[["wt"]])

d <- data.frame(x = 0:3, y = c(1.3, 0, 1.1, 4.2))
c(RSS_x2 = deviance(lm(y ~ I(x^2), d)), RSS_shift = deviance(lm(y ~ I((x - 1)^2), d)),
  RSS_full = deviance(lm(y ~ x + I(x^2), d)), RSS_full_shift = deviance(lm(y ~ I(x - 1) + I((x - 1)^2), d)))

n <- 40; p <- c(3, 5, 9); rss <- c(52, 45, 41.5)
aic <- n * log(rss / n) + 2 * p; bic <- n * log(rss / n) + p * log(n)
aicc <- aic + 2 * (p + 1) * (p + 2) / (n - p - 2)
round(rbind(AIC = aic, BIC = bic, AICc_with_k = aicc), 2)

data(seatpos, package = "faraway")
sf <- lm(hipcenter ~ ., data = seatpos)
c(R2 = summary(sf)$r.squared, F = summary(sf)$fstatistic[[1]],
  p_F = pf(summary(sf)$fstatistic[[1]], 8, 29, lower.tail = FALSE))
round(summary(sf)$coefficients, 3)
round(vif(sf), 1); cor(seatpos$Ht, seatpos$HtShoes)
xs <- as.matrix(seatpos[, names(seatpos) != "hipcenter"]); ys <- seatpos$hipcenter
set.seed(7); fid <- sample(rep(1:10, length.out = nrow(xs)))
cv_err <- function(form) mean(sapply(1:10, function(k) {
  f <- lm(form, seatpos[fid != k, ]); mean((ys[fid == k] - predict(f, seatpos[fid == k, ]))^2) }))
cr <- cv.glmnet(xs, ys, alpha = 0, foldid = fid); cl <- cv.glmnet(xs, ys, alpha = 1, foldid = fid)
c(OLS_full = cv_err(hipcenter ~ .), OLS_small = cv_err(hipcenter ~ Age + Weight + Ht),
  ridge = min(cr$cvm), lasso_min = min(cl$cvm))
round(cbind(ridge = coef(cr, s = "lambda.min")[, 1], lasso = coef(cl, s = "lambda.1se")[, 1]), 3)
