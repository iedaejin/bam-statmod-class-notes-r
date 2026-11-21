source("_common.R")

library(mgcv)
data(airquality)
aq <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])
add <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "REML")

set.seed(1)
k.check(add)

ggplot(aq, aes(Temp, Wind, colour = log(Ozone))) + geom_point(size = 2.2) +
  scale_colour_gradient(low = ie_blue, high = ie_coral) +
  labs(x = "Temperature (F)", y = "Wind (mph)", colour = "log(Ozone)") +
  theme(legend.position = "right")

set.seed(1); d <- data.frame(x = runif(30), z = runif(30))
sx <- smoothCon(s(x, bs = "cr", k = 5), data = d, absorb.cons = FALSE)[[1]]
sz <- smoothCon(s(z, bs = "cr", k = 4), data = d, absorb.cons = FALSE)[[1]]
tp <- smoothCon(te(x, z, bs = "cr", k = c(5, 4)), data = d, absorb.cons = FALSE)[[1]]
byhand <- t(sapply(1:30, function(i) kronecker(sx$X[i, ], sz$X[i, ])))
c(dim(tp$X), max_difference = max(abs(byhand - tp$X)))

library(splines)
kn <- c(rep(0, 4), 0.5, rep(1, 4)); gx <- seq(0, 1, length.out = 100)
Bg <- splineDesign(kn, gx, ord = 4, outer.ok = TRUE)
prs <- list(c(2, 2), c(3, 3), c(4, 2))
dd <- do.call(rbind, lapply(prs, function(p) {
  g <- expand.grid(i = seq_along(gx), j = seq_along(gx))
  data.frame(x = gx[g$i], z = gx[g$j], v = Bg[g$i, p[1]] * Bg[g$j, p[2]],
             lab = sprintf("a%d(x) c%d(z)", p[1], p[2]))
}))
ggplot(dd, aes(x, z, fill = v)) + geom_raster() + facet_wrap(~lab) + coord_equal() +
  scale_fill_gradient(low = "white", high = ie_navy) + theme(legend.position = "right") +
  labs(fill = NULL)

Sx <- kronecker(sx$S[[1]], diag(4)); Sz <- kronecker(diag(5), sz$S[[1]])
ratio <- function(A, B) { i <- abs(B) > 1e-6; range(A[i] / B[i]) }
rbind(x_penalty = ratio(tp$S[[1]], Sx), z_penalty = ratio(tp$S[[2]], Sz))

c(rank_x = qr(tp$S[[1]])$rank, rank_z = qr(tp$S[[2]])$rank,
  rank_sum = qr(tp$S[[1]] + tp$S[[2]])$rank, KL = ncol(tp$X))

aq$TempC <- (aq$Temp - 32) * 5 / 9; aq$Windk <- aq$Wind * 1.609
fits <- list(
  iso_F = gam(log(Ozone) ~ s(Solar.R) + s(Temp, Wind),   data = aq, method = "REML"),
  iso_C = gam(log(Ozone) ~ s(Solar.R) + s(TempC, Windk), data = aq, method = "REML"),
  te_F  = gam(log(Ozone) ~ s(Solar.R) + te(Temp, Wind),   data = aq, method = "REML"),
  te_C  = gam(log(Ozone) ~ s(Solar.R) + te(TempC, Windk), data = aq, method = "REML"))
sapply(fits, AIC)

te1 <- gam(log(Ozone) ~ s(Solar.R) + te(Temp, Wind), data = aq, method = "REML")
summary(te1)
te1$sp

grid <- expand.grid(Temp = seq(56, 97, length.out = 70), Wind = seq(2, 21, length.out = 70))
far <- exclude.too.far(grid$Temp, grid$Wind, aq$Temp, aq$Wind, dist = 0.1)
grid$z <- predict(te1, transform(grid, Solar.R = median(aq$Solar.R)), type = "terms")[, "te(Temp,Wind)"]
grid$z[far] <- NA
ggplot(grid, aes(Temp, Wind)) + geom_raster(aes(fill = z)) +
  geom_contour(aes(z = z), colour = "white", linewidth = .3, bins = 10, na.rm = TRUE) +
  geom_point(data = aq, aes(Temp, Wind), size = .7, colour = "grey15") +
  scale_fill_gradient2(low = ie_blue, mid = "white", high = ie_coral, midpoint = 0, na.value = "transparent") +
  labs(x = "Temperature (F)", y = "Wind (mph)", fill = "effect") + theme(legend.position = "right")

ti_toy <- smoothCon(ti(x, z, k = c(5, 4)), data = d, absorb.cons = TRUE)[[1]]
te_toy <- smoothCon(te(x, z, k = c(5, 4)), data = d, absorb.cons = TRUE)[[1]]
c(ti = ncol(ti_toy$X), `(K-1)(L-1)` = 4 * 3, te = ncol(te_toy$X), `KL-1` = 5 * 4 - 1)

bad <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp) + te(Temp, Wind), data = aq, method = "REML")
round(concurvity(bad, full = TRUE)["worst", ], 3)

ti1 <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp) + ti(Temp, Wind), data = aq, method = "REML")
summary(ti1)$s.table

round(concurvity(ti1, full = TRUE)[c("worst", "estimate"), ], 2)

te1 <- gam(log(Ozone) ~ s(Solar.R) + te(Temp, Wind), data = aq, method = "REML")
tab <- data.frame(model = c("additive", "s(Solar.R) + te(Temp, Wind)", "additive + ti(Temp, Wind)"),
                  logLik = sapply(list(add, te1, ti1), function(m) as.numeric(logLik(m))),
                  df = sapply(list(add, te1, ti1), function(m) attr(logLik(m), "df")),
                  AIC = sapply(list(add, te1, ti1), AIC))
tab$delta <- tab$AIC - min(tab$AIC)
tab$weight <- exp(-tab$delta / 2) / sum(exp(-tab$delta / 2))
knitr::kable(tab, digits = c(0, 2, 2, 2, 2, 3), booktabs = TRUE,
             caption = "AIC, AIC differences and Akaike weights for three models of log(Ozone).")

add_ml <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "ML")
ti_ml  <- update(add_ml, . ~ . + ti(Temp, Wind))
anova(add_ml, ti_ml, test = "F")
anova(add, ti1, test = "F")

set.seed(18); n <- 100
pv <- replicate(200, {
  x1 <- runif(n); x2 <- runif(n); y <- sin(2 * pi * x1) + 4 * (x2 - .5)^2 + rnorm(n, sd = .5)
  dd <- data.frame(y, x1, x2)
  a <- gam(y ~ s(x1) + s(x2), data = dd, method = "ML")
  b <- gam(y ~ s(x1) + s(x2) + ti(x1, x2), data = dd, method = "ML")
  c(anova = anova(a, b, test = "F")[2, "Pr(>F)"], wald = summary(b)$s.table["ti(x1,x2)", "p-value"])
})
rowMeans(pv < 0.05, na.rm = TRUE)

set.seed(2027)
aq$noise1 <- rnorm(nrow(aq)); aq$noise2 <- runif(nrow(aq))
f_big <- log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp) + s(noise1) + s(noise2)
m_no  <- gam(f_big, data = aq, method = "REML")
m_sel <- gam(f_big, data = aq, method = "REML", select = TRUE)
round(cbind(without = summary(m_no)$s.table[, "edf"], select = summary(m_sel)$s.table[, "edf"]), 4)

X <- model.matrix(add); jj <- add$smooth[[2]]; idx <- jj$first.para:jj$last.para
Xj <- X[, idx]; Xo <- X[, -idx]
Qo <- qr.Q(qr(Xo)); Qj <- qr.Q(qr(Xj)); PXj <- Qo %*% crossprod(Qo, Xj)
b <- coef(add)[idx]
hand <- c(worst = max(svd(crossprod(Qo, Qj))$d)^2,
          observed = sum((PXj %*% b)^2) / sum((Xj %*% b)^2),
          estimate = sum(PXj^2) / sum(Xj^2))
rbind(by_hand = hand, concurvity = concurvity(add)[, "s(Wind)"])

f1 <- function(x) sin(2 * pi * x); f2 <- function(x) 2 * (x - .5)^2
conc_sim <- function(s, seed = 25) {
  set.seed(seed); n <- 300; x1 <- runif(n)
  x2 <- if (is.na(s)) runif(n) else pmin(pmax(x1 + rnorm(n, sd = s), 0), 1)
  y <- f1(x1) + f2(x2) + rnorm(n, sd = .3)
  m <- gam(y ~ s(x1) + s(x2), method = "REML")
  se <- predict(m, type = "terms", se.fit = TRUE)$se.fit[, "s(x1)"]
  c(cor = cor(x1, x2), concurvity = concurvity(m, full = TRUE)["estimate", "s(x1)"], mean_se = mean(se))
}
res <- t(sapply(c(NA, .3, .08, .02), conc_sim))
rownames(res) <- c("independent", "s = 0.30", "s = 0.08", "s = 0.02")
round(res, 3)

fit_band <- function(s, seed = 25) {
  set.seed(seed); n <- 300; x1 <- runif(n)
  x2 <- if (is.na(s)) runif(n) else pmin(pmax(x1 + rnorm(n, sd = s), 0), 1)
  y <- f1(x1) + f2(x2) + rnorm(n, sd = .3)
  m <- gam(y ~ s(x1) + s(x2), method = "REML")
  xg <- seq(0, 1, length.out = 150)
  p <- predict(m, data.frame(x1 = xg, x2 = .5), type = "terms", se.fit = TRUE)
  data.frame(x = xg, fit = p$fit[, "s(x1)"], se = p$se.fit[, "s(x1)"], truth = f1(xg) - mean(f1(x1)),
             design = ifelse(is.na(s), "independent", "x2 close to x1 (s = 0.08)"))
}
db <- rbind(fit_band(NA), fit_band(.08))
ggplot(db, aes(x, fit)) + geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), fill = "#D2E6F4") +
  geom_line(aes(y = truth), colour = ie_coral, linetype = 2) + geom_line(colour = ie_navy, linewidth = 1) +
  facet_wrap(~design) + labs(x = expression(x[1]), y = expression(hat(f)[1]))

set.seed(1)
k.check(ti1)

aqd <- airquality; aqd$day <- seq_len(nrow(aqd)); aqd <- na.omit(aqd)
r <- residuals(ti1, type = "deviance")
oka <- diff(aqd$day) == 1
lag1 <- cor(r[-length(r)][oka], r[-1][oka])
mon <- anova(gam(r ~ factor(aqd$Month)))
dd <- rbind(data.frame(x = aqd$day, r = r, what = "day of the season"),
            data.frame(x = fitted(ti1), r = r, what = "fitted log(Ozone)"))
ggplot(dd, aes(x, r)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(colour = ie_navy, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, colour = ie_coral, fill = "grey80") +
  facet_wrap(~what, scales = "free_x") + labs(x = NULL, y = "residual")

nd <- expand.grid(Temp = c(70, 80, 90), Wind = c(5, 10, 15), Solar.R = median(aq$Solar.R))
pr <- predict(ti1, nd, se.fit = TRUE)
nd$median <- exp(pr$fit); nd$lower <- exp(pr$fit - 1.96 * pr$se.fit); nd$upper <- exp(pr$fit + 1.96 * pr$se.fit)
knitr::kable(round(nd, 1), booktabs = TRUE,
             caption = "Predicted median ozone (ppb) with 95% intervals from the interaction model (Solar.R at its median, 207).")

contrast <- function(m, t1, t2, w, S = 207) {
  X <- predict(m, data.frame(Temp = c(t1, t2), Wind = w, Solar.R = S), type = "lpmatrix")
  d <- X[2, ] - X[1, ]; est <- sum(d * coef(m)); se <- sqrt(drop(t(d) %*% vcov(m) %*% d))
  c(log_ratio = est, se = se, ratio = exp(est), lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se))
}
round(rbind(`wind 5` = contrast(ti1, 70, 90, 5), `wind 15` = contrast(ti1, 70, 90, 15)), 3)

gratia::draw(add, residuals = TRUE, nrow = 1)

set.seed(3); dd <- data.frame(x = runif(200), y = runif(200), z = runif(200))
dd$r <- rnorm(200)
cnt <- function(f) { g <- gam(f, data = dd, fit = FALSE); sm <- g$smooth[[1]]
  c(coefficients = sm$last.para - sm$first.para + 1, smoothing_par = length(sm$S)) }
rbind(te_xz = cnt(r ~ te(x, z, k = c(6, 5))), ti_xz = cnt(r ~ ti(x, z, k = c(6, 5))),
      te_xyz = cnt(r ~ te(x, y, z, k = 4)))

data(trees)
gA <- gam(Volume ~ s(Girth) + s(Height), family = Gamma(link = "log"), data = trees, method = "REML")
gB <- update(gA, . ~ . + ti(Girth, Height))
gL <- gam(Volume ~ log(Girth) + log(Height), family = Gamma(link = "log"), data = trees, method = "REML")
AIC(gA, gB, gL)
round(summary(gL)$p.table, 3)
round(summary(gB)$s.table, 4)
round(concurvity(gA, full = TRUE), 2)

t0 <- 1:4; Q <- qr.Q(qr(cbind(1, t0)))
conc1 <- function(q) sum((Q %*% crossprod(Q, q))^2) / sum(q^2)
c(q = conc1(c(-2, -1, 2, 1)), q_prime = conc1(c(1, -1, -1, 1)))

res2 <- t(sapply(c(.5, .15, .04), conc_sim, seed = 26))
rownames(res2) <- paste("s =", c(.5, .15, .04)); round(res2, 3)

set.seed(11)
out <- t(replicate(60, {
  aq$noise3 <- rnorm(nrow(aq))
  m3 <- gam(log(Ozone) ~ s(Solar.R) + s(Wind) + s(Temp) + s(noise3), data = aq,
            method = "REML", select = TRUE)
  summary(m3)$s.table["s(noise3)", c("edf", "p-value")]
}))
c(median_edf = median(out[, 1]), share_edf_below_half = mean(out[, 1] < .5),
  share_p_below_5pct = mean(out[, 2] < .05))

cc <- function(m, unc) {
  X <- predict(m, data.frame(Temp = c(80, 90), Wind = 10, Solar.R = 207), type = "lpmatrix")
  d <- X[2, ] - X[1, ]; est <- sum(d * coef(m)); se <- sqrt(drop(t(d) %*% vcov(m, unconditional = unc) %*% d))
  c(ratio = exp(est), lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se), se = se)
}
rbind(conditional = cc(ti1, FALSE), unconditional = cc(ti1, TRUE))
X5 <- predict(ti1, data.frame(Temp = c(70, 90), Wind = 5, Solar.R = 207), type = "lpmatrix")
X15 <- predict(ti1, data.frame(Temp = c(70, 90), Wind = 15, Solar.R = 207), type = "lpmatrix")
d <- (X5[2, ] - X5[1, ]) - (X15[2, ] - X15[1, ])
est <- sum(d * coef(ti1)); se <- sqrt(drop(t(d) %*% vcov(ti1) %*% d))
c(est = est, se = se, z = est / se, ratio_of_ratios = exp(est), lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se))
