source("_common.R")

library(MASS)
data(mcycle)
x <- mcycle$times
y <- mcycle$accel
n <- length(y)

ggplot(mcycle, aes(times, accel)) +
  geom_point(colour = ie_navy, size = 1.8) +
  labs(x = "Time after impact (ms)", y = "Head acceleration (g)")

degrees <- c(1, 3, 5, 8, 12, 20)
fits <- lapply(degrees, function(q) lm(accel ~ poly(times, q), data = mcycle))
tab <- data.frame(degree = degrees,
                  RSS = sapply(fits, function(f) round(sum(resid(f)^2))),
                  R2 = sapply(fits, function(f) round(summary(f)$r.squared, 3)),
                  AIC = sapply(fits, function(f) round(AIC(f), 1)))
tab

grid <- data.frame(times = seq(min(x), max(x), length.out = 300))
pf <- do.call(rbind, lapply(c(3, 8, 12, 20), function(q) {
  m <- lm(accel ~ poly(times, q), data = mcycle)
  data.frame(times = grid$times, fit = predict(m, grid), degree = paste("degree", q))
}))
pf$degree <- factor(pf$degree, levels = paste("degree", c(3, 8, 12, 20)))
ggplot(pf, aes(times, fit)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey60", size = 1) +
  geom_line(colour = ie_coral, linewidth = 0.9) + facet_wrap(~degree, nrow = 1) +
  coord_cartesian(ylim = c(-150, 100)) + labs(x = "Time (ms)", y = "Acceleration (g)")

lev <- do.call(rbind, lapply(c(5, 12), function(q) {
  m <- lm(accel ~ poly(times, q), data = mcycle)
  data.frame(times = x, h = hatvalues(m), degree = paste("degree", q))
}))
ggplot(lev, aes(times, h, colour = degree)) + geom_point(size = 1.4) +
  scale_colour_manual(values = c(ie_navy, ie_coral), name = NULL) +
  labs(x = "Time (ms)", y = expression("leverage " * h[ii]))
m12 <- lm(accel ~ poly(times, 12), data = mcycle)
c(last = unname(hatvalues(m12)[n]), median = median(hatvalues(m12)), average = 13 / n)

new <- data.frame(times = c(57.6, 60, 65))
sapply(c(3, 8, 12), function(q) round(predict(lm(accel ~ poly(times, q), data = mcycle), new), 1))

signif(sapply(c(3, 5, 8, 10, 12), function(q) kappa(outer(x, 0:q, `^`), exact = TRUE)), 3)
sum(is.na(coef(lm(accel ~ poly(times, 15, raw = TRUE), data = mcycle))))

Q <- poly(x, 3)
Xraw <- cbind(1, x, x^2, x^3)
A <- solve(crossprod(Xraw), crossprod(Xraw, cbind(1, Q)))   # [1, Q] = Xraw %*% A
fit_orth <- lm(y ~ poly(x, 3))
fit_raw  <- lm(y ~ poly(x, 3, raw = TRUE))
cbind(converted = drop(A %*% coef(fit_orth)), raw_fit = coef(fit_raw))
c(max_error_A = max(abs(Xraw %*% A - cbind(1, Q))),
  max_diff_fitted = max(abs(fitted(fit_orth) - fitted(fit_raw))))

set.seed(2027)
folds <- sample(rep(1:5, length.out = n))
cv_degree <- function(q) {
  sse <- 0
  for (k in 1:5) {
    tr <- folds != k
    m <- lm(accel ~ poly(times, q), data = mcycle[tr, ])
    sse <- sse + sum((mcycle$accel[!tr] - predict(m, mcycle[!tr, ]))^2)
  }
  sse / n
}
cv <- sapply(1:15, cv_degree)
round(cv, 0)
which.min(cv)

nw_S <- function(x0, x, h) {                       # Nadaraya-Watson weights, one row per x0
  K <- dnorm(outer(x0, x, "-") / h)
  K / rowSums(K)
}
i0 <- which.min(abs(x - 20))
w_bin <- as.numeric(cut(x, seq(0, 60, 5)) == cut(x[i0], seq(0, 60, 5))); w_bin <- w_bin / sum(w_bin)
w_knn <- numeric(n); w_knn[order(abs(x - x[i0]))[1:15]] <- 1 / 15
w_nw <- nw_S(x[i0], x, 2)[1, ]
d <- abs(x - x[i0]); Delta <- sort(d)[ceiling(0.3 * n)]
w_tri <- ifelse(d < Delta, (1 - (d / Delta)^3)^3, 0); w_tri <- w_tri / sum(w_tri)
wd <- data.frame(times = rep(x, 4), w = c(w_bin, w_knn, w_nw, w_tri),
                 scheme = factor(rep(c("bins (5 ms)", "kNN (k = 15)", "Gaussian (h = 2)", "tricube (span 0.3)"),
                                     each = n), levels = c("bins (5 ms)", "kNN (k = 15)", "Gaussian (h = 2)", "tricube (span 0.3)")))
ggplot(wd, aes(times, w)) + geom_segment(aes(xend = times, yend = 0), colour = ie_navy) +
  geom_vline(xintercept = x[i0], colour = ie_coral, linetype = 2) + facet_wrap(~scheme, nrow = 1) +
  coord_cartesian(xlim = c(5, 35)) + labs(x = "Time (ms)", y = "weight")

kernels <- list(
  Gaussian     = function(u) dnorm(u),
  Epanechnikov = function(u) 0.75 * (1 - u^2) * (abs(u) < 1),
  biweight     = function(u) 15 / 16 * (1 - u^2)^2 * (abs(u) < 1),
  tricube      = function(u) 70 / 81 * (1 - abs(u)^3)^3 * (abs(u) < 1),
  uniform      = function(u) 0.5 * (abs(u) < 1))
kmom <- t(sapply(kernels, function(K) c(
  integral = integrate(K, -Inf, Inf)$value,
  mu2 = integrate(function(u) u^2 * K(u), -Inf, Inf)$value,
  R = integrate(function(u) K(u)^2, -Inf, Inf)$value)))
round(kmom, 4)

gridx <- seq(min(x), max(x), length.out = 300)
hs3 <- c(0.5, 2, 6)
kf <- do.call(rbind, lapply(hs3, function(h)
  data.frame(times = gridx, fit = drop(nw_S(gridx, x, h) %*% y), h = paste("h =", h))))
ggplot(kf, aes(times, fit, colour = h)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey65", size = 1, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = c(ie_coral, ie_navy, ie_green), name = NULL) +
  labs(x = "Time (ms)", y = "Acceleration (g)")

h <- 2
b <- h / (0.25 / qnorm(0.75))                      # bandwidth argument that gives sd h
ks <- ksmooth(x, y, kernel = "normal", bandwidth = b, x.points = x)
c(b = b, max_abs_diff = max(abs(ks$y - drop(nw_S(x, x, h) %*% y))))

ll_S <- function(x0, x, h) {                      # local linear smoother matrix (see "Local polynomial regression")
  vapply(seq_along(x0), function(i) {
    d <- x - x0[i]; w <- dnorm(d / h)
    s0 <- sum(w); s1 <- sum(w * d); s2 <- sum(w * d^2)
    w * (s2 - d * s1) / (s0 * s2 - s1^2)
  }, numeric(length(x)))|> t()
}
set.seed(2027)
xs <- seq(0, 1, length.out = 100); gx <- c(0, 0.05, 0.1, 0.2, 0.5)
SN <- nw_S(gx, xs, 0.1); SL <- ll_S(gx, xs, 0.1)
est <- replicate(300, { ys <- 1 + 3 * xs + rnorm(100, sd = 0.3); c(SN %*% ys, SL %*% ys) })
bias <- data.frame(x = rep(gx, 2), method = rep(c("NW", "local linear"), each = length(gx)),
                   bias = rowMeans(est) - rep(1 + 3 * gx, 2))
ggplot(bias, aes(x, bias, colour = method)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point() + scale_colour_manual(values = c(ie_coral, ie_navy), name = NULL) +
  labs(x = "x (the boundary is at 0)", y = "average estimate minus truth")

cmp <- data.frame(times = rep(gridx, 2),
                  fit = c(drop(nw_S(gridx, x, 3) %*% y), drop(ll_S(gridx, x, 3) %*% y)),
                  method = rep(c("Nadaraya-Watson", "local linear"), each = length(gridx)))
ggplot(cmp, aes(times, fit, colour = method)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey65", size = 1, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = c(ie_navy, ie_coral), name = NULL) +
  labs(x = "Time (ms)", y = "Acceleration (g)")

sp <- c(0.1, 0.3, 0.75)
lf <- do.call(rbind, lapply(sp, function(s) {
  m <- loess(accel ~ times, mcycle, span = s, degree = 2)
  data.frame(times = gridx, fit = predict(m, data.frame(times = gridx)),
             span = sprintf("span = %.2f (tr S = %.1f)", s, m$trace.hat))
}))
ggplot(lf, aes(times, fit, colour = span)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey65", size = 1, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = c(ie_coral, ie_navy, ie_green), name = NULL) +
  coord_cartesian(ylim = c(-150, 100)) + labs(x = "Time (ms)", y = "Acceleration (g)")
lo <- loess(accel ~ times, data = mcycle, span = 0.3, degree = 2)
c(trace_S = lo$trace.hat, enp = lo$enp, RSS = sum(resid(lo)^2))

edf_nw <- sapply(c(0.25, 0.5, 1, 2, 4, 8, 16, 100), function(h) sum(diag(nw_S(x, x, h))))
round(edf_nw, 2)

xt <- 1:5; yt <- c(2, 4, 3, 7, 9)
Kt <- pmax(1 - abs(outer(xt, xt, "-")) / 2, 0)
St <- Kt / rowSums(Kt)
fractions <- MASS::fractions
fractions(St)
ft <- drop(St %*% yt); rt <- yt - ft
rbind(fitted = ft, residual = rt, loo = rt / (1 - diag(St)))
c(edf = sum(diag(St)), RSS = sum(rt^2),
  LOOCV = mean((rt / (1 - diag(St)))^2), GCV = 5 * sum(rt^2) / (5 - sum(diag(St)))^2)

f_true <- drop(ll_S(x, x, 2.5) %*% y)
hg <- seq(0.5, 5, by = 0.05); sigma <- 23
bv <- t(sapply(hg, function(h) {
  S <- nw_S(x, x, h)
  c(bias2 = mean((S %*% f_true - f_true)^2), var = sigma^2 * mean(rowSums(S^2)))
}))
mse <- rowSums(bv)
bvl <- data.frame(h = rep(hg, 3), value = c(bv, mse),
                  part = rep(c("squared bias", "variance", "MSE"), each = length(hg)))
ggplot(bvl, aes(h, value, colour = part)) + geom_line(linewidth = 0.9) +
  geom_vline(xintercept = hg[which.min(mse)], linetype = 3) +
  scale_colour_manual(values = c(ie_coral, ie_navy, ie_blue), name = NULL) +
  labs(x = "bandwidth h (ms)", y = "average over design points")
hg[which.min(mse)]

h <- 2
S <- nw_S(x, x, h)
r <- y - drop(S %*% y)
shortcut <- mean((r / (1 - diag(S)))^2)
brute <- mean(sapply(1:n, function(i) {
  w <- dnorm((x[-i] - x[i]) / h)
  (y[i] - sum(w * y[-i]) / sum(w))^2
}))
c(shortcut = shortcut, brute_force = brute)

crit <- function(S) {
  r <- y - drop(S %*% y); edf <- sum(diag(S))
  c(edf = edf, LOOCV = mean((r / (1 - diag(S)))^2), GCV = n * sum(r^2) / (n - edf)^2)
}
hs <- seq(0.5, 5, by = 0.01)
tab_nw <- t(sapply(hs, function(h) crit(nw_S(x, x, h))))
tab_ll <- t(sapply(hs, function(h) crit(ll_S(x, x, h))))
best <- function(tab, col) c(h = hs[which.min(tab[, col])], tab[which.min(tab[, col]), ])
opt <- rbind(NW_LOOCV = best(tab_nw, "LOOCV"), NW_GCV = best(tab_nw, "GCV"),
             LL_LOOCV = best(tab_ll, "LOOCV"), LL_GCV = best(tab_ll, "GCV"))
cvd <- data.frame(h = rep(hs, 4),
                  value = c(tab_nw[, "LOOCV"], tab_nw[, "GCV"], tab_ll[, "LOOCV"], tab_ll[, "GCV"]),
                  method = rep(c("NW", "local linear"), each = 2 * length(hs)),
                  criterion = rep(rep(c("LOOCV", "GCV"), each = length(hs)), 2))
ggplot(cvd, aes(h, value, colour = method, linetype = criterion)) + geom_line(linewidth = 0.9) +
  geom_point(data = data.frame(h = opt[, "h"], value = opt[, c("LOOCV")] * c(1, NA, 1, NA) + opt[, "GCV"] * c(NA, 1, NA, 1),
                               method = rep(c("NW", "local linear"), each = 2),
                               criterion = rep(c("LOOCV", "GCV"), 2)), size = 2.5) +
  scale_colour_manual(values = c(ie_coral, ie_navy), name = NULL) + coord_cartesian(ylim = c(540, 1000)) +
  labs(x = "bandwidth h (ms)", y = "criterion")
round(opt, 2)

target <- 12
h_nw <- uniroot(function(h) sum(diag(nw_S(x, x, h))) - target, c(0.5, 5))$root
h_ll <- uniroot(function(h) sum(diag(ll_S(x, x, h))) - target, c(0.5, 5))$root
span <- uniroot(function(s) loess(accel ~ times, mcycle, span = s)$trace.hat - target, c(0.2, 0.5))$root
lo12 <- loess(accel ~ times, mcycle, span = span)
ss12 <- smooth.spline(x, y, df = target)
f4 <- cbind(NW = drop(nw_S(x, x, h_nw) %*% y), LL = drop(ll_S(x, x, h_ll) %*% y),
            LOESS = fitted(lo12), spline = predict(ss12, x)$y)
rss4 <- colSums((y - f4)^2)
data.frame(tuning = round(c(h_nw, h_ll, span, NA), 3), RSS = round(rss4),
           GCV = round(n * rss4 / (n - target)^2, 1), min_fit = round(apply(f4, 2, min), 1))
gd <- data.frame(times = rep(gridx, 4),
                 fit = c(drop(nw_S(gridx, x, h_nw) %*% y), drop(ll_S(gridx, x, h_ll) %*% y),
                         predict(lo12, data.frame(times = gridx)), predict(ss12, gridx)$y),
                 method = rep(c("NW", "local linear", "LOESS (quadratic)", "smoothing spline"), each = length(gridx)))
ggplot(gd, aes(times, fit, colour = method)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey65", size = 1, inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = ie_pal[1:4], name = NULL) +
  labs(x = "Time (ms)", y = "Acceleration (g)")

Sl <- ll_S(x, x, h_ll); fl <- drop(Sl %*% y)
sig_hat <- sqrt(sum((y - fl)^2) / (n - sum(diag(Sl))))
se_f <- sig_hat * sqrt(rowSums(Sl^2))
bd <- data.frame(times = x, f = fl, se = se_f, sdr = sd(y - fl))
ggplot(bd, aes(times)) +
  geom_ribbon(aes(ymin = f - 2 * sdr, ymax = f + 2 * sdr), fill = ie_coral, alpha = 0.18) +
  geom_ribbon(aes(ymin = f - 2 * se, ymax = f + 2 * se), fill = ie_navy, alpha = 0.35) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey55", size = 1, inherit.aes = FALSE) +
  geom_line(aes(y = f), colour = ie_navy, linewidth = 0.9) + labs(x = "Time (ms)", y = "Acceleration (g)")
c(sigma_hat = sig_hat, sd_resid = sd(y - fl), se_min = min(se_f), se_max = max(se_f))

xe <- c(-2, -1, 0, 1, 2); ye <- c(1, 3, 2, 5, 4)
rbind(orth_q2 = coef(lm(ye ~ poly(xe, 2))), orth_q3 = coef(lm(ye ~ poly(xe, 3)))[1:3])
rbind(raw_q2 = coef(lm(ye ~ poly(xe, 2, raw = TRUE))), raw_q3 = coef(lm(ye ~ poly(xe, 3, raw = TRUE)))[1:3])

Kt10 <- pmax(1 - abs(outer(xt, xt, "-")) / 10, 0)
S10 <- Kt10 / rowSums(Kt10)
c(edf = sum(diag(S10)), edf_h2 = sum(diag(St)))

x4 <- 1:4; y4 <- c(2, 3, 5, 4)
loc <- function(x0, h = 2) { w <- pmax(1 - abs(x4 - x0) / h, 0); d <- x4 - x0
  c(NW = sum(w * y4) / sum(w), LL = unname(coef(lm(y4 ~ d, weights = w))[1])) }
rbind(`x0 = 1` = loc(1), `x0 = 2.5` = loc(2.5))

B <- 1; V <- 0.25 / (2 * sqrt(pi))
res <- sapply(c(100, 1600), function(n) { h <- (V / (4 * B^2 * n))^(1 / 5); c(n = n, h_opt = h, AMSE = 5 * B^2 * h^4) })
round(t(res), 4)

s2 <- summary(lm(accel ~ poly(times, 12), data = mcycle))$sigma^2
cp <- sapply(1:15, function(q) sum(resid(lm(accel ~ poly(times, q), data = mcycle))^2) / n + 2 * s2 * (q + 1) / n)
round(cp, 1)
which.min(cp)

m5 <- lm(accel ~ poly(times, 5), data = mcycle)
loo_formula <- mean((resid(m5) / (1 - hatvalues(m5)))^2)
loo_brute <- mean(sapply(1:n, function(i) {
  mi <- lm(accel ~ poly(times, 5), data = mcycle[-i, ])
  (y[i] - predict(mi, mcycle[i, ]))^2
}))
c(formula = loo_formula, brute_force = loo_brute)

xk <- 1:6; yk <- c(2, 4, 3, 7, 9, 8); k <- 3
Sk <- t(sapply(xk, function(i) { w <- numeric(6); w[order(abs(xk - i), xk)[1:k]] <- 1 / k; w }))
fk <- drop(Sk %*% yk)
short <- (yk - fk) / (1 - diag(Sk))
loo <- sapply(1:6, function(i) { o <- setdiff(1:6, i); nb <- o[order(abs(xk[o] - xk[i]), xk[o])[1:k]]; yk[i] - mean(yk[nb]) })
rbind(shortcut = short, true_loo = loo)

data(exa, package = "faraway")
xe <- exa$x; ye <- exa$y; me <- exa$m; ne <- length(ye)
critE <- function(S) { r <- ye - drop(S %*% ye); edf <- sum(diag(S))
  c(edf = edf, LOOCV = mean((r / (1 - diag(S)))^2), GCV = ne * sum(r^2) / (ne - edf)^2) }
hE <- seq(0.01, 0.08, by = 0.002)
cN <- t(sapply(hE, function(h) critE(nw_S(xe, xe, h))))
cL <- t(sapply(hE, function(h) critE(ll_S(xe, xe, h))))
h_n <- hE[which.min(cN[, "LOOCV"])]; h_l <- hE[which.min(cL[, "LOOCV"])]
spE <- seq(0.05, 0.3, by = 0.01)
cLo <- sapply(spE, function(s) { m <- loess(ye ~ xe, span = s, degree = 2); ne * sum(resid(m)^2) / (ne - m$trace.hat)^2 })
s_lo <- spE[which.min(cLo)]
fE <- cbind(NW = drop(nw_S(xe, xe, h_n) %*% ye), LL = drop(ll_S(xe, xe, h_l) %*% ye),
            LOESS = fitted(loess(ye ~ xe, span = s_lo, degree = 2)))
reg <- ifelse(xe < 0.1 | xe > 0.9, "boundary", "interior")
mseE <- rbind(overall = colMeans((fE - me)^2),
              boundary = colMeans((fE[reg == "boundary", ] - me[reg == "boundary"])^2),
              interior = colMeans((fE[reg == "interior", ] - me[reg == "interior"])^2))
list(tuning = c(h_NW = h_n, h_LL = h_l, span = s_lo), mse = round(mseE * 1000, 2))

Sl <- ll_S(x, x, h_ll); fl <- drop(Sl %*% y)
sig_hat <- sqrt(sum((y - fl)^2) / (n - sum(diag(Sl))))
se_f <- sig_hat * sqrt(rowSums(Sl^2))
c(inside_2sd = mean(abs(y - fl) < 2 * sd(y - fl)), sigma_hat = sig_hat,
  se_min = min(se_f), se_max = max(se_f), halfwidth_naive = 2 * sd(y - fl))

cv_loess <- function(span, folds) {
  sse <- 0
  for (k in unique(folds)) {
    te <- folds == k
    m <- loess(accel ~ times, data = mcycle[!te, ], span = span, degree = 2,
               control = loess.control(surface = "direct"))
    sse <- sse + sum((mcycle$accel[te] - predict(m, mcycle[te, ]))^2)
  }
  sse / n
}
spans <- seq(0.15, 0.6, by = 0.05)
best_shuf <- sapply(1:10, function(s) { set.seed(s); f <- sample(rep(1:10, length.out = n))
  spans[which.min(sapply(spans, cv_loess, folds = f))] })
fsort <- rep(1:10, each = 14)[1:n]
cv_sorted <- sapply(spans, cv_loess, folds = fsort)
set.seed(1); f1 <- sample(rep(1:10, length.out = n)); cv_shuf1 <- sapply(spans, cv_loess, folds = f1)
list(shuffled_best_spans = best_shuf, shuffled_min_seed1 = min(cv_shuf1),
     sorted_best_span = spans[which.min(cv_sorted)], sorted_min = min(cv_sorted))
