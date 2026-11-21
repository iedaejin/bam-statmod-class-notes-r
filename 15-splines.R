source("_common.R")

library(MASS)
library(splines)
library(mgcv)
data(mcycle)
x <- mcycle$times
y <- mcycle$accel
n <- length(y)

gridx <- seq(min(x), max(x), length.out = 300)
rs <- do.call(rbind, lapply(c(3, 8, 20), function(K) {
  m <- lm(accel ~ bs(times, df = K + 3), data = mcycle)
  data.frame(times = gridx, fit = predict(m, data.frame(times = gridx)), K = paste(K, "knots"))
}))
rs$K <- factor(rs$K, levels = paste(c(3, 8, 20), "knots"))
ggplot(rs, aes(times, fit)) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey60", size = 1) +
  geom_line(colour = ie_coral, linewidth = 0.9) + facet_wrap(~K, nrow = 1) +
  coord_cartesian(ylim = c(-150, 100)) + labs(x = "Time (ms)", y = "Acceleration (g)")

dfs <- c(5, 8, 10, 12, 15)
aic <- sapply(dfs, function(df) c(
  bs = AIC(lm(accel ~ bs(times, df = df), data = mcycle)),
  ns = AIC(lm(accel ~ ns(times, df = df), data = mcycle))))
colnames(aic) <- paste0("df=", dfs)
round(aic, 1)

tpower <- function(x, knots, q) outer(x, knots, function(u, k) pmax(u - k, 0)^q)
knots8 <- quantile(x, probs = (1:8) / 9)
X_tp <- cbind(outer(x, 0:3, `^`), tpower(x, knots8, 3))
X_bs <- cbind(1, bs(x, knots = knots8, degree = 3))
c(dim = dim(X_tp), kappa_TP = kappa(X_tp, exact = TRUE), kappa_B = kappa(X_bs, exact = TRUE),
  max_diff_fitted = max(abs(fitted(lm(y ~ X_tp - 1)) - fitted(lm(y ~ X_bs - 1)))))

splineDesign(-3:8, 1.5, ord = 4)[, 2:5]
c(1, 23, 23, 1) / 48

bbase <- function(x, xl = min(x) - 0.01, xr = max(x) + 0.01, nseg = 20, bdeg = 3) {
  dx <- (xr - xl) / nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)       # bdeg extra segments at each end
  B <- splineDesign(knots, x, ord = bdeg + 1, outer.ok = TRUE)
  attr(B, "knots") <- knots
  B
}
B <- bbase(x)
c(dim(B), range(rowSums(B)), max(rowSums(B > 0)))

Bg <- bbase(gridx)
Dm <- diff(diag(ncol(B)), differences = 2)
th <- solve(crossprod(B) + 0.64 * crossprod(Dm), crossprod(B, y))
bl <- data.frame(times = rep(gridx, ncol(B)), j = rep(1:ncol(B), each = length(gridx)),
                 val = as.vector(Bg), wval = as.vector(sweep(Bg, 2, th, "*")))
ggplot(bl, aes(times, wval, group = j)) + geom_line(colour = ie_blue, linewidth = 0.4) +
  geom_line(data = data.frame(times = gridx, fit = drop(Bg %*% th)), aes(times, fit, group = NULL),
            colour = ie_coral, linewidth = 1) +
  geom_point(data = mcycle, aes(times, accel, group = NULL), colour = "grey60", size = 0.9) +
  labs(x = "Time (ms)", y = expression(beta[j] * B[j](x)))

tk <- attr(B, "knots"); dl <- diff(tk)[1]
thv <- drop(th)
f1_direct <- drop(splineDesign(tk, gridx, ord = 4, derivs = 1, outer.ok = TRUE) %*% thv)
f1_diff <- drop(splineDesign(tk[-c(1, length(tk))], gridx, ord = 3, outer.ok = TRUE) %*% (diff(thv) / dl))
c(dim_Bq1 = ncol(splineDesign(tk[-c(1, length(tk))], gridx, ord = 3, outer.ok = TRUE)),
  max_abs_diff = max(abs(f1_direct - f1_diff)))

gx2 <- seq(0, 75, length.out = 300)
cmp <- do.call(rbind, lapply(c("bs", "ns"), function(b) {
  m <- if (b == "bs") lm(accel ~ bs(times, df = 8), mcycle) else lm(accel ~ ns(times, df = 8), mcycle)
  data.frame(times = gx2, fit = suppressWarnings(predict(m, data.frame(times = gx2))), basis = b)
}))
ggplot(cmp, aes(times, fit, colour = basis)) + geom_vline(xintercept = range(x), linetype = 3) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey65", size = 1, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) + scale_colour_manual(values = c(ie_coral, ie_navy), name = NULL) +
  coord_cartesian(ylim = c(-200, 100)) + labs(x = "Time (ms)", y = "Acceleration (g)")

ss <- smooth.spline(x, y)                         # lambda chosen by GCV
ss_cv <- suppressWarnings(smooth.spline(x, y, cv = TRUE))
rbind(GCV = c(df = ss$df, lambda = ss$lambda), LOOCV = c(df = ss_cv$df, lambda = ss_cv$lambda))
c(sum_leverage = sum(ss$lev), df = ss$df)

crossprod(diff(diag(5), differences = 2))

psfit <- function(B, y, lambda, d = 2) {
  D <- diff(diag(ncol(B)), differences = d)
  A <- crossprod(B) + lambda * crossprod(D)          # B'B + lambda D'D
  theta <- solve(A, crossprod(B, y))
  f <- drop(B %*% theta); n <- length(y)
  edf <- sum(diag(solve(A, crossprod(B)))); rss <- sum((y - f)^2)
  list(theta = drop(theta), f = f, edf = edf, rss = rss,
       gcv = n * rss / (n - edf)^2, aic = n * log(rss / n) + 2 * edf, sigma = sqrt(rss / (n - edf)))
}
sapply(c(1e-8, 0.01, 1, 100, 1e4, 1e10), function(l) psfit(B, y, l)$edf)

ylin <- 3 + 2 * x                                    # exactly linear data
errs <- sapply(1:3, function(d) max(abs(psfit(B, ylin, 1000, d = d)$f - ylin)))
names(errs) <- paste0("d=", 1:3); signif(errs, 3)
f2 <- psfit(B, y, 0.64, d = 2)$f
rbind(data = c(sum(y), sum(x * y), sum(x^2 * y)), fit = c(sum(f2), sum(x * f2), sum(x^2 * f2)))

ll <- seq(-6, 8, by = 0.25)
ggplot(data.frame(ll, edf = sapply(ll, function(l) psfit(B, y, 10^l)$edf)), aes(ll, edf)) +
  geom_line(colour = ie_navy, linewidth = 1) + geom_hline(yintercept = c(2, 23), linetype = 3) +
  labs(x = expression(log[10](lambda)), y = "edf")

loglam <- seq(-3, 6, by = 0.05)
crit <- t(sapply(loglam, function(l) { f <- psfit(B, y, 10^l); c(GCV = f$gcv, AIC = f$aic, edf = f$edf) }))
opt <- optimize(function(l) psfit(B, y, 10^l)$gcv, c(-3, 6))
fit <- psfit(B, y, 10^opt$minimum)
c(log10_lambda = opt$minimum, lambda = 10^opt$minimum, edf = fit$edf, GCV = fit$gcv, sigma_hat = fit$sigma,
  lambda_AIC_grid = 10^loglam[which.min(crit[, "AIC"])])
ggplot(data.frame(loglam, GCV = crit[, "GCV"]), aes(loglam, GCV)) + geom_line(colour = ie_navy, linewidth = 1) +
  geom_vline(xintercept = opt$minimum, linetype = 2) + coord_cartesian(ylim = c(550, 1500)) +
  labs(x = expression(log[10](lambda)), y = "GCV")

g <- gam(accel ~ s(times, bs = "ps", k = ncol(B), m = c(2, 2)), data = mcycle,
         knots = list(times = attr(B, "knots")), method = "GCV.Cp")
c(edf_mgcv = sum(g$edf), edf_hand = fit$edf, max_abs_diff_fit = max(abs(fitted(g) - fit$f)))
g10 <- gam(accel ~ s(times, bs = "ps"), data = mcycle, method = "GCV.Cp")   # default k = 10
c(edf_default_k = sum(g10$edf), edf_k23 = sum(g$edf))

A <- crossprod(B) + 10^opt$minimum * crossprod(diff(diag(ncol(B)), differences = 2))
S_l <- B %*% solve(A, t(B))
se_bayes <- fit$sigma * sqrt(diag(S_l)); se_freq <- fit$sigma * sqrt(rowSums(S_l^2))
rbind(bayes = range(se_bayes), frequentist = range(se_freq))
bd <- data.frame(times = x, f = fit$f, lo = fit$f - 2 * se_bayes, hi = fit$f + 2 * se_bayes)
ggplot(bd, aes(times)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = ie_navy, alpha = 0.25) +
  geom_point(data = mcycle, aes(times, accel), colour = "grey60", size = 1) +
  geom_line(aes(y = f), colour = ie_navy, linewidth = 1) + labs(x = "Time (ms)", y = "Acceleration (g)")

q <- 3; K <- 30
knots30 <- quantile(unique(x), probs = (1:K) / (K + 1))
xs <- x / max(x); ks <- knots30 / max(x)                      # rescale for conditioning
Xt <- cbind(outer(xs, 0:q, `^`), tpower(xs, ks, q)); Pt <- diag(c(rep(0, q + 1), rep(1, K)))
tp_fit <- function(lambda) {
  A <- crossprod(Xt) + lambda * Pt; f <- drop(Xt %*% solve(A, crossprod(Xt, y)))
  edf <- sum(diag(solve(A, crossprod(Xt)))); c(edf = edf, GCV = n * sum((y - f)^2) / (n - edf)^2)
}
cr <- t(sapply(10^seq(-10, 0, by = 0.05), tp_fit))
c(cr[which.min(cr[, "GCV"]), ], edf_huge_lambda = unname(tp_fit(1e8)["edf"]))

library(nlme)
cc <- ncol(B); D2 <- diff(diag(cc), differences = 2)
G <- cbind(1, 1:cc)                                          # null space of D2: constants and linear
Xf <- B %*% G; Zr <- B %*% t(D2) %*% solve(D2 %*% t(D2))
dd <- data.frame(y = y, g = factor(1))
m_lme <- lme(y ~ Xf - 1, random = list(g = pdIdent(~ Zr - 1)), data = dd, method = "REML")
s2 <- m_lme$sigma^2; sb2 <- as.numeric(VarCorr(m_lme)[2, 1]); lam_reml <- s2 / sb2
g_reml <- gam(accel ~ s(times, bs = "ps", k = cc, m = c(2, 2)), data = mcycle,
              knots = list(times = attr(B, "knots")), method = "REML")
c(sigma2 = s2, sigma_b2 = sb2, lambda_REML = lam_reml, edf_P_spline = psfit(B, y, lam_reml)$edf,
  edf_mgcv_REML = sum(g_reml$edf), max_diff_fit = max(abs(fitted(m_lme) - psfit(B, y, lam_reml)$f)))

tt <- 0:7; th4 <- c(2, 5, 1, 4)
c(f = sum(splineDesign(tt, 3.7, ord = 4) * th4), fprime = sum(splineDesign(tt, 3.7, ord = 4, derivs = 1) * th4))

M <- crossprod(B); P <- crossprod(diff(diag(cc), differences = 2))
eM <- eigen(M, symmetric = TRUE); Mih <- eM$vectors %*% diag(1 / sqrt(eM$values)) %*% t(eM$vectors)
gam_k <- eigen(Mih %*% P %*% Mih, symmetric = TRUE, only.values = TRUE)$values
c(smallest = sort(gam_k)[1:3], edf_formula = sum(1 / (1 + 0.64 * gam_k)), edf_psfit = psfit(B, y, 0.64)$edf)

xx <- seq(0, 1, length.out = 100); yq <- 1 + 2 * xx + 3 * xx^2
err <- sapply(1:3, function(bd) sapply(1:4, function(d)
  max(abs(psfit(bbase(xx, bdeg = bd), yq, 1e3, d)$f - yq))))
dimnames(err) <- list(paste0("d=", 1:4), paste0("q=", 1:3)); signif(err, 2)

mom <- t(sapply(1:4, function(d) { f <- psfit(B, y, 0.64, d)$f
  sapply(0:3, function(k) sum(x^k * (y - f)) / sum(abs(x^k * y))) }))
dimnames(mom) <- list(paste0("d=", 1:4), paste0("k=", 0:3)); signif(mom, 2)

Pm <- crossprod(diff(diag(cc), differences = 2))
ll <- seq(-3, 3, by = 0.1)
crit_ps <- t(sapply(10^ll, function(l) {
  A <- crossprod(B) + l * Pm; S <- B %*% solve(A, t(B)); r <- y - drop(S %*% y)
  c(GCV = n * sum(r^2) / (n - sum(diag(S)))^2, LOOCV = mean((r / (1 - diag(S)))^2))
}))
grp <- split(seq_len(n), x)
cross <- lapply(grp, function(i) list(i = i, BtB = crossprod(B[-i, ]), Bty = crossprod(B[-i, ], y[-i])))
gcv_grp <- sapply(10^ll, function(l) {
  sum(sapply(cross, function(g) { th <- solve(g$BtB + l * Pm, g$Bty); sum((y[g$i] - B[g$i, , drop = FALSE] %*% th)^2) })) / n })
c(GCV = ll[which.min(crit_ps[, "GCV"])], LOOCV = ll[which.min(crit_ps[, "LOOCV"])], grouped = ll[which.min(gcv_grp)])
c(GCV = min(crit_ps[, "GCV"]), LOOCV = min(crit_ps[, "LOOCV"]), grouped = min(gcv_grp))

res <- lapply(c(10, 20, 40), function(ns) {
  Bn <- bbase(x, nseg = ns); o <- optimize(function(l) psfit(Bn, y, 10^l)$gcv, c(-4, 8)); f <- psfit(Bn, y, 10^o$minimum)
  list(row = c(nseg = ns, columns = ncol(Bn), log10_lambda = o$minimum, edf = f$edf, GCV = f$gcv), f = f$f) })
tab9 <- do.call(rbind, lapply(res, `[[`, "row")); f20 <- res[[2]]$f
cbind(round(tab9, 2), max_diff_vs_20 = round(sapply(res, function(r) max(abs(r$f - f20))), 2))

lam12 <- 10^uniroot(function(l) psfit(B, y, 10^l)$edf - 12, c(-3, 4))$root
ps12 <- psfit(B, y, lam12); ss12 <- smooth.spline(x, y, df = 12)
fit_ss <- predict(ss12, x)$y
c(lambda_PS = lam12, RSS_PS = ps12$rss, RSS_SS = sum((y - fit_ss)^2), max_diff = max(abs(fit_ss - ps12$f)),
  sd_y = sd(y))

D2 <- diff(diag(cc), differences = 2); G2 <- cbind(1, 1:cc)
Xf <- B %*% G2; Zr <- B %*% t(D2) %*% solve(D2 %*% t(D2)); lam <- 0.64
H <- rbind(cbind(crossprod(Xf), crossprod(Xf, Zr)), cbind(crossprod(Zr, Xf), crossprod(Zr) + lam * diag(cc - 2)))
sol <- solve(H, c(crossprod(Xf, y), crossprod(Zr, y)))
f_mm <- drop(Xf %*% sol[1:2] + Zr %*% sol[-(1:2)])
max(abs(f_mm - psfit(B, y, lam)$f))

ed <- sapply(c(10, 23, 40), function(k) sum(gam(accel ~ s(times, bs = "ps", k = k), data = mcycle, method = "GCV.Cp")$edf))
names(ed) <- paste0("k=", c(10, 23, 40)); round(ed, 2)
round(k.check(gam(accel ~ s(times, bs = "ps", k = 10), data = mcycle, method = "GCV.Cp")), 3)
