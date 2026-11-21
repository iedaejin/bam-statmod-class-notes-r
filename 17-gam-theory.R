source("_common.R")

library(mgcv)
library(MASS)
data(airquality)
aq <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp")])
aq$ly <- log(aq$Ozone)
n <- nrow(aq)
# total penalty matrix S_lambda of a fitted gam, padded with zeros to the full length of coef()
Slam_full <- function(g) {
  p <- length(coef(g)); S <- matrix(0, p, p)
  for (sm in g$smooth) {
    ii <- sm$first.para:sm$last.para
    for (l in seq_along(sm$S)) S[ii, ii] <- S[ii, ii] + g$sp[sm$first.sp + l - 1] * sm$S[[l]]
  }
  S
}

sm_raw <- smoothCon(s(Temp, k = 10), data = aq, absorb.cons = FALSE)[[1]]
dim(sm_raw$X); dim(sm_raw$S[[1]])
ev <- eigen(sm_raw$S[[1]], symmetric = TRUE)$values
round(ev, 4)
sum(ev < max(ev) * 1e-10)       # dimension of the null space

f0 <- gam(ly ~ s(Temp), data = aq, method = "REML")
xg <- data.frame(Temp = seq(min(aq$Temp), max(aq$Temp), length.out = 200))
d <- do.call(rbind, lapply(c(1e-6, f0$sp, 1e8), function(l) {
  gg <- gam(ly ~ s(Temp), data = aq, sp = l)
  data.frame(Temp = xg$Temp, fit = predict(gg, xg),
             lab = sprintf("lambda = %s (edf %.1f)", format(l, digits = 2), sum(gg$edf) - 1))
}))
ggplot(d, aes(Temp, fit, colour = lab)) +
  geom_point(aes(Temp, ly), data = aq, inherit.aes = FALSE, colour = "grey60") +
  geom_line(linewidth = 1) + scale_colour_manual(values = c(ie_coral, ie_navy, ie_green)) +
  labs(x = "Temp (F)", y = "log(Ozone)", colour = NULL)

X_raw <- cbind(1, sm_raw$X)
c(columns = ncol(X_raw), rank = qr(X_raw)$rank)           # intercept + 10 basis columns
sum(resid(lm(rep(1, n) ~ sm_raw$X - 1))^2)                # the constant is in the span

sm_con <- smoothCon(s(Temp, k = 10), data = aq, absorb.cons = TRUE)[[1]]
X_con <- cbind(1, sm_con$X)
c(columns = ncol(X_con), rank = qr(X_con)$rank)
range(colSums(sm_con$X))                          # every column sums to (numerically) 0
sum(eigen(sm_con$S[[1]], symmetric = TRUE)$values < 1e-8)   # null space is now 1-dimensional

pima <- rbind(Pima.tr, Pima.te)
mb <- gam(type ~ s(glu) + s(bmi) + s(age) + s(ped), family = binomial,
          data = pima, method = "REML")
X <- model.matrix(mb); yb <- as.numeric(pima$type == "Yes")
Sl <- Slam_full(mb)
beta <- rep(0, ncol(X))
for (it in 1:50) {
  eta <- drop(X %*% beta); mu <- plogis(eta)
  w <- mu * (1 - mu); z <- eta + (yb - mu) / w
  beta_new <- drop(solve(crossprod(X, w * X) + Sl, crossprod(X, w * z)))
  done <- max(abs(beta_new - beta)) < 1e-10; beta <- beta_new
  if (done) break
}
c(iterations = it, max_abs_difference = max(abs(beta - coef(mb))))

g3 <- gam(ly ~ s(Solar.R) + s(Wind) + s(Temp), data = aq, method = "REML")
X3 <- model.matrix(g3); S3 <- Slam_full(g3)
F3 <- solve(crossprod(X3) + S3, crossprod(X3))
c(trF = sum(diag(F3)), sum_edf = sum(g3$edf), trA = sum(g3$hat))
range(diag(F3))
round(range(Re(eigen(F3)$values)), 6)
sapply(g3$smooth, function(sm) sum(diag(F3)[sm$first.para:sm$last.para]))
summary(g3)$edf

c(edf = sum(g3$edf), edf1 = sum(g3$edf1), edf2 = sum(g3$edf2), AIC_df = attr(logLik(g3), "df"))
round(cbind(edf = summary(g3)$edf, Ref.df = summary(g3)$s.table[, "Ref.df"]), 3)
set.seed(1); xs <- runif(100); zs <- runif(100); ys <- sin(2 * pi * xs) + zs + rnorm(100, sd = 0.5)
gs <- gam(ys ~ s(xs) + s(zs), method = "REML")
c(edf = sum(gs$edf), edf1 = sum(gs$edf1), edf2 = sum(gs$edf2), AIC_df = attr(logLik(gs), "df"))

# brute-force leave-one-out with the design and penalty of the REML fit g3 (fixed lambda)
Xf <- X3; Sf <- S3; yy <- aq$ly
loo <- sapply(1:n, function(i) {
  b <- solve(crossprod(Xf[-i, ]) + Sf, crossprod(Xf[-i, ], yy[-i]))
  yy[i] - sum(Xf[i, ] * b)
})
A3 <- Xf %*% solve(crossprod(Xf) + Sf, t(Xf))
short <- (yy - A3 %*% yy) / (1 - diag(A3))
max(abs(loo - short))                                   # Lemma: identical
c(OCV = mean(loo^2), GCV = n * sum((yy - A3 %*% yy)^2) / (n - sum(diag(A3)))^2)

reml_ml <- function(g, y) {
  X <- model.matrix(g); S <- Slam_full(g); b <- coef(g); n <- length(y)
  P <- sum((y - X %*% b)^2) + drop(t(b) %*% S %*% b)
  eg <- eigen(S, symmetric = TRUE); pos <- eg$values > max(eg$values) * 1e-10; M <- sum(!pos)
  H <- crossprod(X) + S; Um <- eg$vectors[, pos]; Hpp <- crossprod(Um, H %*% Um)
  ldS <- sum(log(eg$values[pos])); ldH <- determinant(H)$modulus; ldHpp <- determinant(Hpp)$modulus
  s2R <- P / (n - M); s2M <- P / n
  c(M = M, sigma2_REML = s2R, sigma2_ML = s2M,
    REML = 0.5 * (P / s2R + (n - M) * log(2 * pi * s2R) + ldH - ldS),
    ML   = 0.5 * (P / s2M + n * log(2 * pi * s2M) + ldHpp - ldS))
}
fR <- ly ~ s(Solar.R) + s(Wind) + s(Temp)
gR <- gam(fR, data = aq, method = "REML"); gM <- gam(fR, data = aq, method = "ML")
round(rbind(REML_fit = reml_ml(gR, aq$ly), ML_fit = reml_ml(gM, aq$ly)), 4)
c(mgcv_REML_fit = gR$gcv.ubre, mgcv_ML_fit = gM$gcv.ubre)

mu <- plogis(drop(X %*% coef(mb)))
ll <- sum(dbinom(yb, 1, mu, log = TRUE))
H  <- crossprod(X, (mu * (1 - mu)) * X) + Sl
eg <- eigen(Sl, symmetric = TRUE); pos <- eg$values > max(eg$values) * 1e-10; M <- sum(!pos)
V_laml <- -(ll - 0.5 * drop(t(coef(mb)) %*% Sl %*% coef(mb))) + 0.5 * determinant(H)$modulus -
  0.5 * sum(log(eg$values[pos])) - M / 2 * log(2 * pi)
c(by_hand = as.numeric(V_laml), mgcv = mb$gcv.ubre, M = M)

f0fun <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * (10 * x)^3 * (1 - x)^10
set.seed(5); nrep <- 200; n_s <- 60; xg <- data.frame(x = seq(0, 1, length.out = 200))
sim <- do.call(rbind, lapply(1:nrep, function(i) {
  x <- runif(n_s); y <- f0fun(x) + rnorm(n_s, sd = 2)
  do.call(rbind, lapply(c(GCV = "GCV.Cp", REML = "REML", ML = "ML"), function(m) {
    g <- gam(y ~ s(x, k = 20), method = m)
    data.frame(edf = sum(g$edf) - 1, mse = mean((predict(g, xg) - f0fun(xg$x))^2))
  })) -> r
  r$method <- rownames(r); r$rep <- i; r
}))
sim$method <- factor(sim$method, levels = c("GCV", "REML", "ML"))
tab <- do.call(rbind, lapply(split(sim, sim$method), function(d)
  data.frame(edf_mean = mean(d$edf), edf_q90 = quantile(d$edf, .9), edf_max = max(d$edf),
             share_edf_above_12 = mean(d$edf > 12), mse_mean = mean(d$mse), mse_q90 = quantile(d$mse, .9))))
round(tab, 3)

ggplot(sim, aes(method, edf, fill = method)) + geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c(ie_coral, ie_navy, ie_green)) + labs(x = NULL, y = "edf of s(x)")

Xt <- model.matrix(f0); St <- f0$smooth[[1]]$S[[1]]; pt <- ncol(Xt); yt <- aq$ly
crit <- t(sapply(seq(-7, 14, by = 0.2), function(L) {
  Sl_ <- matrix(0, pt, pt); Sl_[2:pt, 2:pt] <- exp(L) * St
  Hm <- crossprod(Xt) + Sl_; b <- solve(Hm, crossprod(Xt, yt)); rss <- sum((yt - Xt %*% b)^2)
  P <- rss + drop(t(b) %*% Sl_ %*% b); trA <- sum(diag(solve(Hm, crossprod(Xt))))
  ev <- eigen(Sl_, symmetric = TRUE, only.values = TRUE)$values; pos <- ev > max(ev) * 1e-10; Mt <- sum(!pos)
  Um <- eigen(Sl_, symmetric = TRUE)$vectors[, pos]; Hpp <- crossprod(Um, Hm %*% Um)
  s2R <- P / (n - Mt); s2M <- P / n
  c(L = L, GCV = n * rss / (n - trA)^2,
    REML = 0.5 * (P / s2R + (n - Mt) * log(2 * pi * s2R) + determinant(Hm)$modulus - sum(log(ev[pos]))),
    ML = 0.5 * (P / s2M + n * log(2 * pi * s2M) + determinant(Hpp)$modulus - sum(log(ev[pos]))))
}))
cd <- as.data.frame(crit)
dd <- rbind(data.frame(L = cd$L, score = cd$GCV / min(cd$GCV), crit = "GCV (relative to minimum)"),
            data.frame(L = cd$L, score = cd$REML - min(cd$REML), crit = "REML (minus minimum)"),
            data.frame(L = cd$L, score = cd$ML - min(cd$ML), crit = "ML (minus minimum)"))
ggplot(dd, aes(L, score)) + geom_line(colour = ie_navy, linewidth = 1) +
  facet_wrap(~crit, scales = "free_y") + labs(x = expression(log(lambda)), y = NULL)
c(GCV = cd$L[which.min(cd$GCV)], REML = cd$L[which.min(cd$REML)], ML = cd$L[which.min(cd$ML)])

nd <- data.frame(Solar.R = 207, Wind = 9.7, Temp = c(70, 90))
Xp <- predict(g3, nd, type = "lpmatrix")
se_hand <- sqrt(rowSums((Xp %*% vcov(g3)) * Xp))
rbind(by_hand = se_hand, predict = predict(g3, nd, se.fit = TRUE)$se.fit)

set.seed(6)
cover <- replicate(200, {
  x <- runif(100); y <- f0fun(x) + rnorm(100, sd = 1.5)
  g <- gam(y ~ s(x, k = 20), method = "REML")
  p <- predict(g, xg, se.fit = TRUE)
  as.numeric(abs(p$fit - f0fun(xg$x)) <= 1.96 * p$se.fit)
})
c(average_over_function = mean(cover))
round(quantile(rowMeans(cover), c(0, .05, .5, .95, 1)), 3)

wg <- seq(min(aq$Wind), max(aq$Wind), length.out = 100)
ndw <- data.frame(Solar.R = 207, Wind = wg, Temp = 79)
Xw <- predict(g3, ndw, type = "lpmatrix")
idx <- g3$smooth[[2]]$first.para:g3$smooth[[2]]$last.para
set.seed(2024)
br <- rmvn(1000, coef(g3), vcov(g3))
fsim <- Xw[, idx] %*% t(br[, idx]); fhat <- drop(Xw[, idx] %*% coef(g3)[idx])
se <- sqrt(rowSums((Xw[, idx] %*% vcov(g3)[idx, idx]) * Xw[, idx]))
inside <- abs(fsim - fhat) <= 1.96 * se
c(pointwise = mean(inside), whole_curve = mean(apply(inside, 2, all)))
crit_sim <- quantile(apply(abs(fsim - fhat) / se, 2, max), .95)
crit_sim                                                          # simultaneous multiplier (instead of 1.96)
matplot(wg, fhat + fsim[, 1:25] - fhat, type = "l", lty = 1, col = "grey60",
        xlab = "Wind (mph)", ylab = "s(Wind), centred")
lines(wg, fhat, lwd = 3, col = ie_navy)

idw <- idx
r <- sqrt(diag(vcov(g3, unconditional = TRUE)))[idw] / sqrt(diag(vcov(g3)))[idw]
round(range(r), 3)

sb <- summary(mb)
Vp <- vcov(mb); bb <- coef(mb)
tab <- t(sapply(1:4, function(j) {
  jj <- mb$smooth[[j]]$first.para:mb$smooth[[j]]$last.para
  Xj <- X[, jj]; fj <- Xj %*% bb[jj]
  ev <- eigen(Xj %*% Vp[jj, jj] %*% t(Xj), symmetric = TRUE)
  r <- round(sb$edf[j])
  c(edf = sb$edf[j], rank = r, T_by_hand = sum((t(ev$vectors[, 1:r, drop = FALSE]) %*% fj)^2 / ev$values[1:r]),
    Chi.sq_mgcv = sb$s.table[j, "Chi.sq"])
}))
rownames(tab) <- rownames(sb$s.table); round(tab, 2)

set.seed(17)
pv_reml <- replicate(500, { x <- runif(100); yy <- rnorm(100)
  summary(gam(yy ~ s(x), method = "REML"))$s.table[1, "p-value"] })
pv_gcv  <- replicate(500, { x <- runif(100); yy <- rnorm(100)
  summary(gam(yy ~ s(x)))$s.table[1, "p-value"] })
rbind(REML = c(`p<0.05` = mean(pv_reml < .05), `p<0.01` = mean(pv_reml < .01)),
      GCV  = c(`p<0.05` = mean(pv_gcv < .05),  `p<0.01` = mean(pv_gcv < .01)))

set.seed(24)
xs <- runif(200); ys <- 2 * sin(4 * pi * xs) + rnorm(200, sd = .5); ds <- data.frame(x = xs, y = ys)
k5  <- gam(y ~ s(x, k = 5),  data = ds, method = "REML")
k20 <- gam(y ~ s(x, k = 20), data = ds, method = "REML")
set.seed(1); k.check(k5)
set.seed(1); k.check(k20)
xx <- data.frame(x = seq(0, 1, length.out = 200))
dk <- rbind(data.frame(x = xx$x, fit = predict(k5, xx), k = "k = 5"),
            data.frame(x = xx$x, fit = predict(k20, xx), k = "k = 20"))
ggplot(dk, aes(x, fit, colour = k)) + geom_point(aes(x, y), data = ds, inherit.aes = FALSE, colour = "grey70") +
  stat_function(fun = function(x) 2 * sin(4 * pi * x), colour = "black", linetype = 2) +
  geom_line(linewidth = 1) + scale_colour_manual(values = c(ie_coral, ie_navy)) + labs(x = "x", y = "y", colour = NULL)

set.seed(1)
gam.check(g3)

S <- matrix(0, 4, 4)
for (k in 2:3) for (l in 2:3)       # b_k'' = k(k-1) x^(k-2), exponents 0-indexed
  S[k + 1, l + 1] <- integrate(function(x) k * (k - 1) * x^(k - 2) * l * (l - 1) * x^(l - 2), 0, 1)$value
S

c(lambda1 = 111 * 36.0 / (111 - 6)^2, lambda2 = 111 * 34.5 / (111 - 8)^2)

vu <- function(rss, trA, s2 = .30, n = 111) rss / n - s2 + 2 * s2 * trA / n
c(fit1 = vu(36.0, 6), fit2 = vu(34.5, 8))

B <- splines::bs(aq$Temp, df = 5, intercept = TRUE)
C <- matrix(colSums(B), 1); Z <- qr.Q(qr(t(C)), complete = TRUE)[, -1]
c(rank_raw = qr(cbind(1, B))$rank, cols_raw = ncol(cbind(1, B)),
  rank_con = qr(cbind(1, B %*% Z))$rank, cols_con = ncol(cbind(1, B %*% Z)),
  max_CZ = max(abs(C %*% Z)))
coef(lm(aq$ly ~ cbind(1, B %*% Z) - 1))[1]; mean(aq$ly)

s <- c(0, 0, 1, 3, 9); lam <- 2
sum(1 / (1 + lam * s))                       # 1 + 1 + 1/3 + 1/7 + 1/19
sapply(c(0, 1e-8, 1e8), function(l) sum(1 / (1 + l * s)))

round(reml_ml(gR, aq$ly)[c("M", "sigma2_REML", "sigma2_ML")], 4)

nd <- data.frame(Solar.R = median(aq$Solar.R), Wind = median(aq$Wind), Temp = 85)
Xp <- predict(g3, nd, type = "lpmatrix"); est <- drop(Xp %*% coef(g3))
se1 <- sqrt(drop(Xp %*% vcov(g3) %*% t(Xp)))
se2 <- sqrt(drop(Xp %*% vcov(g3, unconditional = TRUE) %*% t(Xp)))
round(c(est = est, lo = est - 1.96 * se1, hi = est + 1.96 * se1,
        predict_se = predict(g3, nd, se.fit = TRUE)$se.fit, se_hand = se1,
        se_unconditional = se2, width_ratio = se2 / se1), 3)

fits <- lapply(c(GCV = "GCV.Cp", REML = "REML", ML = "ML"), function(m) gam(ly ~ s(Temp), data = aq, method = m))
round(rbind(sp = sapply(fits, function(g) g$sp), edf = sapply(fits, function(g) sum(g$edf) - 1)), 3)
tg <- data.frame(Temp = seq(min(aq$Temp), max(aq$Temp), length.out = 200))
dm <- do.call(rbind, lapply(names(fits), function(m) data.frame(tg, fit = predict(fits[[m]], tg), method = m)))
ggplot(dm, aes(Temp, fit, colour = method)) + geom_point(aes(Temp, ly), data = aq, inherit.aes = FALSE, colour = "grey70") +
  geom_line(linewidth = 1) + scale_colour_manual(values = c(ie_coral, ie_navy, ie_green)) +
  labs(x = "Temp (F)", y = "log(Ozone)", colour = NULL)

set.seed(1)
res <- t(sapply(c(5, 10, 20, 40), function(k) {
  g <- gam(y ~ s(x, k = k), data = ds, method = "REML")
  kc <- k.check(g); c(k = k, kprime = kc[1, 1], edf = kc[1, 2], kindex = kc[1, 3], p = kc[1, 4],
                      sigma = sqrt(g$sig2))
}))
round(res, 3)
