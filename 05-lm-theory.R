source("_common.R")

data(savings, package = "faraway")
str(savings)

library(tidyr)
long <- pivot_longer(savings, c(pop15, pop75, dpi, ddpi), names_to = "predictor", values_to = "x")
ggplot(long, aes(x, sr)) +
  geom_point(colour = ie_navy, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = ie_coral, linewidth = 0.7) +
  facet_wrap(~ predictor, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Savings rate sr (%)")

fit <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = savings)
X <- model.matrix(fit)
y <- savings$sr
dim(X)
head(X, 3)

XtX <- crossprod(X); Xty <- crossprod(X, y)
beta_hat <- solve(XtX, Xty)           # solves the normal equations
cbind(by_hand = drop(beta_hat), lm = coef(fit))
e <- y - X %*% beta_hat
max(abs(crossprod(X, e)))             # X'e = 0 up to rounding
sum(e)                                # residuals sum to 0 (intercept in the model)

poly <- data.frame(x = c(0, 4.2, 5.4, 1.2), y = c(0, 0, 1.5, 1.5))
ggplot() +
  geom_polygon(data = poly, aes(x, y), fill = "#E0EEF8") +
  annotate("segment", x = 1.2, y = 0.5, xend = 3.4, yend = 2.9, colour = ie_navy, linewidth = 0.9, arrow = arrow(length = unit(0.18, "cm"))) +
  annotate("segment", x = 1.2, y = 0.5, xend = 3.4, yend = 0.9, colour = ie_navy, linewidth = 0.9, arrow = arrow(length = unit(0.18, "cm"))) +
  annotate("segment", x = 3.4, y = 0.9, xend = 3.4, yend = 2.9, colour = ie_coral, linewidth = 0.9, arrow = arrow(length = unit(0.18, "cm"))) +
  annotate("text", x = 3.55, y = 3.0, label = "y", colour = ie_navy, size = 5) +
  annotate("text", x = 3.8, y = 0.55, label = "hat(y) == Hy", parse = TRUE, colour = ie_navy, size = 4.5) +
  annotate("text", x = 3.9, y = 1.9, label = "e == (I - H) * y", parse = TRUE, colour = ie_coral, size = 4.5) +
  annotate("text", x = 0.8, y = 0.2, label = "C(X)", colour = ie_navy, size = 4.5) +
  coord_equal(xlim = c(-0.2, 6), ylim = c(-0.2, 3.3)) + theme_void()

summary(fit)

n <- nrow(X); p <- ncol(X)
RSS <- sum(resid(fit)^2)
sigma2_hat <- RSS / (n - p)
V <- sigma2_hat * solve(crossprod(X))          # estimated Var(beta_hat)
se <- sqrt(diag(V))
tval <- coef(fit) / se
pval <- 2 * pt(abs(tval), n - p, lower.tail = FALSE)
round(cbind(estimate = coef(fit), se = se, t = tval, p = pval), 4)
c(RSS = RSS, sigma_hat = sqrt(sigma2_hat))
coef(fit)["pop15"] + c(-1, 1) * qt(0.975, n - p) * se["pop15"]     # CI by hand
confint(fit)["pop15", ]

fit0 <- lm(sr ~ pop15 + ddpi, data = savings)     # M0: drops pop75 and dpi
anova(fit0, fit)
RSS0 <- sum(resid(fit0)^2)
Fstat <- ((RSS0 - RSS) / 2) / (RSS / (n - p))
c(F = Fstat, p.value = pf(Fstat, 2, n - p, lower.tail = FALSE))

sigma2_ml <- RSS / n
c(by_hand = -n / 2 * (log(2 * pi * sigma2_ml) + 1), logLik = as.numeric(logLik(fit)))
c(LR = n * log(RSS0 / RSS), from_F = n * log(1 + 2 / (n - p) * Fstat),
  chisq_p = pchisq(n * log(RSS0 / RSS), 2, lower.tail = FALSE))

qrX <- qr(X)
beta_qr <- backsolve(qr.R(qrX), crossprod(qr.Q(qrX), y))
cbind(qr = drop(beta_qr), normal_eq = drop(beta_hat), lm = coef(fit))
H <- X %*% solve(crossprod(X), t(X))
c(symmetric = max(abs(H - t(H))), idempotent = max(abs(H %*% H - H)), trace = sum(diag(H)))
c(kappa_X = kappa(X, exact = TRUE), kappa_XtX = kappa(crossprod(X), exact = TRUE))

# share of the population aged 15-75 is the complement of the other two shares
sav <- transform(savings, pop_mid = 100 - pop15 - pop75)
fit_def <- lm(sr ~ pop15 + pop75 + pop_mid, data = sav)
coef(fit_def)   # pop_mid is aliased: pop15 + pop75 + pop_mid = 100 * (intercept column)

set.seed(505)
beta_true <- coef(fit); sigma_true <- 4; B <- 4000
sim <- replicate(B, {
  ysim <- drop(X %*% beta_true) + rnorm(n, sd = sigma_true)
  f <- lm(ysim ~ X - 1)
  ci <- confint(f)["Xpop15", ]
  c(rss = sum(resid(f)^2), b = coef(f)["Xpop15"],
    cover = unname(ci[1] < beta_true["pop15"] & beta_true["pop15"] < ci[2]))
})
c(mean_RSS = mean(sim["rss", ]), theory = sigma_true^2 * (n - p),
  mean_RSS_over_n = mean(sim["rss", ]) / n, mean_RSS_over_np = mean(sim["rss", ]) / (n - p),
  coverage = mean(sim["cover", ]))
c(cor_beta_RSS = cor(sim["b.Xpop15", ], sim["rss", ]),
  var_beta = var(sim["b.Xpop15", ]), theory = sigma_true^2 * solve(crossprod(X))["pop15", "pop15"])

x5 <- 1:5; y5 <- c(2.1, 3.9, 6.2, 7.8, 10.1)
f5 <- lm(y5 ~ x5); summary(f5)$coefficients; c(RSS = deviance(f5), R2 = summary(f5)$r.squared)
summary(f5)$fstatistic

m0 <- lm(sr ~ pop15 + ddpi, data = savings); m1 <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = savings)
a <- anova(m0, m1); LR <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
c(F = a$F[2], p_F = a$`Pr(>F)`[2], LR = LR, p_LR = pchisq(LR, 2, lower.tail = FALSE))

data(trees)
Xt <- cbind(1, trees$Girth, trees$Height); yt <- trees$Volume
b_hand <- solve(crossprod(Xt), crossprod(Xt, yt))
ft <- lm(Volume ~ Girth + Height, data = trees)
cbind(hand = drop(b_hand), lm = coef(ft))
summary(ft)$coefficients; confint(ft)["Girth", ]
flog <- lm(log(Volume) ~ log(Girth) + log(Height), data = trees)
summary(flog)$coefficients; confint(flog)
c(R2_linear = summary(ft)$r.squared, R2_log = summary(flog)$r.squared)

set.seed(11)
one <- function(skew) {
  x <- runif(20, 0, 10)
  eps <- if (skew) sqrt(2) * (rnorm(20)^2 - 1) else rnorm(20, sd = 2)
  f <- lm(I(1 + 0.5 * x + eps) ~ x)
  ci <- confint(f)["x", ]
  c(s2 = sigma(f)^2, s2ml = deviance(f) / 20, cover = unname(ci[1] < 0.5 & 0.5 < ci[2]))
}
rbind(normal = rowMeans(replicate(5000, one(FALSE))), skewed = rowMeans(replicate(5000, one(TRUE))))

d <- data.frame(g = factor(rep(c("A", "B", "C"), each = 4)), y = c(4, 5, 6, 5, 7, 8, 9, 8, 3, 2, 4, 3))
Xd <- cbind(1, model.matrix(~ g - 1, d)); qr(Xd)$rank
coef(lm(y ~ g, d)); coef(lm(y ~ I(1 * (g == "A")) + I(1 * (g == "B")) + I(1 * (g == "C")), d))   # intercept + all three dummies
tapply(d$y, d$g, mean)
