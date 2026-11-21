source("_common.R")

y <- as.numeric(discoveries)
tab <- table(factor(y, levels = 0:12))
dd <- data.frame(count = 0:12, observed = as.numeric(tab) / length(y),
                 poisson = dpois(0:12, mean(y)))
p1 <- ggplot(data.frame(year = 1860:1959, y), aes(year, y)) +
  geom_col(fill = ie_navy, width = 0.8) + labs(x = "Year", y = "Discoveries")
p2 <- ggplot(dd, aes(count)) + geom_col(aes(y = observed), fill = "#D2E6F4") +
  geom_line(aes(y = poisson), colour = ie_coral) + geom_point(aes(y = poisson), colour = ie_coral) +
  labs(x = "Count", y = "Relative frequency")
gridExtra::grid.arrange(p1, p2, ncol = 2)
c(n = length(y), mean = mean(y), variance = var(y), ratio = var(y) / mean(y))

d <- data.frame(
  Distribution = c("Bernoulli$(\\pi)$", "Binomial$(m,\\pi)$", "Poisson$(\\lambda)$",
                   "Neg. binomial$(\\mu,k)$", "Normal$(\\mu,\\sigma^2)$", "Exponential$(\\lambda)$",
                   "Gamma$(\\alpha,\\beta)$", "Chi-squared$(\\nu)$"),
  Support = c("$\\{0,1\\}$", "$\\{0,\\dots,m\\}$", "$\\{0,1,\\dots\\}$", "$\\{0,1,\\dots\\}$",
              "$\\mathbb{R}$", "$(0,\\infty)$", "$(0,\\infty)$", "$(0,\\infty)$"),
  Mean = c("$\\pi$", "$m\\pi$", "$\\lambda$", "$\\mu$", "$\\mu$", "$1/\\lambda$", "$\\alpha/\\beta$", "$\\nu$"),
  Variance = c("$\\pi(1-\\pi)$", "$m\\pi(1-\\pi)$", "$\\lambda$", "$\\mu+\\mu^2/k$", "$\\sigma^2$",
               "$1/\\lambda^2$", "$\\alpha/\\beta^2$", "$2\\nu$"),
  R = c("`dbinom(y,1,p)`", "`dbinom(y,m,p)`", "`dpois(y,lambda)`", "`dnbinom(y,mu=,size=k)`",
        "`dnorm(y,mean,sd)`", "`dexp(y,rate)`", "`dgamma(y,shape,rate)`", "`dchisq(y,df)`"))
knitr::kable(d, format = "latex", booktabs = TRUE, escape = FALSE,
             caption = "Distributions used in the course, with R names. Rate $\\beta$ means density $\\propto y^{\\alpha-1}e^{-\\beta y}$.",
             linesep = "") |>
  kableExtra::kable_styling(latex_options = "hold_position", font_size = 9)

mu_hat <- mean(y); s2 <- var(y)
k_hat <- mu_hat^2 / (s2 - mu_hat)
c(mu = mu_hat, k = k_hat)
# Expected relative frequencies under the two models
rbind(poisson = dpois(0:3, mu_hat), negbin = dnbinom(0:3, mu = mu_hat, size = k_hat),
      observed = as.numeric(table(factor(y, levels = 0:3))) / length(y))

dpois(3, lambda = 3.1)                       # P(Y = 3)
ppois(3, lambda = 3.1)                       # P(Y <= 3)
qpois(0.5, lambda = 3.1)                     # median
set.seed(1); rpois(5, 3.1)                   # five random counts
pnorm(1.96); qnorm(0.975)                    # normal: Phi and Phi^{-1}
pgamma(2, shape = 3, rate = 2)               # rate parametrisation
pgamma(2, shape = 3, scale = 1/2)            # same, with scale = 1/rate
pexp(1, rate = 2) - pgamma(1, shape = 1, rate = 2)   # exponential is Gamma(1, rate): 0
integrate(function(t) t^2 * dgamma(t, 3, 2), 0, Inf)$value   # E(Y^2) = alpha(alpha+1)/beta^2 = 3

mu <- c(1, 2, 3)
Sigma <- matrix(c(4, 2, 0,  2, 5, 1,  0, 1, 3), 3, byrow = TRUE)
A <- matrix(c(1, -1, 0,  1, 1, 1), 2, byrow = TRUE)
A %*% mu
A %*% Sigma %*% t(A)
set.seed(2027)
Ysim <- MASS::mvrnorm(1e5, mu, Sigma)          # 100000 x 3 matrix, one draw per row
AY <- Ysim %*% t(A)                            # apply A to every draw
colMeans(AY); cov(AY)

set.seed(2027)
rho <- 0.6; s1 <- 15; s2 <- 10; m1 <- 100; m2 <- 50
S <- matrix(c(s1^2, rho*s1*s2, rho*s1*s2, s2^2), 2)
xy <- as.data.frame(MASS::mvrnorm(2000, c(m1, m2), S)); names(xy) <- c("y1", "y2")
b <- rho * s1 / s2
ggplot(xy, aes(y2, y1)) + geom_point(alpha = 0.3, colour = ie_navy, size = 0.8) +
  geom_abline(intercept = m1 - b * m2, slope = b, colour = ie_coral) +
  geom_abline(intercept = m1 - (s1/s2) * m2, slope = s1/s2, linetype = 2)
# Condition on Y2 in a narrow window around 65
sub <- subset(xy, abs(y2 - 65) < 1)
c(n = nrow(sub), mean = mean(sub$y1), theory_mean = m1 + b * (65 - m2),
  sd = sd(sub$y1), theory_sd = s1 * sqrt(1 - rho^2))

data(mtcars)
X_ok  <- model.matrix(~ wt + hp + am, data = mtcars)
mtcars$manual <- 1 - mtcars$am
X_bad <- cbind(X_ok, manual = mtcars$manual)
c(p_ok = ncol(X_ok), rank_ok = qr(X_ok)$rank, p_bad = ncol(X_bad), rank_bad = qr(X_bad)$rank)
# the null vector v: intercept - am - manual = 0
range(X_bad %*% c(1, 0, 0, -1, -1))
coef(lm(mpg ~ wt + hp + am + manual, data = mtcars))   # R drops one column: NA

X <- model.matrix(~ wt + hp, data = mtcars); ymp <- mtcars$mpg
XtX <- crossprod(X)                    # t(X) %*% X, faster and symmetric
Xty <- crossprod(X, ymp)                # t(X) %*% y
H <- X %*% solve(XtX) %*% t(X)         # hat matrix (inverting is fine for a demonstration)
c(trace = sum(diag(H)), idempotent = max(abs(H %*% H - H)))
eigen(XtX)$values                      # all positive: X'X is positive definite
chol(XtX)                              # upper triangular R with R'R = X'X
qr(X)$rank                             # numerical rank
kappa(X, exact = TRUE)                 # condition number: large values signal near-collinearity

fit <- lm(mpg ~ wt + hp, data = mtcars)
beta_hat <- solve(XtX, Xty)                 # solves X'X b = X'y without inverting
cbind(by_hand = drop(beta_hat), lm = coef(fit))
e <- ymp - X %*% beta_hat
crossprod(X, e)                             # X'e = 0 (up to rounding): residuals orthogonal to columns
rss <- sum(e^2); c(RSS = rss, sigma2_hat = rss / (nrow(X) - ncol(X)))

f <- function(x) x^2 * exp(-x)
taylor2 <- function(x) 4 * exp(-2) - exp(-2) * (x - 2)^2
xs <- seq(0, 5, length.out = 300)
ggplot(data.frame(xs), aes(xs)) + stat_function(fun = f, colour = ie_navy, linewidth = 1) +
  stat_function(fun = taylor2, colour = ie_coral, linetype = 2) + coord_cartesian(ylim = c(-0.2, 0.7)) +
  labs(x = "x", y = "f(x)")
c(max_value = f(2), four_exp_minus2 = 4 * exp(-2), f2_exact = f(2.3), f2_taylor = taylor2(2.3))

g <- function(p) p[1] * exp(-(p[1]^2 + p[2]^2))
opt <- optim(c(0.3, 0.3), g, method = "BFGS", control = list(fnscale = -1), hessian = TRUE)
rbind(numerical = opt$par, exact = c(1/sqrt(2), 0))
opt$hessian                                   # negative definite
eigen(opt$hessian)$values

set.seed(2027)
sim <- sapply(c(2, 5, 20, 80), function(lambda) {
  ybar <- replicate(5000, mean(rpois(30, lambda)))
  c(lambda = lambda, var_sqrt = var(sqrt(ybar)), delta = 1 / (4 * 30))
})
round(t(sim), 5)

f1 <- function(x) (2 * x - x^2) * exp(-x)           # f'
f2 <- function(x) (2 - 4 * x + x^2) * exp(-x)       # f''
newton <- function(x0, iter = 6) {
  x <- numeric(iter + 1); x[1] <- x0
  for (t in seq_len(iter)) x[t + 1] <- x[t] - f1(x[t]) / f2(x[t])
  x
}
rbind(from_1.5 = newton(1.5), from_3 = newton(3))

S <- function(l) sum(y) / l - length(y)
I <- function(l) sum(y) / l^2
l <- 1; path <- l
for (t in 1:8) { l <- l + S(l) / I(l); path <- c(path, l) }
round(path, 6)
c(lambda_hat = mean(y), se = sqrt(mean(y) / length(y)), se_from_I = 1 / sqrt(I(mean(y))))

optimize(function(x) x^2 * exp(-x), c(0, 6), maximum = TRUE)           # 2 and 4 exp(-2)
uniroot(function(x) x^2 - 2, c(1, 2), tol = 1e-10)$root                # sqrt(2)
D(quote(x^2 * exp(-x)), "x")                                           # symbolic derivative
# Rosenbrock's banana function: minimum at (1, 1)
rosen <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
optim(c(-1.2, 1), rosen, method = "BFGS")[c("par", "value", "counts")]

c(poisson = dpois(0, mean(y)), negbin = dnbinom(0, mu = mean(y), size = mean(y)^2 / (var(y) - mean(y))),
  observed = mean(y == 0))

Sg <- matrix(c(4,1,0, 1,3,1, 0,1,2), 3, byrow = TRUE)
Am <- rbind(c(1,1,0), c(1,-1,0), c(0,1,1))
Am %*% c(2, 0, 1); Am %*% Sg %*% t(Am)

set.seed(2027)
Sg2 <- matrix(c(100, 0.6*80, 0.6*80, 64), 2)
d2 <- MASS::mvrnorm(4e5, c(70, 170), Sg2)
w <- d2[abs(d2[, 2] - 186) < 1, 1]
c(n = length(w), mean = mean(w), sd = sd(w), p85 = mean(w > 85),
  theory_p85 = 1 - pnorm(0.375), unconditional = 1 - pnorm(1.5))

D <- (length(y) - 1) * var(y) / mean(y)
c(D = D, df = length(y) - 1, p_value = pchisq(D, length(y) - 1, lower.tail = FALSE))

Xc <- cbind(1, sapply(c(4, 6, 8), function(v) as.numeric(mtcars$cyl == v)))
c(p = ncol(Xc), rank = qr(Xc)$rank)
dd4 <- data.frame(mpg = mtcars$mpg, c4 = Xc[, 2], c6 = Xc[, 3], c8 = Xc[, 4])
coef(lm(mpg ~ c4 + c6 + c8, data = dd4))

Xh <- cbind(1, 1:3); Hh <- Xh %*% solve(crossprod(Xh)) %*% t(Xh)
Hh * 6; drop(Hh %*% c(1, 3, 2))

yy <- mtcars$mpg
b_ridge <- solve(crossprod(X_bad) + 1e-4 * diag(ncol(X_bad)), crossprod(X_bad, yy))
fit_ok <- fitted(lm(mpg ~ wt + hp + am, data = mtcars))
c(max_abs_diff_fitted = max(abs(X_bad %*% b_ridge - fit_ok)))
round(drop(b_ridge), 3)

fc <- function(p) p[1]^3 - 3 * p[1] * p[2] + p[2]^3
o <- optim(c(1.3, 0.7), fc, method = "BFGS", hessian = TRUE)
round(c(x = o$par[1], y = o$par[2], value = o$value), 4); round(eigen(o$hessian)$values, 3)

x <- 1; for (t in 1:3) { x <- (x + 2 / x) / 2; cat(sprintf("x%d = %.7f  error = %.2e\n", t, x, x - sqrt(2))) }

n <- length(y); ybar <- mean(y)
ci_wald <- ybar + c(-1, 1) * 1.96 * sqrt(ybar / n)
ci_log  <- exp(log(ybar) + c(-1, 1) * 1.96 / sqrt(n * ybar))
rbind(wald = ci_wald, log_scale = ci_log)
