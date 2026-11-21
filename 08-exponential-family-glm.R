source("_common.R")
suppressPackageStartupMessages({library(MASS); library(faraway)})

bw <- with(birthwt, data.frame(low = low, age = age, lwt = lwt,
  race = factor(race, labels = c("white", "black", "other")),
  smoke = smoke, ptd = as.integer(ptl > 0), ht = ht, ui = ui))
table(bw$low)
sh <- subset(ships, service > 0)          # 6 of the 40 cells have no service
sh$year <- factor(sh$year); sh$period <- factor(sh$period)
c(cells = nrow(sh), incidents = sum(sh$incidents), months = sum(sh$service))

lpm <- lm(low ~ age + lwt + race + smoke + ptd + ht + ui, data = bw)
logit_fit <- glm(low ~ age + lwt + race + smoke + ptd + ht + ui, family = binomial, data = bw)
c(negative = sum(fitted(lpm) < 0), min = min(fitted(lpm)), max = max(fitted(lpm)))
range(fitted(logit_fit))

d <- rbind(data.frame(model = "Linear model", fitted = fitted(lpm)),
           data.frame(model = "Logistic GLM", fitted = fitted(logit_fit)))
ggplot(d, aes(fitted)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = ie_navy, colour = "white") +
  geom_vline(xintercept = c(0, 1), colour = ie_coral, linetype = 2) +
  facet_wrap(~model, ncol = 1) + labs(x = "fitted value of P(low = 1)", y = "births")

d1 <- function(b, th, h = 1e-4) (b(th + h) - b(th - h)) / (2 * h)
d2 <- function(b, th, h = 1e-4) (b(th + h) - 2 * b(th) + b(th - h)) / h^2
c(mean = d1(exp, log(3)), var = d2(exp, log(3)))                  # Poisson, mu = 3
b_gamma <- function(th) -log(-th)                                    # gamma, mu = 2, nu = 4
c(mean = d1(b_gamma, -1/2), var = d2(b_gamma, -1/2) * (1/4))      # a(phi) = 1/nu
set.seed(1); x <- rgamma(1e5, shape = 4, rate = 2)               # mean 2, variance 1
c(mean = mean(x), var = var(x))

eta <- seq(-4, 4, length.out = 300)
d <- data.frame(eta = rep(eta, 3), mu = c(plogis(eta), pnorm(eta), 1 - exp(-exp(eta))),
                link = rep(c("logit", "probit", "cloglog"), each = 300))
ggplot(d, aes(eta, mu, colour = link)) + geom_line(linewidth = 1) +
  scale_colour_manual(values = ie_pal[1:3]) + labs(x = expression(eta), y = expression(mu), colour = NULL)

set.seed(2027)
mus <- c(1, 2, 4, 6, 8, 10)
sim_var <- function(fam, mu) var(switch(fam,
  Normal = rnorm(4000, mu, sqrt(2)), Poisson = rpois(4000, mu),
  Gamma = rgamma(4000, shape = 2, rate = 2 / mu)))
pts <- expand.grid(family = c("Normal", "Poisson", "Gamma"), mu = mus, stringsAsFactors = FALSE)
pts$var <- mapply(sim_var, pts$family, pts$mu)
mm <- seq(0.5, 10.5, length.out = 100)
lines_df <- rbind(data.frame(family = "Normal", mu = mm, var = 2),
                  data.frame(family = "Poisson", mu = mm, var = mm),
                  data.frame(family = "Gamma", mu = mm, var = mm^2 / 2))
ggplot(pts, aes(mu, var, colour = family)) + geom_line(data = lines_df, linewidth = 1) + geom_point(size = 2.5) +
  scale_colour_manual(values = ie_pal[c(3, 1, 2)]) + labs(x = expression(mu), y = "Var(Y)", colour = NULL)

set.seed(2027)
res <- replicate(200, { x <- runif(200); y <- rpois(200, exp(1 + 0.5 * x))
  c(`lm(log(y+1))` = unname(coef(lm(log(y + 1) ~ x))[2]),
    `Poisson glm` = unname(coef(glm(y ~ x, family = poisson))[2])) })
round(rbind(mean = rowMeans(res), sd = apply(res, 1, sd), truth = c(0.5, 0.5)), 3)

o <- glm(cbind(damage, 6 - damage) ~ temp, family = binomial, data = orings)
round(summary(o)$coefficients, 4)

pr <- predict(o, data.frame(temp = 31), type = "link", se.fit = TRUE)
c(prob = plogis(pr$fit), lower = plogis(pr$fit - 1.96 * pr$se.fit), upper = plogis(pr$fit + 1.96 * pr$se.fit))
range(orings$temp)

fs <- glm(incidents ~ type + year + period + offset(log(service)), family = poisson, data = sh)
round(coef(summary(fs)), 4)

fg <- glm(Volume ~ log(Girth) + log(Height), family = Gamma(link = "log"), data = trees)
fl <- lm(log(Volume) ~ log(Girth) + log(Height), data = trees)
round(rbind(`gamma, log link` = coef(fg), `lm on log(y)` = coef(fl)), 3)
c(phi_hat = summary(fg)$dispersion, nu_hat = 1 / summary(fg)$dispersion)
round(confint(fg), 3)

fam <- poisson()
c(fam$family, fam$link)
fam$linkfun(3); fam$linkinv(log(3)); fam$mu.eta(log(3)); fam$variance(3)
binomial()$linkfun(0.3); Gamma()$link; binomial(link = "probit")$link

# glm(low ~ smoke + age, family = binomial, data = bw)                     # 0/1 response
# glm(cbind(damage, 6 - damage) ~ temp, family = binomial, data = orings)  # counts out of 6
# glm(damage / 6 ~ temp, family = binomial, weights = rep(6, 23), data = orings)  # same fit
# glm(incidents ~ type + offset(log(service)), family = poisson, data = sh)
# glm(Volume ~ log(Girth), family = Gamma(link = "log"), data = trees)

y <- c(2, 3, 6, 7, 10, 12, 15); x <- c(0, 0, 0, 0, 1, 1, 1)
coef(glm(y ~ x, family = poisson))
c(log(mean(y[x == 0])), log(mean(y[x == 1]) / mean(y[x == 0])))

fg <- glm(Volume ~ log(Girth) + log(Height), family = Gamma(link = "log"), data = trees)
fl <- lm(log(Volume) ~ log(Girth) + log(Height), data = trees)
est <- coef(summary(fg))["log(Girth)", ]
c(z = (est[1] - 2) / est[2])
nd <- data.frame(Girth = 15, Height = 75)
c(gamma = predict(fg, nd, type = "response"), lm_backtransformed = exp(predict(fl, nd)),
  lm_smearing = exp(predict(fl, nd)) * mean(exp(residuals(fl))))

b <- coef(o); V <- vcov(o)
ld50 <- -b[1] / b[2]
g <- c(-1 / b[2], b[1] / b[2]^2)                    # gradient of -b0/b1
c(ld50 = ld50, se = sqrt(drop(t(g) %*% V %*% g)))
(1 - plogis(b[1] + b[2] * 70))^6

sim <- function(b0, reps = 200) {
  set.seed(2027)
  res <- replicate(reps, { x <- runif(200); y <- rpois(200, exp(b0 + 0.5 * x))
    c(unname(coef(lm(log(y + 1) ~ x))[2]), unname(coef(glm(y ~ x, family = poisson))[2])) })
  setNames(rowMeans(res), c("lm(log(y+1))", "Poisson glm"))
}
round(rbind(`intercept 1` = sim(1), `intercept 4` = sim(4)), 3)
