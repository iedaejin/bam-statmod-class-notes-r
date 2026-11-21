source("_common.R")

data(nhtemp)
nh <- data.frame(year = as.numeric(time(nhtemp)), temp = as.numeric(nhtemp))
n <- nrow(nh)
ggplot(nh, aes(year, temp)) + geom_line(colour = ie_navy) + geom_point(colour = ie_navy, size = 1.2) +
  labs(x = "Year", y = "Mean annual temperature (deg F)")
summary(nh$temp)

m1 <- lm(temp ~ 1, data = nh)
m2 <- lm(temp ~ year, data = nh)
ybar <- mean(nh$temp); s <- sd(nh$temp)
Stt <- sum((nh$year - mean(nh$year))^2); Sty <- sum((nh$year - mean(nh$year)) * (nh$temp - ybar))
c(ybar = ybar, s = s, se_ybar = s / sqrt(n))
c(Stt = Stt, Sty = Sty, b1 = Sty / Stt, b0 = ybar - Sty / Stt * mean(nh$year))
cbind(by_hand = c(ybar - Sty / Stt * mean(nh$year), Sty / Stt), lm = coef(m2))
coef(lm(temp ~ I(year - 1912), data = nh))          # intercept = expected temperature in 1912
10 * coef(m2)["year"]                                 # degrees F per decade

ggplot(nh, aes(year, temp)) + geom_point(colour = ie_navy, size = 1.3) +
  geom_hline(yintercept = ybar, linetype = 2, colour = ie_coral) +
  geom_line(aes(y = fitted(m2)), colour = ie_blue, linewidth = 1) +
  labs(x = "Year", y = "Temperature (deg F)")

rss1 <- sum(resid(m1)^2); rss2 <- sum(resid(m2)^2)
c(RSS1 = rss1, RSS2 = rss2, sigma1 = summary(m1)$sigma, sigma2 = summary(m2)$sigma, R2 = 1 - rss2 / rss1)
anova(m1, m2)

res_df <- rbind(data.frame(year = nh$year, e = resid(m1), model = "Model 1 (constant)"),
                data.frame(year = nh$year, e = resid(m2), model = "Model 2 (trend)"))
p1 <- ggplot(res_df, aes(year, e)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(colour = ie_navy, size = 1) + facet_wrap(~model) + labs(x = "Year", y = "Residual")
p2 <- ggplot(data.frame(e = resid(m2)), aes(sample = e)) + stat_qq(colour = ie_navy) +
  stat_qq_line(colour = ie_coral) + labs(x = "Theoretical quantiles", y = "Residuals")
a2 <- acf(resid(m2), plot = FALSE, lag.max = 15)
p3 <- ggplot(data.frame(lag = a2$lag[-1], acf = a2$acf[-1]), aes(lag, acf)) +
  geom_col(fill = ie_navy, width = 0.3) + geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(n), linetype = 2, colour = ie_coral) +
  labs(x = "Lag (years)", y = "ACF")
gridExtra::grid.arrange(p1, gridExtra::arrangeGrob(p2, p3, ncol = 2), ncol = 1)
c(lag1_m1 = acf(resid(m1), plot = FALSE)$acf[2], lag1_m2 = acf(resid(m2), plot = FALSE)$acf[2])
shapiro.test(resid(m2))

set.seed(2027)
pos <- sample(1:9, 1)
sig2 <- summary(m2)$sigma
lu <- do.call(rbind, lapply(1:9, function(k) {
  data.frame(panel = k, year = nh$year,
             temp = if (k == pos) nh$temp else fitted(m2) + rnorm(n, 0, sig2))
}))
ggplot(lu, aes(year, temp)) + geom_line(colour = ie_navy, linewidth = 0.4) + facet_wrap(~panel, nrow = 3) +
  labs(x = NULL, y = NULL) + theme(axis.text = element_text(size = 6))

set.seed(2027)
B <- 1000
sig1 <- summary(m1)$sigma
slope1 <- replicate(B, coef(lm(rnorm(n, ybar, sig1) ~ nh$year))[2])
slope2 <- replicate(B, coef(lm(fitted(m2) + rnorm(n, 0, sig2) ~ nh$year))[2])
b1_obs <- unname(coef(m2)[2])
c(observed = b1_obs, max_under_m1 = max(slope1), share_m1_at_least_obs = mean(slope1 >= b1_obs))
c(mean_m2 = mean(slope2), sd_m2 = sd(slope2), se_from_lm = summary(m2)$coefficients["year", "Std. Error"])
sl <- data.frame(slope = 10 * c(slope1, slope2), model = rep(c("Model 1: no trend", "Model 2: trend"), each = B))
ggplot(sl, aes(slope, fill = model)) + geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 10 * b1_obs, linewidth = 1) + scale_fill_manual(values = c(ie_coral, ie_navy)) +
  labs(x = "Slope (deg F per decade)", y = "Count", fill = NULL)

set.seed(2027)
ns <- c(10, 30, 100, 300, 1000)
vu <- t(sapply(ns, function(nn) {
  sims <- replicate(500, { z <- rnorm(nn, ybar, s); c(mean(z), sd(z)) })
  c(n = nn, variability = mean(sims[2, ]), uncertainty = sd(sims[1, ]), theory = s / sqrt(nn))
}))
round(vu, 3)
vl <- data.frame(n = rep(ns, 2), value = c(vu[, "variability"], vu[, "uncertainty"]),
                 what = rep(c("variability (SD of data)", "uncertainty (SE of mean)"), each = length(ns)))
ggplot(vl, aes(n, value, colour = what)) + geom_line() + geom_point() + scale_x_log10() +
  scale_colour_manual(values = c(ie_coral, ie_navy)) + labs(x = "Sample size n", y = "deg F", colour = NULL)

tq <- qt(0.975, n - 1)
ci1 <- ybar + c(-1, 1) * tq * s / sqrt(n); pi1 <- ybar + c(-1, 1) * tq * s * sqrt(1 + 1 / n)
rbind(CI_mu = ci1, PI_y0 = pi1)
new <- data.frame(year = 1972)
rbind(fit = predict(m2, new), predict(m2, new, interval = "confidence")[, -1], predict(m2, new, interval = "prediction")[, -1])

set.seed(2027)
N <- 1e5
x <- rnorm(N, 0, 1); yy <- 2 + 1 * x + rnorm(N, 0, 1); w <- x + rnorm(N, 0, 1)   # sigma_x^2 = sigma_u^2 = 1
c(slope_true_x = unname(coef(lm(yy ~ x))[2]), slope_noisy_x = unname(coef(lm(yy ~ w))[2]), theory = 1 / (1 + 1))
y_noisy <- yy + rnorm(N, 0, 1)
c(slope = unname(coef(lm(y_noisy ~ x))[2]), sigma_clean = summary(lm(yy ~ x))$sigma, sigma_noisy = summary(lm(y_noisy ~ x))$sigma)

data(airquality)
c(rows = nrow(airquality), complete = sum(complete.cases(airquality)))
colSums(is.na(airquality))
miss <- with(airquality, tapply(is.na(Ozone), Month, sum)); days <- with(airquality, tapply(Ozone, Month, length))
rbind(missing = miss, days = days)
oz_month <- with(airquality, tapply(Ozone, Month, mean, na.rm = TRUE))
c(complete_case = mean(airquality$Ozone, na.rm = TRUE), month_weighted = sum(oz_month * days) / sum(days))
nobs(lm(Ozone ~ Solar.R + Wind + Temp, data = airquality))     # lm silently drops incomplete rows

sel <- subset(nh, temp >= 50)
c(n_all = nrow(nh), n_selected = nrow(sel))
c(mean_all = mean(nh$temp), mean_selected = mean(sel$temp))
c(slope_all = 10 * unname(coef(lm(temp ~ year, data = nh))[2]), slope_selected = 10 * unname(coef(lm(temp ~ year, data = sel))[2]))

train <- subset(nh, year <= 1961); test <- subset(nh, year > 1961)
h1 <- lm(temp ~ 1, data = train); h2 <- lm(temp ~ year, data = train)
rmse <- function(fit) sqrt(mean((test$temp - predict(fit, test))^2))
c(rmse_model1 = rmse(h1), rmse_model2 = rmse(h2))
c(slope_train = unname(coef(h2)[2]), se = summary(h2)$coef[2, 2])
ggplot(nh, aes(year, temp)) + annotate("rect", xmin = 1961.5, xmax = 1971.5, ymin = -Inf, ymax = Inf, fill = "#D2E6F4") +
  geom_point(colour = ie_navy, size = 1.2) +
  geom_line(data = test, aes(y = predict(h1, test), colour = "Model 1")) +
  geom_line(data = test, aes(y = predict(h2, test), colour = "Model 2")) +
  geom_line(data = train, aes(y = fitted(h2), colour = "Model 2")) + geom_line(data = train, aes(y = fitted(h1), colour = "Model 1")) +
  scale_colour_manual(values = c(ie_coral, ie_blue)) + labs(x = "Year", y = "deg F", colour = NULL)

set.seed(2027); a <- rnorm(3)
set.seed(2027); b <- rnorm(3)
identical(a, b)                       # same seed, same numbers
identical(a, rnorm(3))                # no reset: different numbers
# A check that a stored copy of the data equals the package data
nh_csv <- tempfile(fileext = ".csv"); write.csv(nh, nh_csv, row.names = FALSE)
all.equal(read.csv(nh_csv)$temp, as.numeric(nhtemp))
R.version.string

# # A typical layout (relative paths from the project root):
# # project/
# #   data/nhtemp.csv        # raw data, never edited by hand
# #   R/01-fit-models.R      # analysis scripts
# #   report.Rmd             # report: reads data/, sources R/
# #   output/                # figures and tables written by code
# nh <- read.csv("data/nhtemp.csv")      # relative path: works on any machine
# stopifnot(nrow(nh) == 60, all(names(nh) == c("year", "temp")))
# sessionInfo()                          # last line of every report

yy2 <- nh$temp; i <- which.max(yy2); yy3 <- yy2; yy3[i] <- yy3[i] + 10
c(mean_change = mean(yy3) - mean(yy2), median_change = median(yy3) - median(yy2))

se0 <- 1.266 / sqrt(60); f <- 1 + 2 * 0.3 * 59 / 60
c(se0 = se0, se_dep = se0 * sqrt(f), ratio = sqrt(f), coverage = 2 * pnorm(1.96 / sqrt(f)) - 1)

c(se = 1.266 / sqrt(60), n_needed = (1.266 / 0.05)^2, limit_halfwidth = qnorm(0.975) * 1.266)

xo <- 1:5; yo <- c(2.1, 3.9, 6.2, 7.8, 10.1)
fo <- lm(yo ~ 0 + xo); c(beta = coef(fo), RSS = sum(resid(fo)^2), sigma2 = sigma(fo)^2, se = summary(fo)$coef[1, 2])

Stt <- sum((nh$year - mean(nh$year))^2); sg <- sigma(m2)
c(Stt = Stt, sigma = sg, se_b1 = sg / sqrt(Stt), se_b0 = sg * sqrt(1 / n + mean(nh$year)^2 / Stt))
summary(m2)$coef[, 2]
m2c <- lm(temp ~ I(year - 1912), data = nh); summary(m2c)$coef[, 1:2]

r <- cor(nh$year, nh$temp); R2 <- r^2
c(r = r, R2 = R2, R2_lm = summary(m2)$r.squared, F = R2 * (n - 2) / (1 - R2), t2 = summary(m2)$coef[2, 3]^2)

set.seed(2027); N <- 2e5
xs <- rnorm(N, 0, 2); ys <- 1 + 2 * xs + rnorm(N); ws <- xs + rnorm(N)
c(slope_on_w = unname(coef(lm(ys ~ ws))[2]), theory = 1.6, slope_on_x = unname(coef(lm(ys ~ xs))[2]))

set.seed(2027); z <- rnorm(1e6); obs <- z[z <= 1]
c(simulated = mean(obs), theory = -dnorm(1) / pnorm(1), share_observed = mean(z <= 1))

set.seed(2027); N <- 1e5
xs <- rnorm(N); ys <- 1 + 2 * xs + rnorm(N)
c(all = unname(coef(lm(ys ~ xs))[2]), x_positive = unname(coef(lm(ys[xs > 0] ~ xs[xs > 0]))[2]),
  y_above_1 = unname(coef(lm(ys[ys > 1] ~ xs[ys > 1]))[2]))

c(full_b1 = unname(coef(m2)[2]), full_ratio = unname(coef(m2)[2]) / (sigma(m2) / sqrt(Stt)))
Stt_tr <- sum((train$year - mean(train$year))^2)
c(train_b1 = unname(coef(h2)[2]), train_ratio = unname(coef(h2)[2]) / (sigma(h2) / sqrt(Stt_tr)))

set.seed(2027)
tt <- nh$year; b1_true <- unname(coef(m2)[2]); mu_t <- fitted(m2); sg <- sigma(m2)
cover <- function(ar) mean(replicate(2000, {
  e <- if (ar == 0) rnorm(n, 0, sg) else as.numeric(arima.sim(list(ar = ar), n, sd = sg * sqrt(1 - ar^2)))
  f <- lm((mu_t + e) ~ tt); ci <- confint(f)[2, ]
  ci[1] <= b1_true && b1_true <= ci[2]
}))
c(iid = cover(0), ar1 = cover(0.5))

predict(m2, data.frame(year = 2070), interval = "prediction")
predict(m2, data.frame(year = 2070), interval = "confidence")
