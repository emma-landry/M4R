library(rstanarm)
library(epidemia)
library(dplyr)

data("EuropeCovid2")
data <- EuropeCovid2$data
head(data)
data <- filter(data, date > date[which(cumsum(deaths) > 10)[1] - 30])
data <- filter(data, date < as.Date("2020-05-05"))
rt <- epirt(
  formula = R(country, date) ~ 0 + (1 + public_events + schools_universities + self_isolating_if_ill + social_distancing_encouraged + lockdown || country) + public_events + schools_universities + self_isolating_if_ill + social_distancing_encouraged + lockdown,
  prior = shifted_gamma(shape=1/6, scale = 1, shift = log(1.05)/6),
  prior_covariance = decov(shape = c(2, rep(0.5, 5)),scale=0.25),
  link = scaled_logit(6.5)
)
inf <- epiinf(gen = EuropeCovid$si, seed_days = 6)
deaths <- epiobs(
  formula = deaths ~ 1,
  i2o = EuropeCovid2$inf2death,
  prior_intercept = normal(0,0.2),
  link = scaled_logit(0.02)
)



args <- list(rt=rt, inf=inf, obs=deaths, data=data, seed=12345)
args$algorithm <- "fullrank"
args$iter <- 1e4
args$tol_rel_obj <- 1e-3

#args$algorithm <- "sampling"
#args$chains <- 0

fm <- do.call(epim, args)