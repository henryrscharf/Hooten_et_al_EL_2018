## libraries ----
library(smover); library(spatstat)
library(reDyn)
## load mtl data ----
load("../data/mountain_lion/mtlAM80.RData")
ml.df <- ml.df[251:400, ] ## JASA paper subset
## medium rasters ----
ELE.rast <- raster("../data/mountain_lion/mtlELE.grd")
ELE.rast <- crop(ELE.rast, extent(c(435000, 478000, 4347000, 4397000)))
## covariates ----
killsite <- c(458404, 4384135)
dist2kill <- raster(ext = extent(ELE.rast), res = res(ELE.rast), crs = crs(ELE.rast))
dist2kill <- distanceFromPoints(object = dist2kill, xy = killsite)
dist2kill <- dist2kill
names(dist2kill) <- "dist2kill"
X <- stack(list(
  ELE.rast,
  terrain(ELE.rast, opt = "slope"),
  terrain(ELE.rast, opt = "aspect"),
  terrain(ELE.rast, opt = "roughness"),
  dist2kill))
X_scaled <- scale(X)
gradX <- sapply(1:nlayers(X_scaled), function(layer) ctmcmove::rast.grad(X_scaled[[layer]])
                [c('rast.grad.x', 'rast.grad.y')], simplify = F)
gradX_norm_consts <- unlist(lapply(gradX, function(gradX_i){
  sqrt(sum(values(gradX_i$rast.grad.x)^2, values(gradX_i$rast.grad.y)^2) /
         (2 * length(gradX_i$rast.grad.x) - 1))
}))
gradX_scaled <- sapply(1:length(gradX), function(layer){
  list("rast.grad.x" = gradX[[layer]]$rast.grad.x / gradX_norm_consts[layer],
       "rast.grad.y" = gradX[[layer]]$rast.grad.y / gradX_norm_consts[layer])
}, simplify = F)

intercept <- raster(X_scaled)
values(intercept) <- 1
W <- stack(list(
  "intercept" = intercept,
  "killsite" = X$dist2kill < 1000, ## up to a kilometer away
  "elevation" = X_scaled$elevation,
  "slope" = X_scaled$slope,
  "elevation_slope" = X_scaled$elevation * X_scaled$slope
))
rm(X)
## set up model fit to data ----
TIMES <- nrow(ml.df); FITTIMES <- TIMES * 2
obstimes <- ml.df$datetime
fit_times <- sort(unique(c(seq(min(obstimes), max(obstimes), l=FITTIMES + 2), obstimes)))
## initial ----
mu_proposal_cov <- smover::build.H(kernel = smover::kernel.samp()$kernel.exp,
                                   phi = 3*mean(diff(as.numeric(fit_times))),
                                   warped.times = as.numeric(fit_times), knots = as.numeric(fit_times),
                                   row.normalize = T, l = 1)
mu_proposal_cov_chol <- chol(mu_proposal_cov)
mu_proposal_var <- 1e1
mu_start <- cbind(approx(x = obstimes, y = ml.df[, 1], xout = fit_times)$y,
                  approx(x = obstimes, y = ml.df[, 2], xout = fit_times)$y) +
  rnorm(length(fit_times) * 2, sd = sqrt(mu_proposal_var))
initial <- list(
  "mu" = mu_start, "sigsq_s" = 1e1, "z" = sample(1:0, size = TIMES, replace = T),
  "beta" = rep(0, nlayers(X_scaled)), "sigsq_mu" = 5e6, "tausq" = 5e6,
  "alpha" = 1, "theta" = c(-1/2, 1, rep(0, nlayers(W) - 2)), "g0" = 0, "gmax" = Inf
)
initial_pr_z <- get_fc_z(mu = initial$mu, sigsq_mu = initial$sigsq_mu, tausq = initial$tausq,
                         beta = initial$beta, theta = initial$theta, alpha = initial$alpha,
                         g0 = initial$g0, gmax = initial$gmax,
                         times = fit_times, gradX = gradX, W = W)
initial$z <- sapply(initial_pr_z, function(pr_t) sample(1:0, 1, prob = c(pr_t, 1 - pr_t)))
## fixed ----
fixed <- list(
  "alpha", "gmax" ## used in more general model; held constant here
)
## tuning ----
tuning_init <- list(
  "mu" = mu_proposal_var,
  "beta" = 5e-1, "sigsq_mu" = 1e5, "tausq" = 1e5,
  "alpha" = 1e0, "theta" = 5e-2, "g0" = 1e-2, "gmax" = 1e-1
)
tuning_mu_chol <- TRUE
tuning_cov_init <- list()
## priors ----
mode_s <- 1; var_s <- 1e1
alpha_s <- uniroot(f = function(alpha, mode = mode_s, var = var_s){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e5))$root
beta_s <- mode_s * (alpha_s + 1)
##
mode_mu <- 1e7
var_mu <- 1e16
alpha_mu <- uniroot(f = function(alpha, mode = mode_mu, var = var_mu){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_mu <- mode_mu * (alpha_mu + 1)
##
mode_tau <- 1e7
var_tau <- 1e16
alpha_tau <- uniroot(f = function(alpha, mode = mode_tau, var = var_tau){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_tau <- mode_tau * (alpha_tau + 1)
##
log_priors <- list(
  "sigsq_s" = c(alpha_s, beta_s),
  "beta" = function(beta) sum(dnorm(beta, 0, 1e10, log = T)),
  "sigsq_mu" = function(sigsq_mu) dinvgamma(sigsq_mu, alpha_mu, beta_mu),
  "tausq" = function(tausq) dinvgamma(tausq, alpha_tau, beta_tau),
  "theta" = function(theta) sum(dnorm(theta, 0, 2, log = T)),
  "g0" = function(g0) dnorm(x = g0, mean = 0, sd = 2, log = T)
)
## N_iterations + other ----
data <- list("W" = W, "gradX" = gradX_scaled, "s" = as.matrix(ml.df[, 1:2]),
             "obstimes" = obstimes, "times" = fit_times)
N_iterations <- 1e5; batch_size <- 1e2; block_updates <- TRUE; alpha_blk_reps <- 5
verbose_mcmc <- F; save_tuning_covs <- F; save_gradX <- F; seed <- 0
save_out_rate <- N_iterations + 1; save_out_file <- paste0("../data/mtl/partial_fits/",
                                                           format(Sys.time(), "%Y%m%d_%H%M%S"),
                                                           "_mtl_reDyn_mcmc.RData")
used_iterations <- seq(N_iterations / 2, N_iterations + 1, 10)
## start mcmc algorithm ----
system.time({
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         tuning_cov_init = tuning_cov_init, tuning_mu_chol = tuning_mu_chol,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_gradX = save_gradX,
                         save_out_rate = save_out_rate, save_out_file = save_out_file,
                         used_iterations = used_iterations, seed = seed)
})
## save / load ----
save(out_mcmc, file = paste0("../data/mountain_lion/model_fits/",
                             out_mcmc$time_finished,
                             "_mtl_reDyn_mcmc.RData"))
## plots ----
source("mountain_lion/plot_mtl.R")
