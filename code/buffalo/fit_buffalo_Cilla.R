## libraries ----
library(smover)
library(reDyn)
library(doParallel)
## load buffalo data ----
load("../data/buffalo/buffalo_Cilla.RData")
# buffalo_proj <- buffalo_proj[buffalo_proj$POSIX > "2005-09-14" &
#                                buffalo_proj$POSIX < "2005-10-15", ]
buffalo_proj <- buffalo_proj[buffalo_proj$POSIX > "2005-09-30" &
                               buffalo_proj$POSIX < "2005-10-15", ]
obstimes <- as.numeric(buffalo_proj$POSIX) / 3600
## rasters + sabie river ----
ELE_rast <- raster("../data/buffalo/buffaloELE_Cilla.grd")
load("../data/buffalo/dist2sabie_ext.RData")
## covariates ----
## X
X <- stack(list(
  # "elevation" = ELE_rast,
  # "slope" = terrain(ELE_rast, opt = "slope"),
  "dist2sabie" = dist2sabie
  # "sin_aspect" = sin(terrain(ELE_rast, opt = "aspect")),
  # "cos_aspect" = cos(terrain(ELE_rast, opt = "aspect"))
))
X_names <- names(X)
# X$slope[is.na(X$slope)] <- mean(values(X$slope), na.rm = T)
# X$sin_aspect[is.na(X$sin_aspect)] <- mean(values(X$sin_aspect), na.rm = T)
# X$cos_aspect[is.na(X$cos_aspect)] <- mean(values(X$cos_aspect), na.rm = T)
X_scaled <- scale(X)
## standardize X based on slopes of gradients
mean_slopes <- unlist(sapply(1:nlayers(X), function(layer){
  mean(values(terrain(X[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
}))
X_grad_scaled <- X / mean_slopes
unlist(sapply(1:nlayers(X_grad_scaled), function(layer){
  mean(values(terrain(X_grad_scaled[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
}))
gradX_scaled <- sapply(1:nlayers(X_grad_scaled), function(layer) rast_grad(X_grad_scaled[[layer]])
                       [c('rast.grad.x', 'rast.grad.y')], simplify = F)
## W
intercept <- raster(X_scaled)
values(intercept) <- 1
W <- stack(list(
  "intercept" = intercept,
  # "elevation" = X_scaled$elevation,
  # "slope" = X_scaled$slope,
  "near_sabie" = dist2sabie < 0.5e3
  # "sin_aspect" = X_scaled$sin_aspect,
  # "cos_aspect" = X_scaled$cos_aspect
  # "elevation_slope" = X_scaled$elevation * X_scaled$slope
))
W_names <- names(W)
## set up timing ----
TIMES <- nrow(buffalo_proj)
n_pts_btw <- 3
fit_times <- sort(c(rep(1, n_pts_btw) %x% obstimes[-TIMES] + (1:n_pts_btw) %x% diff(obstimes) / (n_pts_btw + 1), obstimes))
# fit_times <- obstimes
# get_intermediates <- function(x){
#   (x[-1] + x[-length(x)]) / 2
# }
# fit_times <- obstimes
# for(i in 1:2){
#   fit_times <- sort(c(fit_times, get_intermediates(fit_times)))
# }
FITTIMES <- length(fit_times)
## initial mu ----
mu_proposal_var <- 1
mu_start <- cbind(approx(x = obstimes, y = buffalo_proj@coords[, 1], xout = fit_times)$y,
                  approx(x = obstimes, y = buffalo_proj@coords[, 2], xout = fit_times)$y) +
  rnorm(length(fit_times) * 2, sd = sqrt(mu_proposal_var))
print("distribution of step lengths in meters:")
print(quantile(sqrt(rowSums((apply(mu_start, 2, diff))^2))))
## orthogonalize W based on initial mu ----
W_ortho <- W
W_path <- extract(x = W, y = matrix(mu_start, ncol = 2))
W_tilde <- apply(W_path * c(0, diff(fit_times)), 2, cumsum)
W_tilde_svd <- svd(W_tilde)
W_tilde_proj_mat <- W_tilde_svd$v %*% diag(W_tilde_svd$d^(-1))
W_mat <- as.matrix(W)
W_mat_proj <- W_mat %*% W_tilde_proj_mat
for(layer in 1:ncol(W_mat)){
  values(W_ortho[[layer]]) <- W_mat_proj[, layer]
  names(W_ortho[[layer]]) <- paste0("svd", layer)
}
# ## check orthogonality of W_tilde ----
# layout(matrix(1:4, 2, 2))
# W_path <- extract(x = W_ortho, y = matrix(mu_start, ncol = 2))
# W_tilde <- apply(W_path * c(0, diff(fit_times)), 2, cumsum)
# matplot(fit_times, W_path, type = "l")
# corrplot::corrplot(t(W_path) %*% W_path, is.corr = F)
# matplot(fit_times, W_tilde, type = "l")
# corrplot::corrplot(t(W_tilde) %*% W_tilde, is.corr = F, addCoef.col = T)
## remove large leftover objects ----
rm(X_grad_scaled, X_scaled, W_mat_proj, W_mat, W, dist2sabie, intercept)
## other initial values ----
initial <- list(
  "mu" = mu_start, "sigsq_s" = mu_proposal_var, "z" = sample(1:0, size = TIMES, replace = T),
  "beta" = rep(0, nlayers(X)), "sigsq_mu" = 3e2^2, "tausq" = 3e2^2,
  "alpha" = 1, "theta" = rep(0, nlayers(W_ortho)), "g0" = 0, "gmax" = Inf
)
initial_pr_z <- get_fc_z(mu = initial$mu, sigsq_mu = initial$sigsq_mu, tausq = initial$tausq,
                         beta = initial$beta, theta = initial$theta, alpha = initial$alpha,
                         g0 = initial$g0, gmax = initial$gmax,
                         times = fit_times, gradX = gradX_scaled, W = W_ortho)
initial$z <- sapply(initial_pr_z, function(pr_t) sample(1:0, 1, prob = c(pr_t, 1 - pr_t)))
## fixed ----
fixed <- c(
  # "mu",
  # "sigsq_s",
  "alpha", "gmax"
)
## tuning ----
tuning_init <- list(
  "mu" = mu_proposal_var,
  "beta" = 1e2, "sigsq_mu" = 1e2^2, "tausq" = 1e2^2,
  "theta" = 1, "g0" = 1
)
tuning_cov_init <- list()
## priors ----
mode_s <- initial$sigsq_s; var_s <- initial$sigsq_s * 5^2
alpha_s <- uniroot(f = function(alpha, mode = mode_s, var = var_s){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e5))$root
beta_s <- mode_s * (alpha_s + 1)
##
mode_mu <- initial$sigsq_mu
var_mu <- mode_mu^2
alpha_mu <- uniroot(f = function(alpha, mode = mode_mu, var = var_mu){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_mu <- mode_mu * (alpha_mu + 1)
##
mode_tau <- initial$tausq
var_tau <- mode_tau^2
alpha_tau <- uniroot(f = function(alpha, mode = mode_tau, var = var_tau){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_tau <- mode_tau * (alpha_tau + 1)
##
log_priors <- list(
  "sigsq_s" = c(alpha_s, beta_s),
  "beta" = function(beta) sum(dnorm(x = beta, mean = 0, sd = 1e2, log = T)),
  "sigsq_mu" = function(sigsq_mu) dinvgamma(sigsq_mu, alpha_mu, beta_mu),
  "tausq" = function(tausq) dinvgamma(tausq, alpha_tau, beta_tau),
  "theta" = function(theta) sum(dnorm(x = theta, mean = 0, sd = sqrt(diff(range(obstimes))), log = T)),
  "g0" = function(g0) dnorm(x = g0, mean = 0, sd = 1, log = T),
  "alpha" = function(x) 0,
  "gmax" = function(x) 0
)
# ## check variance priors ----
# layout(matrix(1:3, 1, 3))
# plot(sqrt(seq(mode_s / 10, mode_s * 10, l=2e2)),
#      exp(dinvgamma(seq(mode_s / 10, mode_s * 10, l=2e2), alpha_s, beta_s)),
#      xlab = expression(sigma[s]), ylab = "density", type = "l")
# plot(sqrt(seq(mode_mu / sqrt(10), mode_mu * 10, l=2e2)),
#      exp(dinvgamma(seq(mode_mu / 3, mode_mu * 3, l=2e2), alpha_mu, beta_mu)),
#      xlab = expression(sigma[mu]), ylab = "density", type = "l")
# plot(sqrt(seq(mode_tau / sqrt(10), mode_tau * 10, l=2e2)),
#      exp(dinvgamma(seq(mode_tau / 3, mode_tau * 3, l=2e2), alpha_tau, beta_tau)),
#      xlab = expression(sigma[tau]), ylab = "density", type = "l")
## k-fold indicies ----
set.seed(1)
nfolds <- 8
pts_per_fold <- ceiling((length(obstimes) - 2) / nfolds)
left_over <- pts_per_fold * nfolds - (length(obstimes) - 2)
holdout_indices <- apply(matrix(c(sample(2:(length(obstimes) - 1), length(obstimes) - 2),
                                  rep(0, left_over)), ncol = nfolds), 2, sort)
## N_iterations + other ----
ncores <- nfolds
data <- list("W" = W_ortho, "gradX" = gradX_scaled, "s" = as.matrix(buffalo_proj@coords),
             "obstimes" = obstimes, "times" = fit_times)
N_iterations <- 1e5; batch_size <- 2e2; block_updates <- TRUE; alpha_blk_reps <- 5
verbose_mcmc <- F; save_tuning_covs <- F; save_gradX <- F; save_W <- F; seed <- 0
save_out_rate <- N_iterations + 1
# save_out_rate <- 100
save_out_file <- paste0("../data/buffalo/partial_fits/",
                        format(Sys.time(), "%Y%m%d_%H%M%S"),
                        "_buffalo_Cilla_reDyn_mcmc_partial.RData")
used_iterations <- seq(N_iterations / 2, N_iterations + 1, 10)
## ----
## k-fold fits ----
## register parallel backend ----
doParallel::registerDoParallel(ncores)
## ----
# ## start mcmc algorithms for FULL model ----
# init_time <- Sys.time()
# out_mcmcs_FULL <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
#   data$s <- data$s[-holdout, ]
#   data$obstimes <- data$obstimes[-holdout]
#   out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
#                          log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
#                          batch_size = batch_size, block_updates = block_updates,
#                          alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
#                          save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
#                          save_out_rate = save_out_rate, save_out_file = save_out_file,
#                          used_iterations = used_iterations, seed = seed)
# }
# out_mcmcs_FULL$time_started <- format(init_time, "%Y%m%d_%H%M%S")
# out_mcmcs_FULL$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
# class(out_mcmcs_FULL) <- "reDyn_mcmcs"
# Sys.time() - init_time
# ## save FULL ----
# save(out_mcmcs_FULL, X_names, W_names,
#      W_ortho, X,
#      log_priors, W_tilde_proj_mat, mean_slopes,
#      holdout_indices,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_FULL$time_finished,
#                    "_buffalo_Cilla_reDyn_mcmc_FULL.RData"))
# ## score FULL ----
# system.time({
#   lds_FULL <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
#     get_score(reDyn_mcmc = out_mcmcs_FULL[[i]],
#               holdout = data$s[holdout_indices[, i], ],
#               holdout_times = data$obstimes[holdout_indices[, i]])
#   }
# })
# lds_FULL_mat <- matrix(c(unlist(lds_FULL), rep(NA, left_over)), ncol = 8)
# mean(lds_FULL_mat, na.rm = T)
# ## save FULL ----
# save(lds_FULL_mat,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_FULL$time_finished,
#                    "_buffalo_Cilla_reDyn_scores_FULL.RData"))
# ## ----
# ## prep SUB model (M0) ----
# ## fixed ----
# fixed <- list(
#   "z", "g0", "theta",
#   # "mu",
#   # "sigsq_s",
#   "alpha", "gmax"
# )
# initial$z <- rep(0, FITTIMES)
# ## start mcmc algorithms for SUB model ----
# init_time <- Sys.time()
# out_mcmcs_M0 <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
#   data$s <- data$s[-holdout, ]
#   data$obstimes <- data$obstimes[-holdout]
#   out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
#                          log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
#                          batch_size = batch_size, block_updates = block_updates,
#                          alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
#                          save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
#                          save_out_rate = save_out_rate, save_out_file = save_out_file,
#                          used_iterations = used_iterations, seed = seed)
# }
# out_mcmcs_M0$time_started <- format(init_time, "%Y%m%d_%H%M%S")
# out_mcmcs_M0$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
# class(out_mcmcs_M0) <- "reDyn_mcmcs"
# Sys.time() - init_time
# ## save M0 ----
# save(out_mcmcs_M0, log_priors,
#      holdout_indices,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_M0$time_finished,
#                    "_buffalo_Cilla_reDyn_mcmc_M0.RData"))
# ## score M0 ----
# system.time({
#   lds_M0 <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
#     get_score(reDyn_mcmc = out_mcmcs_M0[[i]],
#               holdout = data$s[holdout_indices[, i], ],
#               holdout_times = data$obstimes[holdout_indices[, i]])
#   }
# })
# lds_M0_mat <- matrix(c(unlist(lds_M0), rep(NA, left_over)), ncol = 8)
# mean(lds_M0_mat, na.rm = T)
# ## save M0 ----
# save(lds_M0_mat,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_M0$time_finished,
#                    "_buffalo_Cilla_reDyn_scores_M0.RData"))
# ## ----
# ## prep SUB model (M1) ----
# ## fixed ----
# fixed <- list(
#   "z", "g0", "theta",
#   # "mu",
#   # "sigsq_s",
#   "alpha", "gmax"
# )
# initial$z <- rep(1, FITTIMES)
# ## start mcmc algorithms for SUB model ----
# init_time <- Sys.time()
# out_mcmcs_M1 <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
#   data$s <- data$s[-holdout, ]
#   data$obstimes <- data$obstimes[-holdout]
#   out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
#                          log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
#                          batch_size = batch_size, block_updates = block_updates,
#                          alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
#                          save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
#                          save_out_rate = save_out_rate, save_out_file = save_out_file,
#                          used_iterations = used_iterations, seed = seed)
# }
# out_mcmcs_M1$time_started <- format(init_time, "%Y%m%d_%H%M%S")
# out_mcmcs_M1$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
# class(out_mcmcs_M1) <- "reDyn_mcmcs"
# Sys.time() - init_time
# ## save M1 ----
# save(out_mcmcs_M1, X_names, X,
#      log_priors, mean_slopes, holdout_indices,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_M1$time_finished,
#                    "_buffalo_Cilla_reDyn_mcmc_M1.RData"))
# ## score M1 ----
# system.time({
#   lds_M1 <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
#     get_score(reDyn_mcmc = out_mcmcs_M1[[i]],
#               holdout = data$s[holdout_indices[, i], ],
#               holdout_times = data$obstimes[holdout_indices[, i]])
#   }
# })
# lds_M1_mat <- matrix(c(unlist(lds_M1), rep(NA, left_over)), ncol = 8)
# mean(lds_M1_mat, na.rm = T)
# ## save M1 ----
# save(lds_M1_mat,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_M1$time_finished,
#                    "_buffalo_Cilla_reDyn_scores_M1.RData"))
# ## ----
## single fit ----
## ----
## start mcmc algorithm ----
system.time({
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_gradX = save_gradX,
                         save_out_rate = save_out_rate, save_out_file = save_out_file,
                         used_iterations = used_iterations, seed = seed)
})
## save / load ----
save(out_mcmc, X_names, W_names,
     W_ortho, X,
     log_priors, W_tilde_proj_mat, mean_slopes,
     file = paste0("../data/buffalo/model_fits/",
                   out_mcmc$time_finished,
                   "_buffalo_Cilla_reDyn_mcmc.RData"))
## plots ----
source("buffalo/plot_buffalo.R")
