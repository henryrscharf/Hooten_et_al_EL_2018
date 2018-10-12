## libraries ----
library(Matrix); library(spatstat); library(maptools)
library(reDyn); library(parallel); library(doParallel)
## simulate covariates ----
# ## GP-type ----
# set.seed(20)
# S <- 4e2; n_covariates <- 1
# calS <- seq(0, 1, l=S)
# H_x <- smover::build.H(kernel = smover::kernel.samp()$kernel.gauss, phi = 1/8,
#                        warped.times = calS, knots = calS, row.normalize = T)
# H_y <- smover::build.H(kernel = smover::kernel.samp()$kernel.gauss, phi = 1/8,
#                        warped.times = calS, knots = calS, row.normalize = T)
# x_list <- vector("list", n_covariates)
# for(cov_i in 1:n_covariates){
#   x_list[[cov_i]] <- raster(ncol = S, nrow = S, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
#   values(x_list[[cov_i]]) <- H_x %*% t(H_y %*% matrix(rnorm(S^2), S, S))
# }
# X <- stack(x_list) ## RasterStack
# X <- -1*X
# W <- (X > quantile(X, 0.9))*2 - 1
# # plot(X)
# # lines(rasterToContour(X - quantile(X, 0.95), levels = 0))
## patch-type ----
# plot(c(0.2, 0.8), c(0.2, 0.8), type = "n")
# locations <- locator(12)
load("../data/simulation/locations.RData")
locations <- lapply(locations, function(x) c(x, x[1]))
W_owin <- owin(xrange = c(-0.25, 1.25), yrange = c(-0.25, 1.25), poly = locations)
X <- stack(raster(distmap(W_owin, eps = 1/4e2)))
names(X) <- "dist2poly"
X_names <- names(X)
projection(X) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
mean_slopes <- unlist(sapply(1:nlayers(X), function(layer){
  mean(values(terrain(X[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
}))
X_grad_scaled <- X / mean_slopes
gradX_scaled <- sapply(1:nlayers(X_grad_scaled), function(layer) rast_grad(X_grad_scaled[[layer]])
                       [c('rast.grad.x', 'rast.grad.y')], simplify = F)
unlist(sapply(1:nlayers(X_grad_scaled), function(layer){
  mean(values(terrain(X_grad_scaled[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
}))
intercept <- raster(X)
values(intercept) <- 1
W <- stack(list("intercept" = intercept,
                "polygon" = (X == 0)))
W_names <- names(W)
# W <- stack(list("not_polygon" = X != 0,
#                 "polygon" = X == 0))
## simulate path ----
seed <- 3
gradX <- sapply(1:nlayers(X), function(layer) ctmcmove::rast.grad(X[[layer]]), simplify = F)
method <- "position"; innovation = "BM"
if(method == "position"){
  sigsq_mu <- 2e-2; tausq <- 3e-2; alpha <- 1
  beta <- maxValue(gradX[[1]]$rast.grad.x)^(-1) * c(-1/2)
  theta <- c(-1, 4)
  # theta <- c(1, 2)
  sigsq_s <- 1e-5
  if(innovation == "IBM"){
    beta <- beta * 1
    tausq <- tausq * 1e-1
    sigsq_mu <- sigsq_mu * 1e-1
  }
  gamma <- 0
}
if(method == "velocity"){
  sigsq_mu <- 1e-2; tausq <- 1e-2; alpha <- 6
  beta <- maxValue(gradX[[1]]$rast.grad.x)^(-1) * c(-5e-1)
  theta <- 3 * c(-1, 2)
  sigsq_s <- 1e-5
  gamma <- 0.01
}
mu0 <- c(1, 1); g0 <- -1; gmax <- Inf
TIMES <- 3e3; TIME <- 20
set.seed(seed)
times <- sort(c(0, TIME, TIME * runif(TIMES - 2)))
OBSTIMES <- min(2^9 + 2, TIMES)
obstimes <- sort(c(0, TIME, sample(times[-c(1, TIMES)], OBSTIMES - 2)))
path_z <- sim_path_z(sigsq_mu = sigsq_mu, tausq = tausq, beta = beta,
                     theta = theta, alpha = alpha, sigsq_s = sigsq_s,
                     gamma = gamma, mu0 = mu0, g0 = g0, gmax = gmax,
                     times = times, obstimes = obstimes, gradX = gradX, W = W,
                     method = method, innovation = innovation, seed = seed)
## orthogonalize W ----
W_ortho <- W
W_path <- extract(x = W, y = matrix(path_z$mu, ncol = 2))
W_tilde <- apply(W_path * c(0, diff(times)), 2, cumsum)
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
# W_path <- extract(x = W_ortho, y = matrix(path_z$mu, ncol = 2))
# W_tilde <- apply(W_path * c(0, diff(times)), 2, cumsum)
# matplot(times, W_path, type = "l")
# corrplot::corrplot(t(W_path) %*% W_path, is.corr = F)
# matplot(times, W_tilde, type = "l")
# corrplot::corrplot(t(W_tilde) %*% W_tilde, is.corr = F, addCoef.col = T)
# ## orthogonalize W (old) ----
# w_polygon <- extract(W$polygon, path_z$mu[-TIMES, ])
# w_til_polygon <- cumsum(w_polygon * diff(times)) ## accumulated w_polygon = \int w(mu(t))dt
# w_til_1 <- times[-1] ## accumulated w_int = \int 1 dt = t
# g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon) %*% path_z$theta)
# P_til_1 <- w_til_1 %*% (w_til_1 %*% w_til_1)^(-1) %*% w_til_1
# w_til_polygon_orth <- w_til_polygon - P_til_1 %*% w_til_polygon
# unique(round(diff(w_til_polygon_orth) / diff(times[-1]), 3))
# W$polygon <- W$polygon + min(diff(w_til_polygon_orth) / diff(times[-1]))
# # W$polygon <- scale(W$polygon)
# theta[1] <- path_z$theta[1] - min(diff(w_til_polygon_orth) / diff(times[-1])) * path_z$theta[2]
# # theta[2] <- path_z$theta[2] / sqrt(sum(values(W$polygon)^2) /
# #                                      (prod(ncol(W$polygon), nrow(W$polygon)) - 1))
# g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon_orth) %*% theta)
# w_polygon_test <- extract(W$polygon, path_z$mu[-TIMES, ])
# w_til_polygon_test <- cumsum(w_polygon_test * diff(times)) ## accumulated w_polygon = \int w(mu(t))dt
# g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon_test) %*% theta)
## set up timing ----
set.seed(0)
FITTIMES <- TIMES * 0.5
fit_times <- sort(c(obstimes, sample(times[!(times %in% obstimes)],
                                     size = FITTIMES - OBSTIMES)))
## initial ----
mu_proposal_var <- 1e-4
mu_start <- path_z$mu[times %in% fit_times, ]
initial <- list(
  "mu" = mu_start, "sigsq_s" = sigsq_s, "z" = path_z$z[times %in% fit_times],
  "beta" = beta * mean_slopes, "sigsq_mu" = sigsq_mu,  "tausq" = tausq,
  "alpha" = alpha, "theta" = c(solve(W_tilde_proj_mat) %*% theta), "g0" = g0, "gmax" = gmax
)
## fixed ----
fixed <- c(
  # "z",
  # "theta", "g0",
  # "mu",
  # "sigsq_s",
  # "beta", "sigsq_mu", "tausq",
  "alpha", "gmax"
)
initial$z <- rep(0, length(fit_times))
## tuning ----
tuning_init <- list(
  "mu" = mu_proposal_var,
  "beta" = 1, "sigsq_mu" = (4e-2)^2, "tausq" = (4e-2)^2,
  "alpha" = 1, "theta" = 1e-1, "gmax" = 1e-2, "g0" = 1
)
tuning_cov_init <- list()
## priors ----
mode_s <- sigsq_s; var_s <- 1e-2
alpha_s <- uniroot(f = function(alpha, mode = mode_s, var = var_s){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e1))$root
beta_s <- mode_s * (alpha_s + 1)
##
mode_mu <- sigsq_mu; var_mu <- 1
alpha_mu <- uniroot(f = function(alpha, mode = mode_mu, var = var_mu){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e1))$root
beta_mu <- mode_mu * (alpha_mu + 1)
##
mode_tau <- tausq; var_tau <- 1
alpha_tau <- uniroot(f = function(alpha, mode = mode_tau, var = var_tau){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e1))$root
beta_tau <- mode_tau * (alpha_tau + 1)
##
log_priors <- list(
  "sigsq_s" = c(alpha_s, beta_s),
  "beta" = function(beta) sum(dnorm(beta, 0, 50, log = T)),
  "sigsq_mu" = function(sigsq_mu) dinvgamma(sigsq_mu, alpha_mu, beta_mu),
  "tausq" = function(tausq) dinvgamma(tausq, alpha_tau, beta_tau),
  "theta" = function(theta) sum(dnorm(theta, 0, sqrt(diff(range(path_z$times)) * 50), log = T)),
  "g0" = function(g0) dnorm(x = g0, mean = 0, sd = 1, log = T),
  "alpha" = function(x) 0,
  "gmax" = function(x) 0
)
## k-fold indicies ----
set.seed(1)
nfolds <- 8
pts_per_fold <- ceiling((length(obstimes) - 2) / nfolds)
left_over <- pts_per_fold * nfolds - (length(obstimes) - 2)
holdout_indices <- apply(matrix(c(sample(2:(length(obstimes) - 1), length(obstimes) - 2),
                                  rep(0, left_over)), ncol = nfolds), 2, sort)
## N_iterations + other ----
ncores <- nfolds
data <- list("W" = W_ortho, "gradX" = gradX_scaled, "s" = path_z$s,
             "obstimes" = obstimes, "times" = fit_times)
N_iterations <- 1e5; batch_size <- 1e2; block_updates <- T; alpha_blk_reps <- 5
verbose_mcmc <- F; save_tuning_covs <- F; save_gradX <- F; save_W <- F; RBZ <- F; seed <- 0
save_out_rate <- N_iterations + 1
save_out_file <- paste0("../data/simulation/partial_fits/",
                        format(Sys.time(), "%Y%m%d_%H%M%S"),
                        "_reDyn_mcmc.RData")
used_iterations <- seq(1, N_iterations + 1, 10)
## ----
## k-fold fits ----
## register parallel backend ----
doParallel::registerDoParallel(ncores)
## start mcmc algorithms for FULL model ----
init_time <- Sys.time()
out_mcmcs_FULL <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
  data$s <- data$s[-holdout, ]
  data$obstimes <- data$obstimes[-holdout]
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
                         save_out_rate = save_out_rate, save_out_file = save_out_file,
                         used_iterations = used_iterations, seed = seed)
}
out_mcmcs_FULL$time_started <- format(init_time, "%Y%m%d_%H%M%S")
out_mcmcs_FULL$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
class(out_mcmcs_FULL) <- "reDyn_mcmcs"
Sys.time() - init_time
## save FULL ----
save(out_mcmcs_FULL, X_names, W_names,
     W_ortho, X_grad_scaled,
     log_priors, W_tilde_proj_mat, mean_slopes,
     holdout_indices,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_FULL$time_finished,
                   "_simulation_reDyn_mcmc_FULL.RData"))
## score FULL ----
system.time({
  lds_FULL <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
    get_score(reDyn_mcmc = out_mcmcs_FULL[[i]],
              holdout = data$s[holdout_indices[, i], ],
              holdout_times = data$obstimes[holdout_indices[, i]])
  }
})
lds_FULL_mat <- matrix(c(unlist(lds_FULL), rep(NA, left_over)), ncol = 8)
mean(lds_FULL_mat, na.rm = T)
## save FULL score ----
save(lds_FULL_mat,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_FULL$time_finished,
                   "_simulation_reDyn_score_FULL.RData"))
## ----
## prep SUB model (M0) ----
## fixed ----
fixed <- list(
  "z", "g0", "theta",
  # "mu",
  # "sigsq_s",
  "alpha", "gmax"
)
## initial ----
initial$theta <- rep(0, length(theta))
initial$z <- rep(0, FITTIMES)
## start mcmc algorithms for M0 model ----
init_time <- Sys.time()
out_mcmcs_M0 <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
  data$s <- data$s[-holdout, ]
  data$obstimes <- data$obstimes[-holdout]
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
                         save_out_rate = save_out_rate, save_out_file = save_out_file,
                         used_iterations = used_iterations, seed = seed)
}
out_mcmcs_M0$time_started <- format(init_time, "%Y%m%d_%H%M%S")
out_mcmcs_M0$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
class(out_mcmcs_M0) <- "reDyn_mcmcs"
Sys.time() - init_time
## save M0 ----
save(out_mcmcs_M0, log_priors, holdout_indices,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_M0$time_finished,
                   "_simulation_reDyn_mcmc_M0.RData"))
## score M0 ----
system.time({
  lds_M0 <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
    get_score(reDyn_mcmc = out_mcmcs_M0[[i]],
              holdout = data$s[holdout_indices[, i], ],
              holdout_times = data$obstimes[holdout_indices[, i]])
  }
})
lds_M0_mat <- matrix(c(unlist(lds_M0), rep(NA, left_over)), ncol = 8)
mean(lds_M0_mat, na.rm = T)
## save M0 score ----
save(lds_M0_mat,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_M0$time_finished,
                   "_simulation_reDyn_score_M0.RData"))
## ----
## prep SUB model (M1) ----
## fixed ----
fixed <- list(
  "z", "g0", "theta",
  # "mu",
  # "sigsq_s",
  "alpha", "gmax"
)
initial$z <- rep(1, FITTIMES)
## start mcmc algorithms for M1 model ----
init_time <- Sys.time()
out_mcmcs_M1 <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
  data$s <- data$s[-holdout, ]
  data$obstimes <- data$obstimes[-holdout]
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_gradX = save_gradX, save_W = save_W,
                         save_out_rate = save_out_rate, save_out_file = save_out_file,
                         used_iterations = used_iterations, seed = seed)
}
out_mcmcs_M1$time_started <- format(init_time, "%Y%m%d_%H%M%S")
out_mcmcs_M1$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
class(out_mcmcs_M1) <- "reDyn_mcmcs"
Sys.time() - init_time
## save M1 ----
save(out_mcmcs_M1, X_names, X_grad_scaled,
     log_priors, mean_slopes, holdout_indices,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_M1$time_finished,
                   "_simulation_reDyn_mcmc_M1.RData"))
## score M1 ----
system.time({
  lds_M1 <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
    get_score(reDyn_mcmc = out_mcmcs_M1[[i]],
              holdout = data$s[holdout_indices[, i], ],
              holdout_times = data$obstimes[holdout_indices[, i]])
  }
})
lds_M1_mat <- matrix(c(unlist(lds_M1), rep(NA, left_over)), ncol = 8)
mean(lds_M1_mat, na.rm = T)
## save M1 score ----
save(lds_M1_mat,
     file = paste0("../data/simulation/multi_model_fits/",
                   out_mcmcs_M1$time_finished,
                   "_simulation_reDyn_score_M1.RData"))
## ----
# ## start mcmc algorithms for SUB model ----
# fixed <- c(
#   "z",
#   # "beta", "sigsq_mu",
#   "theta", "g0",
#   "alpha", "gmax"
# )
# initial$z <- rep(0, FITTIMES)
# init_time <- Sys.time()
# out_mcmcs_SUB <- foreach::foreach(holdout = iterators::iter(holdout_indices, by = 'col')) %dopar% {
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
# out_mcmcs_SUB$time_started <- format(init_time, "%Y%m%d_%H%M%S")
# out_mcmcs_SUB$time_finished <- format(Sys.time(), "%Y%m%d_%H%M%S")
# class(out_mcmcs_SUB) <- "reDyn_mcmcs"
# Sys.time() - init_time
# ## save SUB ----
# save(out_mcmcs_SUB, X_names, W_names,
#      W_ortho, X_grad_scaled,
#      log_priors, W_tilde_proj_mat, mean_slopes,
#      holdout_indices,
#      file = paste0("../data/buffalo/model_fits/",
#                    out_mcmcs_SUB$time_finished,
#                    "_simulation_reDyn_mcmc_SUB.RData"))
# ## score SUB ----
# system.time({
#   lds_SUB <- foreach(i = 1:ncol(holdout_indices)) %dopar% {
#     get_score(reDyn_mcmc = out_mcmcs_SUB[[i]],
#               holdout = data$s[holdout_indices[, i], ],
#               holdout_times = data$obstimes[holdout_indices[, i]])
#   }
# })
# lds_SUB_mat <- matrix(c(unlist(lds_SUB), rep(NA, left_over)), ncol = 8)
# mean(lds_SUB_mat, na.rm = T)
# ## save SUB score ----
# save(lds_SUB_mat,
#      file = paste0("../data/simulation/multi_model_fits/",
#                    out_mcmcs_SUB$time_finished,
#                    "_simulation_reDyn_score_SUB.RData"))
## ----
# ## plots ----
# out_mcmc <- out_mcmcs[[1]]
# holdout_index <- holdout_indices[, 1]
# source("simulate/plot_simulation.R")
