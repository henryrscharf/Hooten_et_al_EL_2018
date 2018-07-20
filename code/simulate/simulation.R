## libraries ----
library(Matrix); library(spatstat); library(maptools)
library(reDyn)
## simulate covariates ----
## patch-type ----
load("../data/simulation/locations.RData")
locations <- lapply(locations, function(x) c(x, x[1]))
W_owin <- owin(xrange = c(-0.25, 1.25), yrange = c(-0.25, 1.25), poly = locations)
X <- stack(raster(distmap(W_owin, eps = 1/4e2)))
projection(X) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
intercept <- raster(X)
values(intercept) <- 1
W <- stack(list("intercept" = intercept,
                "polygon" = (X == 0)))
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
OBSTIMES <- min(5e2, TIMES)
obstimes <- sort(c(0, TIME, sample(times[-c(1, TIMES)], OBSTIMES - 2)))
path_z <- sim_path_z(sigsq_mu = sigsq_mu, tausq = tausq, beta = beta,
                     theta = theta, alpha = alpha, sigsq_s = sigsq_s,
                     gamma = gamma, mu0 = mu0, g0 = g0, gmax = gmax,
                     times = times, obstimes = obstimes, gradX = gradX, W = W,
                     method = method, innovation = innovation, seed = seed)
## orthogonalize W ----
w_polygon <- extract(W$polygon, path_z$mu[-TIMES, ])
w_til_polygon <- cumsum(w_polygon * diff(times)) ## accumulated w_polygon = \int w(mu(t))dt
w_til_1 <- times[-1] ## accumulated w_int = \int 1 dt = t
g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon) %*% path_z$theta)
P_til_1 <- w_til_1 %*% (w_til_1 %*% w_til_1)^(-1) %*% w_til_1
w_til_polygon_orth <- w_til_polygon - P_til_1 %*% w_til_polygon
unique(round(diff(w_til_polygon_orth) / diff(times[-1]), 3))
W$polygon <- W$polygon + min(diff(w_til_polygon_orth) / diff(times[-1]))
theta[1] <- path_z$theta[1] - min(diff(w_til_polygon_orth) / diff(times[-1])) * path_z$theta[2]
g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon_orth) %*% theta)
w_polygon_test <- extract(W$polygon, path_z$mu[-TIMES, ])
w_til_polygon_test <- cumsum(w_polygon_test * diff(times)) ## accumulated w_polygon = \int w(mu(t))dt
g_test <- path_z$g0 + c(0, cbind(w_til_1, w_til_polygon_test) %*% theta)
## set up model fit to data ----
set.seed(0)
FITTIMES <- TIMES * 0.5
fit_times <- sort(c(obstimes, sample(times[!(times %in% obstimes)],
                                     size = FITTIMES - OBSTIMES)))
## initial ----
mu_proposal_var <- 1e-4
mu_start <- path_z$mu[times %in% fit_times, ]
initial <- list(
  "mu" = mu_start, "sigsq_s" = sigsq_s, "z" = path_z$z[times %in% fit_times],
  "beta" = beta, "sigsq_mu" = sigsq_mu,  "tausq" = tausq,
  "alpha" = alpha, "theta" = theta, "g0" = g0, "gmax" = gmax
)
## fixed ----
fixed <- c(
  "alpha", "gmax"
)
## tuning ----
tuning_init <- list(
  "mu" = mu_proposal_var,
  "beta" = 1e4, "sigsq_mu" = (4e-2)^2, "tausq" = (4e-2)^2,
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
  "beta" = function(beta) sum(dnorm(beta, 0, 1e7, log = T)),
  "sigsq_mu" = function(sigsq_mu) dinvgamma(sigsq_mu, alpha_mu, beta_mu),
  "tausq" = function(tausq) dinvgamma(tausq, alpha_tau, beta_tau),
  "theta" = function(theta) sum(dnorm(theta, 0, 2, log = T)),
  "g0" = function(g0) dnorm(x = g0, mean = 0, sd = 2, log = T)
)
## N_iterations + other ----
data <- list("W" = W, "gradX" = gradX, "s" = path_z$s,
             "obstimes" = obstimes, "times" = fit_times)
N_iterations <- 1e5; batch_size <- 1e2; block_updates <- T; alpha_blk_reps <- 5
verbose_mcmc <- F; save_tuning_covs <- F; save_gradX <- T; RBZ <- F; seed <- 0
save_out_rate <- N_iterations + 1
save_out_file <- paste0("../data/simulation/partial_fits/",
                        format(Sys.time(), "%Y%m%d_%H%M%S"),
                        "_reDyn_mcmc.RData")
used_iterations <- seq(1, N_iterations + 1, 10)
## start mcmc algorithm ----
system.time({
  out_mcmc <- reDyn_mcmc(N_iterations = N_iterations, data = data, initial = initial,
                         log_priors = log_priors, fixed = fixed, tuning_init = tuning_init,
                         tuning_cov_init = tuning_cov_init, tuning_mu_chol = tuning_mu_chol,
                         batch_size = batch_size, block_updates = block_updates,
                         alpha_blk_reps = alpha_blk_reps, verbose_mcmc = verbose_mcmc,
                         save_tuning_covs = save_tuning_covs, save_out_rate = save_out_rate,
                         save_out_file = save_out_file, used_iterations = used_iterations, seed = seed,
                         save_gradX = save_gradX, RBZ = RBZ)
})
## save ----
save(out_mcmc, path_z, X, file = paste0("../data/simulation/model_fits/",
                                        out_mcmc$time_finished,
                                        "_reDyn_mcmc.RData"))
## plots ----
source("simulate/plot_simulation.R")
