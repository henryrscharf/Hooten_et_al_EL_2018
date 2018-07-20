## libraries ----
library(reDyn); library(smover); library(spatstat)
# ## load fit ----
# load("../data/mountain_lion/model_fits/XXX_mtl_reDyn_mcmc.RData")
## device ----
pdf(file = paste0("../fig/mountain_lion/", out_mcmc$time_finished, "_mtl_reDyn_mcmc.pdf"),
    width = 9)
## plot priors ----
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
layout(matrix(1:(ceiling((length(out_mcmc$log_priors) +
                            nrow(out_mcmc$chains$beta) +
                            nrow(out_mcmc$chains$theta) - 4) / 3) * 3), ncol = 3))
par(mar = c(2, 1, 2, 1))
for(param in names(out_mcmc$log_priors)){
  if(param %in% c("z", "mu", "sigsq_s") || param %in% out_mcmc$fixed) next
  if(param %in% c("sigsq_s", "sigsq_mu", "tausq", "alpha")){
    density <- density(out_mcmc$chains[[param]][1, ], from = 0)
    plot(density, type = "l", main = param)
    lines(density$x, exp(unlist(sapply(density$x, out_mcmc$log_priors[[param]]))), lty = 2)
  } else {
    for(d in 1:nrow(out_mcmc$chains[[param]])){
      density <- density(out_mcmc$chains[[param]][d, ])
      plot(density, type = "l", main = paste(param, d))
      lines(density$x, exp(unlist(sapply(density$x, out_mcmc$log_priors[[param]]))), lty = 2)
    }
  }
}
if(!("sigsq_s" %in% out_mcmc$fixed)){
  density <- density(out_mcmc$chains$sigsq_s[1, ], from = 0)
  plot(density, type = "l", main = "sigsq_s")
  lines(density$x, exp(dinvgamma(density$x, out_mcmc$log_priors$sigsq_s[1],
                                 out_mcmc$log_priors$sigsq_s[2])), lty = 2)
}
## plot traces ----
plot(out_mcmc)
plot(out_mcmc, type = "path")
## plot energetics ----
plot(out_mcmc, type = "energetics", background = X_scaled$elevation)
## correlations ----
chains_df <- sapply(names(out_mcmc$chains)[!(names(out_mcmc$chains) %in%
                                               c('z', 'mu', 'beta', 'theta')) &
                                             !(names(out_mcmc$chains) %in% out_mcmc$fixed)],
                    function(param){
                      return(matrix(c(out_mcmc$chains[[param]]), ncol = 1))
                    })
if(!"beta" %in% out_mcmc$fixed){
  chains_df <- cbind(chains_df, t(out_mcmc$chains$beta))
  colnames(chains_df)[(ncol(chains_df) - nrow(out_mcmc$chains$beta) + 1):ncol(chains_df)] <-
    paste0("beta", 1:nrow(out_mcmc$chains$beta))
}
if(!"theta" %in% out_mcmc$fixed){
  chains_df <- cbind(chains_df, t(out_mcmc$chains$theta))
  colnames(chains_df)[(ncol(chains_df) - nrow(out_mcmc$chains$theta) + 1):ncol(chains_df)] <-
    paste0("theta", 1:nrow(out_mcmc$chains$theta))
}
layout(1)
corrplot::corrplot(cor(chains_df), addCoef.col = T)
## posterior summaries ----
probs <- c(0.025, 0.05, 0.5, 0.95, 0.975)
params_df <- as.data.frame(t(rbind(out_mcmc$chains$beta,
                                   out_mcmc$chains$theta,
                                   out_mcmc$chains$sigsq_s,
                                   out_mcmc$chains$sigsq_mu,
                                   out_mcmc$chains$tausq,
                                   out_mcmc$chains$alpha,
                                   out_mcmc$chains$g0,
                                   out_mcmc$chains$gmax)))
names(params_df) <- c(paste0("X_", names(X_scaled)), paste0("W_", names(out_mcmc$data$W)),
                      "sigsq_s", "sigsq_mu", "tausq", "alpha", "g0", "gmax")
CI_df <- as.data.frame(cbind(
  apply(out_mcmc$chains$beta, 1, quantile, probs = probs),
  apply(out_mcmc$chains$theta, 1, quantile, probs = probs),
  quantile(out_mcmc$chains$sigsq_s, prob = probs),
  quantile(out_mcmc$chains$sigsq_mu, prob = probs),
  quantile(out_mcmc$chains$tausq, prob = probs),
  quantile(out_mcmc$chains$alpha, prob = probs),
  quantile(out_mcmc$chains$g0, prob = probs),
  quantile(out_mcmc$chains$gmax, prob = probs)
))
names(CI_df) <- names(params_df)
t(CI_df)
layout(1)
par(mar = c(3.3, 2, 2, 1))
shown_params <- 1:(sum(nlayers(X_scaled), nlayers(out_mcmc$data$W)))
plot(shown_params, main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, length(shown_params) + 1.5),
     ylim = range(scale(params_df[, shown_params], center = F)))
for(i in shown_params){
  density <- density(scale(params_df[, i], center = F))
  x <- density$y / diff(range(density$y)) * 0.33
  if(i > nrow(out_mcmc$chains$beta)){
    i = i + 1
  }
  polygon(x = c(i - x, i + rev(x)), c(density$x, rev(density$x)),
          col = "gray")
}
axis(2, 0, 0)
abline(h = 0, lty = 2)
text(x = 3, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]), labels = "behavior", font = 2, xpd = T)
text(x = 7, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]), labels = "recharge", font = 2, xpd = T)
names_X <- names(X_scaled)
names_X[names_X == 'dist2kill'] <- "distance to killsite"
names_W <- names(W)
names_W[names_W == 'killsite'] <- "<1km to killsite"
text(x = c(1:nlayers(X_scaled), 1:nlayers(W) + 1 + nlayers(X_scaled)),
     y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 45, xpd = T,
     labels = c(names_X, names_W))
## dev.off() ----
dev.off()

