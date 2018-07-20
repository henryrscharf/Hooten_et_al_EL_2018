## libraries ----
library(reDyn); library(Matrix); library(spatstat); library(maptools)
library(raster); library(ctmcmove)
# ## load ----
# load("../data/simulation/model_fits/XXX_reDyn_mcmc.RData")
## device ----
pdf(file = paste0("../fig/simulation/", out_mcmc$time_finished, "_reDyn_mcmc.pdf"))
## truth ----
truth <- out_mcmc$initial
truth$mu <- path_z$mu
truth$times <- path_z$times
truth$g <- path_z$g
## plot priors ----
## priors
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
layout(matrix(1:length(out_mcmc$log_priors[!(names(out_mcmc$log_priors) %in% c('mu', 'z'))]), ncol = 2))
par(mar = c(2, 1, 2, 1))
for(param in names(out_mcmc$log_priors)){
  if(param %in% c("z", "mu", "sigsq_s") || param %in% out_mcmc$fixed) next
  if(param %in% c("sigsq_s", "sigsq_mu", "tausq", "alpha")){
    density <- density(out_mcmc$chains[[param]][1, ], from = 0)
    plot(density, type = "l", main = param)
    lines(density$x, exp(unlist(sapply(density$x, out_mcmc$log_priors[[param]]))), lty = 2)
    abline(v = truth[param])
  } else {
    for(d in 1:nrow(out_mcmc$chains[[param]])){
      density <- density(out_mcmc$chains[[param]][d, ])
      plot(density, type = "l", main = paste(param, d))
      lines(density$x, exp(unlist(sapply(density$x, out_mcmc$log_priors[[param]]))), lty = 2)
      abline(v = truth[[param]][d])
    }
  }
}
if(!("sigsq_s" %in% out_mcmc$fixed)){
  density <- density(out_mcmc$chains$sigsq_s[1, ], from = 0)
  plot(density, type = "l", main = "sigsq_s")
  lines(density$x, exp(dinvgamma(density$x, out_mcmc$log_priors$sigsq_s[1],
                                 out_mcmc$log_priors$sigsq_s[2])), lty = 2)
  abline(v = truth$sigsq_s)
}
log_priors <- out_mcmc$log_priors
## plot traces ----
plot(out_mcmc, truth = truth)
plot(out_mcmc, type = "path", truth = truth)
## plot energetics ----
plot(out_mcmc, type = "energetics", background = X$layer, truth = truth)
## correlations ----
chains_df <- sapply(names(out_mcmc$chains)[!(names(out_mcmc$chains) %in%
                                               c('z', 'mu', 'beta', 'theta')) &
                                             !(names(out_mcmc$chains) %in% out_mcmc$fixed)],
                    function(param){
                      return(matrix(c(out_mcmc$chains[[param]]), ncol = 1))
                    })
if(!"beta" %in% out_mcmc$fixed){
  chains_df <- cbind(chains_df, t(out_mcmc$chains$beta))
  colnames(chains_df)[ncol(chains_df):(ncol(chains_df) + nrow(out_mcmc$chains$beta) - 1)] <-
    paste0("beta", 1:nrow(out_mcmc$chains$beta))
}
if(!"theta" %in% out_mcmc$fixed){
  chains_df <- cbind(chains_df, t(out_mcmc$chains$theta))
  colnames(chains_df)[ncol(chains_df):(ncol(chains_df) + nrow(out_mcmc$chains$theta) - 1) - 1] <-
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
names(params_df) <- c(paste0("X_", names(X)), paste0("W_", names(out_mcmc$data$W)),
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
shown_params <- 1:(sum(nlayers(X), nlayers(out_mcmc$data$W)))
plot(shown_params, main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, length(shown_params) + 0.5),
     ylim = range(scale(params_df[, shown_params], center = F)))
for(i in shown_params){
  density <- density(scale(params_df[, i], center = F))
  x <- density$y / diff(range(density$y)) * 0.33
  polygon(x = c(i - x, i + rev(x)), c(density$x, rev(density$x)),
          col = "gray")
}
axis(2, 0, 0)
abline(h = 0, lty = 2)
points(1:3, c(truth$beta / sqrt(sum(out_mcmc$chains$beta[1, ]^2) / ncol(out_mcmc$chains$beta)),
              truth$theta / sqrt(sum(out_mcmc$chains$theta[1, ]^2) / ncol(out_mcmc$chains$theta))),
       col = "darkred", pch = 19)
text(x = 1, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]), labels = "behavior", font = 2, xpd = T)
text(x = 2.5, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]), labels = "recharge", font = 2, xpd = T)
text(x = c(1, 2, 3), y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 45, xpd = T,
     labels = c("dist. to polygon", "intercept", "inside polygon"))
## dev.off ----
dev.off()
