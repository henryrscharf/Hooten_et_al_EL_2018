## libraries ----
library(reDyn); library(Matrix); library(spatstat); library(maptools)
library(raster); library(ctmcmove)
## device ----
pdf(file = paste0("../fig/simulation/", out_mcmc$time_finished, "_sim_reDyn_mcmc.pdf"))
## truth ----
truth <- out_mcmc$initial
truth$mu <- path_z$mu
truth$times <- path_z$times
truth$g <- path_z$g
## plot priors ----
## priors
mode_s <- path_z$sigsq_s; var_s <- 1e-2
alpha_s <- uniroot(f = function(alpha, mode = mode_s, var = var_s){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e1))$root
beta_s <- mode_s * (alpha_s + 1)
##
mode_mu <- path_z$sigsq_mu; var_mu <- 1
alpha_mu <- uniroot(f = function(alpha, mode = mode_mu, var = var_mu){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e1))$root
beta_mu <- mode_mu * (alpha_mu + 1)
##
mode_tau <- path_z$tausq; var_tau <- 1
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
plot(out_mcmc, type = "energetics", background = X[[1]],
     used_iterations = sort(sample(1:length(out_mcmc$used_iterations),
                                   min(length(out_mcmc$used_iterations), 5e2))), W = W_ortho)
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
params_df <- as.data.frame(t(rbind(out_mcmc$chains$beta / mean_slopes,
                                   W_tilde_proj_mat %*% out_mcmc$chains$theta,
                                   out_mcmc$chains$sigsq_s,
                                   out_mcmc$chains$sigsq_mu,
                                   out_mcmc$chains$tausq,
                                   out_mcmc$chains$alpha,
                                   out_mcmc$chains$g0,
                                   out_mcmc$chains$gmax)))
names(params_df) <- c(paste0("X_", X_names), paste0("W_", W_names),
                      "sigsq_s", "sigsq_mu", "tausq", "alpha", "g0", "gmax")
CI_df <- apply(params_df, 2, quantile, probs = probs, na.rm = T)
CDF_X <- apply(out_mcmc$chains$beta, 1, function(x) sum(x > 0) / length(x))
CDF_W <- apply(W_tilde_proj_mat %*% out_mcmc$chains$theta, 1, function(x) sum(x > 0) / length(x))
names(CI_df) <- names(params_df)
t(CI_df)
layout(1)
par(mar = c(3.7, 2, 2, 1))
shown_params <- 1:(sum(nrow(out_mcmc$chains$beta), nrow(out_mcmc$chains$theta)))
CDF_shown <- c(CDF_X, CDF_W)
## plot separate ----
layout(matrix(1:2, 1, 2), widths = c(1, 1))
par(mar = c(3.9, 3.8, 2, 1), oma = c(0, 0, 0, 0))
X_scale_factor <- sqrt(sum(params_df[, 1:nrow(out_mcmc$chains$beta)]^2) /
                         (length(out_mcmc$chains$beta) - 1))
W_scale_factor <- sqrt(sum(params_df[, 1:nrow(out_mcmc$chains$theta) +
                                       nrow(out_mcmc$chains$beta)]^2) /
                         (length(out_mcmc$chains$theta) - 1))
## X
plot(1:nrow(out_mcmc$chains$beta), main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, nrow(out_mcmc$chains$beta) + 0.5),
     ylim = quantile(c(as.numeric(as.matrix(params_df[, 1:nrow(out_mcmc$chains$beta)] / X_scale_factor)),
                       as.numeric(as.matrix(params_df[, nrow(out_mcmc$chains$beta) +
                                                        1:nrow(out_mcmc$chains$theta)]) /
                                    W_scale_factor)),
                     probs = c(0.001, 0.999)
                     # probs = c(0, 1)
     )
)
abline(h = 0, lty = 2, col = "darkgray")
for(i in 1:nrow(out_mcmc$chains$beta)){
  density <- density(params_df[, i] / X_scale_factor)
  x <- density$y / diff(range(density$y)) * 0.33
  i_plot <- i
  polygon(x = c(i_plot - x, i_plot + rev(x)), c(density$x, rev(density$x)), col = "gray",
          lwd = 1 + 1.2 * as.numeric(min(CDF_shown[i], 1 - CDF_shown[i]) < 0.05), xpd = F)
}
axis_points <- c(-2, -1, 0, 1, 2) * signif(X_scale_factor, 1)
axis(side = 2, at = axis_points / X_scale_factor, labels = axis_points, xpd = T)
mtext(text = expression(beta), side = 2, line = 2.5, at = 0)
text(x = (nrow(out_mcmc$chains$beta) + 1)/2, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]),
     labels = "behavior", font = 2, xpd = T)
text(x = 1:nrow(out_mcmc$chains$beta),
     y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 45, xpd = T,
     labels = "distance to polygon", font = 1 + sapply(CDF_shown[1:nrow(out_mcmc$chains$beta)],
                                                       function(x) min(x, 1 - x) < 0.05))
mtext(text = "a)", side = 1, line = -37.5, at = -0.14, cex = 1.3)
points(1, truth$beta / mean_slopes / X_scale_factor, col = "darkred", pch = 19)
## W
par(mar = c(3.9, 1, 2, 3.8))
plot(1:nrow(out_mcmc$chains$theta), main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, nrow(out_mcmc$chains$theta) + 0.5),
     ylim = quantile(c(as.numeric(as.matrix(params_df[, 1:nrow(out_mcmc$chains$beta)] / X_scale_factor)),
                       as.numeric(as.matrix(params_df[, nrow(out_mcmc$chains$beta) + 1:nrow(out_mcmc$chains$theta)]) /
                                    W_scale_factor)),
                     probs = c(0.001, 0.999)
                     # probs = c(0, 1)
     )
)
abline(h = 0, lty = 2, col = "darkgray")
for(i in 1:nrow(out_mcmc$chains$theta)){
  density <- density(params_df[, i + nrow(out_mcmc$chains$beta)] / W_scale_factor)
  x <- density$y / diff(range(density$y)) * 0.33
  i_plot <- i
  polygon(x = c(i_plot - x, i_plot + rev(x)), c(density$x, rev(density$x)), col = "gray",
          lwd = 1 + 1.2 * as.numeric(min(CDF_shown[i + nrow(out_mcmc$chains$beta)],
                                         1 - CDF_shown[i + nrow(out_mcmc$chains$beta)]) < 0.05))
}
axis_points <- c(-2, -1, 0, 1, 2) * signif(W_scale_factor, 1)
axis(side = 4, at = axis_points / W_scale_factor, labels = axis_points, xpd = T)
mtext(text = expression(theta), side = 4, line = 2.5, at = -0.05)
text(x = (nrow(out_mcmc$chains$theta) + 1)/2, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]),
     labels = "recharge", font = 2, xpd = T)
text(x = 1:nrow(out_mcmc$chains$theta),
     y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 45, xpd = T,
     labels = W_names, font = 1 + sapply(CDF_shown[nrow(out_mcmc$chains$beta) + 1:nrow(out_mcmc$chains$theta)],
                                         function(x) min(x, 1 - x) < 0.05))
mtext(text = "b)", side = 1, line = -37.5, at = 0.4, cex = 1.3)
points(1:2, (W_tilde_proj_mat %*% truth$theta) / W_scale_factor, col = "darkred", pch = 19)
## dev.off ----
dev.off()
