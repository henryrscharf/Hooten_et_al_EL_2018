## libraries ----
library(reDyn); library(Matrix); library(spatstat); library(maptools)
## load buffalo data ----
load("../data/buffalo/buffalo_Cilla.RData")
# buffalo_proj <- buffalo_proj[buffalo_proj$POSIX > "2005-09-14" &
#                                buffalo_proj$POSIX < "2005-10-15", ]
buffalo_proj <- buffalo_proj[buffalo_proj$POSIX > "2005-09-30" &
                               buffalo_proj$POSIX < "2005-10-15", ]
obstimes <- as.numeric(buffalo_proj$POSIX) / 3600
## DEM ----
ELE_rast <- raster("../data/buffalo/buffaloELE_Cilla.grd")
# ## sabie river ----
# load("../data/buffalo/dist2sabie_ext.RData")
# ## covariates ----
# ## X
# X <- stack(list(
#   "elevation" = ELE_rast,
#   "slope" = terrain(ELE_rast, opt = "slope"),
#   "dist2sabie" = dist2sabie,
#   "sin_aspect" = sin(terrain(ELE_rast, opt = "aspect")),
#   "cos_aspect" = cos(terrain(ELE_rast, opt = "aspect")),
#   "roughness" = terrain(ELE_rast, opt = "roughness")
# ))
# X_names <- names(X)
# out_mcmc$slope[is.na(out_mcmc$slope)] <- mean(values(out_mcmc$slope), na.rm = T)
# out_mcmc$sin_aspect[is.na(out_mcmc$sin_aspect)] <- mean(values(out_mcmc$sin_aspect), na.rm = T)
# out_mcmc$cos_aspect[is.na(out_mcmc$cos_aspect)] <- mean(values(out_mcmc$cos_aspect), na.rm = T)
# out_mcmc$roughness[is.na(out_mcmc$roughness)] <- mean(values(out_mcmc$roughness), na.rm = T)
# X_scaled <- scale(X)
# ## standardize X based on slopes of gradients
# mean_slopes <- unlist(sapply(1:nrow(out_mcmc$chains$beta), function(layer){
#   mean(values(terrain(X[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
# }))
# X_grad_scaled <- X / mean_slopes
# unlist(sapply(1:nlayers(X_grad_scaled), function(layer){
#   mean(values(terrain(X_grad_scaled[[layer]], opt = "slope", unit = "tangent")), na.rm = T)
# }))
# gradX_scaled <- sapply(1:nlayers(X_grad_scaled), function(layer) rast_grad(X_grad_scaled[[layer]])
#                        [c('rast.grad.x', 'rast.grad.y')], simplify = F)
# ## W
# intercept <- raster(X_scaled)
# values(intercept) <- 1
# W <- stack(list(
#   "intercept" = intercept,
#   "elevation" = X_scaled$elevation,
#   "slope" = X_scaled$slope,
#   "near_sabie" = dist2sabie < 1e3,
#   "sin_aspect" = X_scaled$sin_aspect,
#   "cos_aspect" = X_scaled$cos_aspect,
#   "roughness" = X_scaled$roughness,
#   "elevation_slope" = X_scaled$elevation * X_scaled$slope
# ))
# W_names <- names(W)
## load fit ----
# load("../data/buffalo/model_fits/20181003_083829_buffalo_Cilla_reDyn_mcmc.RData")
# load("../data/buffalo/model_fits/extended_sabie/20181005_114612_buffalo_Cilla_reDyn_mcmc.RData")
load("../data/buffalo/model_fits/extended_sabie/sabie_only/20181006_043926_buffalo_Cilla_reDyn_mcmc_FULL.RData")
out_mcmc <- out_mcmcs_FULL[[1]]
## device ----
pdf(file = paste0("../fig/buffalo/priors_buffalo_reDyn_mcmc.pdf"))
## priors ----
mode_s <- out_mcmc$initial$sigsq_s; var_s <- out_mcmc$initial$sigsq_s * 1e4
alpha_s <- uniroot(f = function(alpha, mode = mode_s, var = var_s){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e5))$root
beta_s <- mode_s * (alpha_s + 1)
##
mode_mu <- out_mcmc$initial$sigsq_mu
var_mu <- mode_mu^2
alpha_mu <- uniroot(f = function(alpha, mode = mode_mu, var = var_mu){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_mu <- mode_mu * (alpha_mu + 1)
##
mode_tau <- out_mcmc$initial$tausq
var_tau <- mode_tau^2
alpha_tau <- uniroot(f = function(alpha, mode = mode_tau, var = var_tau){
  mode^2 * (alpha + 1)^2 / (alpha - 1)^2 / (alpha - 2) - var
}, interval = c(2, 1e10))$root
beta_tau <- mode_tau * (alpha_tau + 1)
##
log_priors <- list(
  "sigsq_s" = c(alpha_s, beta_s),
  "beta" = function(beta) sum(dnorm(beta, 0, 1e8, log = T)),
  "sigsq_mu" = function(sigsq_mu) dinvgamma(sigsq_mu, alpha_mu, beta_mu),
  "tausq" = function(tausq) dinvgamma(tausq, alpha_tau, beta_tau),
  "theta" = function(theta) sum(dnorm(theta, 0, 2, log = T)),
  "g0" = function(g0) dnorm(x = g0, mean = 0, sd = 2, log = T)
)
## plot priors ----
layout(matrix(1:(ceiling((length(out_mcmc$log_priors) +
                   nrow(out_mcmc$chains$beta) +
                   nrow(out_mcmc$chains$theta) - length(out_mcmc$fixed) - 4) / 3) * 3), ncol = 3))
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
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../fig/buffalo/traces_buffalo_reDyn_mcmc.pdf"))
## plot traces ----
plot(out_mcmc)
# plot(out_mcmc, type = "path")
## dev.off ----
dev.off()
## plot energetics [manual] ----
## device (map) ----
pdf(file = paste0("../fig/buffalo/energetics_map_buffalo_reDyn_mcmc.pdf"), height = 8)
## plot (map) ----
x <- out_mcmc
background <- ELE_rast
truth <- list()
used_iterations <- out_mcmc$used_iterations
plot_used_iterations <- used_iterations
used_iterations <- 1:dim(out_mcmc$chains$alpha)[2]
TIMES <- length(out_mcmc$data$times)
## calc
beta_median <- apply(matrix(out_mcmc$chains$beta[, used_iterations], nrow = nrow(out_mcmc$chains$beta)), 1, median)
theta_median <- apply(matrix(out_mcmc$chains$theta[, used_iterations], nrow = nrow(out_mcmc$chains$theta)), 1, median)
W_theta_median <- calc(W_ortho, function(x) x %*% theta_median)
mu_median <- apply(out_mcmc$chains$mu[, , used_iterations], 1:2, median)
mu_pwCI <- apply(out_mcmc$chains$mu, 1:2, quantile, probs = c(0.025, 0.975))
if(!exists("gs")){
  system.time({
    gs <- sapply(X = used_iterations,
                 FUN = function(iter){
                   g <- get_g_all(g0 = out_mcmc$chains$g0[iter], mu = out_mcmc$chains$mu[, , iter], W = W_ortho,
                                  theta = out_mcmc$chains$theta[, iter], alpha = out_mcmc$chains$alpha[iter],
                                  times = out_mcmc$data$times, gmax = out_mcmc$chains$gmax[iter])
                   return(g)
                 })
  })
} ## around 24 seconds for 10,000 used_iterations
g_pwCI <- apply(gs, 1, quantile, probs = c(0.025, 0.5, 0.975))
g_median <- g_pwCI[2, ]
g_color_range <- range(g_median)
g_color_range <- g_color_range / max(abs(g_color_range)) / 2 + 0.5
g_color_fn <- colorRamp(brewer.pal(n = 11, name = "RdYlBu"))
g_colors <- rgb(g_color_fn(g_median / max(abs(range(g_median))) / 2 + 0.5), maxColorValue = 255)
recharge_patch_col <- "#4d9221"
recharge_patch_cols <- RColorBrewer::brewer.pal(n = 4, name = "YlGn")[-1]
recharge_patch_alpha <- 0.6
## layout
layout(1)
par(mar = c(0.1, 0.1, 0.1, 3.5))
## background
image(background,
      col = gray.colors(1e2, start = 0.2),
      # col = NA,
      asp = 1, bty = "n", axes = F, main = "")
# title(main = names(background), line = 0, cex.main = 2)
image(W_theta_median > 0,
      col = c(NA, alpha(recharge_patch_col, recharge_patch_alpha)), add = T)
# image(W_theta_median > quantile(W_theta_median, probs = 0.95),
#       col = c(NA, alpha(recharge_patch_cols[1], recharge_patch_alpha)), add = T)
# image(W_theta_median > quantile(W_theta_median, probs = 0.99),
#       col = c(NA, alpha(recharge_patch_cols[3], recharge_patch_alpha)), add = T)
## path
lines(truth$mu, col = "darkred", lwd = 2)
segments(x0 = mu_median[-nrow(mu_median), 1], y0 = mu_median[-nrow(mu_median), 2],
         x1 = mu_median[-1, 1], y1 = mu_median[-1, 2], col = g_colors, lwd = 3)
## legend
legend_color_inds <- seq(g_color_range[2], g_color_range[1], l = 7)
legend_raster <- as.raster(matrix(rgb(g_color_fn(legend_color_inds), maxColorValue = 255), ncol = 1))
ext_X <- extent(W_theta_median)
rasterImage(legend_raster, xleft = ext_X@xmax + 0.01 * (ext_X@xmax - ext_X@xmin),
            xright = ext_X@xmax + 0.03 * (ext_X@xmax - ext_X@xmin),
            ybottom = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin),
            ytop = ext_X@ymax - 0.25 * (ext_X@ymax - ext_X@ymin), interpolate = F, xpd = T)
text(x = ext_X@xmax + 0.03 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymax - 0.20 * (ext_X@ymax - ext_X@ymin), labels = "g(t)", font = 2, xpd = T)
text(x = ext_X@xmax + 0.06 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymax - 0.25 * (ext_X@ymax - ext_X@ymin), labels = signif(max(g_median), 3),
     cex = 0.8, xpd = T)
text(x = ext_X@xmax + 0.06 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin), labels = signif(min(g_median), 3),
     cex = 0.8, xpd = T)
# text(x = ext_X@xmax + 0.045 * (ext_X@xmax - ext_X@xmin),
#      y = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin) +
#        (0.5 - g_color_range[1]) / diff(g_color_range) * 0.5 * (ext_X@ymax - ext_X@ymin),
#      labels = 0, cex = 0.8, xpd = T)
## dev.off ----
dev.off()
## device (time series) ----
pdf(file = paste0("../fig/buffalo/energetics_ts_buffalo_reDyn_mcmc.pdf"),
    width = 10, height = 6)
## plot (time series) ----
labels <- c("easting [m]", "northing [m]")
layout(matrix(1:4, 4, 1))
par(mar = c(2, 4, 1, 3), oma = c(0, 1.9, 0, 0))
for(d in 1:2){
  plot(out_mcmc$data$times, out_mcmc$chains$mu[, d, 1], type = "n", xlab = "", ylab = labels[d],
       main = "", bty = "n", xaxt = "n")
  # rect(par()$usr[1], par()$usr[2], par()$usr[3], par()$usr[4], col = scales::alpha("lightgray"), border = NA)
  marginal_bg <- extract(background, mu_median)
  unit_marginal_bg <- (marginal_bg - min(marginal_bg)) / diff(range(marginal_bg))
  scaled_marginal_bg <- unit_marginal_bg * diff(par()$usr[3:4]) + par()$usr[3]
  grad_bg <- diff(marginal_bg) / sqrt(diff(mu_median[, 1])^2 + diff(mu_median[, 2])^2)
  unit_grad_bg <- grad_bg / max(abs(grad_bg)) / 2 + 0.5
  lines(out_mcmc$data$times, scaled_marginal_bg, lwd = 2, col = "darkgray")
  if(d == 1){
    rug(out_mcmc$data$times[extract(W_theta_median, mu_median) > 0], lwd = 4, line = 1.75,
        col = recharge_patch_col)
  }
  polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)),
          y = c(mu_pwCI[1, , d], rev(mu_pwCI[2, , d])), col = scales::alpha("gray", 0.5), border = NA)
  segments(x0 = out_mcmc$data$times[-TIMES], y0 = mu_median[-TIMES, d],
           x1 = out_mcmc$data$times[-1], y1 = mu_median[-1, d], col = alpha(g_colors, 1),
           lwd = 3)
  # mu_subset <- sort(sample(used_iterations, min(length(used_iterations), 1e2)))
  # for(mu_i in mu_subset){
  #   segments(x0 = out_mcmc$data$times[-TIMES], y0 = out_mcmc$chains$mu[-TIMES, d, mu_i],
  #            x1 = out_mcmc$data$times[-1], y1 = out_mcmc$chains$mu[-1, d, mu_i], col = alpha(g_colors, 0.5),
  #            lwd = 2)
  # }
  points(out_mcmc$data$obstimes, out_mcmc$data$s[, d], pch = 19, cex = 0.5)
}
## g
plot(out_mcmc$data$times, g_median, type = "n", xlab = "time", ylab = "g", ylim = range(0, g_pwCI),
     bty = "n", xaxt = "n")
abline(h = 0, lty = 2, lwd = 2, col = "darkgray")
polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)), y = c(g_pwCI[1, ], rev(g_pwCI[3, ])),
        col = scales::alpha("gray", 0.5), border = NA, lty = 2)
segments(x0 = out_mcmc$data$times[-TIMES], y0 = g_median[-TIMES],
         x1 = out_mcmc$data$times[-1], y1 = g_median[-1], col = g_colors, lwd = 3)
## rho
plot(out_mcmc$data$times, pnorm(g_median), type = "n", xlab = "time", ylab = expression(rho), ylim = 0:1,
     bty = "n", xaxt = "n", yaxt = "n")
axis(2, c(0, 0.5, 1))
abline(h = 0.5, lty = 2, lwd = 2, col = "darkgray")
polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)), y = c(pnorm(-g_pwCI[1, ]), pnorm(-rev(g_pwCI[3, ]))),
        col = scales::alpha("gray", 0.5), border = NA, lty = 2)
segments(x0 = out_mcmc$data$times[-TIMES], y0 = pnorm(-g_median[-TIMES]),
         x1 = out_mcmc$data$times[-1], y1 = pnorm(-g_median[-1]), col = g_colors, lwd = 3)
points(out_mcmc$data$times, apply(out_mcmc$chains$z, 1, mean), pch = 16, cex = 0.6)
pretty_dates <- pretty(buffalo_proj$POSIX)
pretty_datetimes <- as.numeric(pretty_dates - min(pretty_dates)) /
  diff(range(as.numeric(pretty_dates))) *
  diff(range(out_mcmc$data$obstimes)) + min(out_mcmc$data$obstimes)
axis(1, pretty_datetimes, pretty_dates)
mtext("a)", side = 1, line = -43, xpd = T, at = 15087, cex = 1.3)
mtext("b)", side = 1, line = -31.8, xpd = T, at = 15087, cex = 1.3)
mtext("c)", side = 1, line = -20.2, xpd = T, at = 15087, cex = 1.3)
mtext("d)", side = 1, line = -8.8, xpd = T, at = 15087, cex = 1.3)
## dev.off ----
dev.off()
## device (map + time series) ----
pdf(file = paste0("../fig/buffalo/energetics_map_ts_buffalo_reDyn_mcmc.pdf"),
    width = 7, height = 7)
## plot (map + time series) ----
## layout
layout(matrix(c(5, 1:4), ncol = 1), height = c(3.3, 1, 1, 1, 1))
# layout(matrix(c(1, 2, 3, 4, rep(5, 4)), ncol = 2), width = c(3, 2))
labels <- c("easting [km]", "northing [km]")
par(mar = c(2, 4, 1, 3), oma = c(0, 1.9, 0, 0))
for(d in 1:2){
  plot(out_mcmc$data$times, out_mcmc$chains$mu[, d, 1] / 1e3, type = "n", xlab = "", ylab = labels[d],
       main = "", bty = "n", xaxt = "n")
  # rect(par()$usr[1], par()$usr[2], par()$usr[3], par()$usr[4], col = scales::alpha("lightgray"), border = NA)
  marginal_bg <- extract(background, mu_median)
  unit_marginal_bg <- (marginal_bg - min(marginal_bg)) / diff(range(marginal_bg))
  scaled_marginal_bg <- unit_marginal_bg * diff(par()$usr[3:4]) + par()$usr[3]
  grad_bg <- diff(marginal_bg) / sqrt(diff(mu_median[, 1])^2 + diff(mu_median[, 2])^2)
  unit_grad_bg <- grad_bg / max(abs(grad_bg)) / 2 + 0.5
  lines(out_mcmc$data$times, scaled_marginal_bg, lwd = 2, col = "darkgray")
  if(d == 1){
    rug(out_mcmc$data$times[extract(W_theta_median, mu_median) > 0], lwd = 4, line = 1.75,
        col = recharge_patch_col)
  }
  polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)),
          y = c(mu_pwCI[1, , d], rev(mu_pwCI[2, , d])) / 1e3, col = scales::alpha("gray", 0.5), border = NA)
  segments(x0 = out_mcmc$data$times[-TIMES], y0 = mu_median[-TIMES, d] / 1e3,
           x1 = out_mcmc$data$times[-1], y1 = mu_median[-1, d] / 1e3, col = alpha(g_colors, 1),
           lwd = 3)
  # mu_subset <- sort(sample(used_iterations, min(length(used_iterations), 1e2)))
  # for(mu_i in mu_subset){
  #   segments(x0 = out_mcmc$data$times[-TIMES], y0 = out_mcmc$chains$mu[-TIMES, d, mu_i],
  #            x1 = out_mcmc$data$times[-1], y1 = out_mcmc$chains$mu[-1, d, mu_i], col = alpha(g_colors, 0.5),
  #            lwd = 2)
  # }
  points(out_mcmc$data$obstimes, out_mcmc$data$s[, d] / 1e3, pch = 19, cex = 0.5)
}
## g
plot(out_mcmc$data$times, g_median, type = "n", xlab = "time", ylab = "g", ylim = range(0, g_pwCI),
     bty = "n", xaxt = "n")
abline(h = 0, lty = 2, lwd = 2, col = "darkgray")
polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)), y = c(g_pwCI[1, ], rev(g_pwCI[3, ])),
        col = scales::alpha("gray", 0.5), border = NA, lty = 2)
segments(x0 = out_mcmc$data$times[-TIMES], y0 = g_median[-TIMES],
         x1 = out_mcmc$data$times[-1], y1 = g_median[-1], col = g_colors, lwd = 3)
## rho
plot(out_mcmc$data$times, pnorm(g_median), type = "n", xlab = "time", ylab = expression(rho), ylim = 0:1,
     bty = "n", xaxt = "n", yaxt = "n")
axis(2, c(0, 0.5, 1))
abline(h = 0.5, lty = 2, lwd = 2, col = "darkgray")
polygon(x = c(out_mcmc$data$times, rev(out_mcmc$data$times)), y = c(pnorm(-g_pwCI[1, ]), pnorm(-rev(g_pwCI[3, ]))),
        col = scales::alpha("gray", 0.5), border = NA, lty = 2)
segments(x0 = out_mcmc$data$times[-TIMES], y0 = pnorm(-g_median[-TIMES]),
         x1 = out_mcmc$data$times[-1], y1 = pnorm(-g_median[-1]), col = g_colors, lwd = 3)
points(out_mcmc$data$times, apply(out_mcmc$chains$z, 1, mean), pch = 16, cex = 0.6)
pretty_dates <- as.POSIXct(paste0("2005-10-", c(1, 5, 10, 15)))
pretty_datetimes <- as.numeric(pretty_dates - min(pretty_dates)) /
  diff(range(as.numeric(pretty_dates))) *
  diff(range(out_mcmc$data$obstimes)) + min(out_mcmc$data$obstimes)
axis(1, pretty_datetimes, pretty_dates)
mtext("a)", side = 1, line = -50.5, xpd = T, at = 313285, cex = 1.1)
mtext("b)", side = 1, line = -26, xpd = T, at = 313285, cex = 1.1)
mtext("c)", side = 1, line = -19.5, xpd = T, at = 313285, cex = 1.1)
mtext("d)", side = 1, line = -11.5, xpd = T, at = 313285, cex = 1.1)
mtext("e)", side = 1, line = -4.5, xpd = T, at = 313285, cex = 1.1)
##
par(mar = c(0.1, 0.1, 0.1, 3.5))
## background
image(background, col = gray.colors(1e2), asp = 1, bty = "n", axes = F, main = "")
# title(main = names(background), line = 0, cex.main = 2)
image(W_theta_median > 0, col = c(NA, alpha(recharge_patch_col, 0.6)), add = T)
# lines(rasterToContour(W_theta_median, levels = 0), lty = 1, lwd = 1.5)
## path
segments(x0 = mu_median[-nrow(mu_median), 1], y0 = mu_median[-nrow(mu_median), 2],
         x1 = mu_median[-1, 1], y1 = mu_median[-1, 2], col = g_colors, lwd = 3)
## legend
legend_color_inds <- seq(g_color_range[2], g_color_range[1], l = 7)
legend_raster <- as.raster(matrix(rgb(g_color_fn(legend_color_inds), maxColorValue = 255), ncol = 1))
ext_X <- extent(W_theta_median)
rasterImage(legend_raster, xleft = ext_X@xmax + 0.01 * (ext_X@xmax - ext_X@xmin),
            xright = ext_X@xmax + 0.03 * (ext_X@xmax - ext_X@xmin),
            ybottom = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin),
            ytop = ext_X@ymax - 0.25 * (ext_X@ymax - ext_X@ymin), interpolate = F, xpd = T)
text(x = ext_X@xmax + 0.07 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymax - 0.20 * (ext_X@ymax - ext_X@ymin), labels = "g(t)", font = 2, xpd = T)
text(x = ext_X@xmax + 0.10 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymax - 0.25 * (ext_X@ymax - ext_X@ymin), labels = signif(max(g_median), 3),
     cex = 0.8, xpd = T)
text(x = ext_X@xmax + 0.09 * (ext_X@xmax - ext_X@xmin),
     y = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin), labels = signif(min(g_median), 3),
     cex = 0.8, xpd = T)
# text(x = ext_X@xmax + 0.045 * (ext_X@xmax - ext_X@xmin),
#      y = ext_X@ymin + 0.25 * (ext_X@ymax - ext_X@ymin) +
#        (0.5 - g_color_range[1]) / diff(g_color_range) * 0.5 * (ext_X@ymax - ext_X@ymin),
#      labels = 0, cex = 0.8, xpd = T)
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../fig/buffalo/data_map_only_buffalo_reDyn_mcmc.pdf"),
    height = 9, width = 6)
## plot data (map only) ----
## layout
par(mar = c(0.1, 0.1, 0.1, 0.1))
## map
image(out_mcmc$elevation, col = gray.colors(1e2, start = 0.2), asp = 1,
      bty = "n", axes = F, main = "")
points(buffalo_proj@coords, pch = 19, cex = 0.6)
# arrows(x0 = buffalo_proj@coords[-361, 1], y0 = buffalo_proj@coords[-361, 2],
#        x1 = buffalo_proj@coords[-1, 1], y1 = buffalo_proj@coords[-1, 2], length = 0.1)
# points(buffalo_proj@coords, pch = 19, cex = 0.6, col = colorRampPalette(RColorBrewer::brewer.pal(5, "YlOrBr")[-1])(nrow(buffalo_proj)))
# points(t(buffalo_proj@coords[1, ]), pch = 19, col = "darkgreen", cex = 2)
# points(t(buffalo_proj@coords[361, ]), pch = 18, col = "darkblue", cex = 2)
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../fig/buffalo/data_buffalo_reDyn_mcmc.pdf"),
    height = 11 * 0.48, width = 11)
## plot data (with marginals) ----
## layout
layout(matrix(c(2, 3, 1, 1), 2, 2), widths = c(4, 3))
par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0, 2, 0, 0))
## map
image(out_mcmc$elevation, col = gray.colors(1e2), asp = 1, bty = "n", axes = F, main = "")
points(out_mcmc$data$s, pch = 19, cex = 0.6)
mtext("c)", side = 1, line = -31.5, xpd = T, at = 434000, cex = 1.3)
par(mar = c(3, 4, 0.1, 1))
## marginals
plot(out_mcmc$data$obstimes, out_mcmc$data$s[, 1], xlab = "", ylab = "",
     xaxt = "n", type = "n", bty = "n")
mtext("easting [m]", 2, line = 3)
marginal_bg <- extract(out_mcmc$elevation, out_mcmc$data$s)
unit_marginal_bg <- (marginal_bg - min(marginal_bg)) / diff(range(marginal_bg))
scaled_marginal_bg <- unit_marginal_bg * diff(par()$usr[3:4]) + par()$usr[3]
lines(out_mcmc$data$obstimes, scaled_marginal_bg, lwd = 2, col = "darkgray")
points(out_mcmc$data$obstimes, out_mcmc$data$s[, 1], pch = 19, cex = 0.6)
plot(out_mcmc$data$obstimes, out_mcmc$data$s[, 2], xlab = "", ylab = "",
     xaxt = "n", type = "n", bty = "n")
mtext("northing [m]", 2, line = 3)
scaled_marginal_bg <- unit_marginal_bg * diff(par()$usr[3:4]) + par()$usr[3]
lines(out_mcmc$data$obstimes, scaled_marginal_bg, lwd = 2, col = "darkgray")
points(out_mcmc$data$obstimes, out_mcmc$data$s[, 2], pch = 19, cex = 0.6)
pretty_dates <- pretty(buffalo_proj$POSIX)
pretty_dates <- c(as.POSIXct("2005-10-01"), as.POSIXct("2005-10-05"), as.POSIXct("2005-10-10"), as.POSIXct("2005-10-14"))
pretty_datetimes <- as.numeric(pretty_dates - min(pretty_dates)) /
  diff(range(as.numeric(pretty_dates))) *
  diff(range(out_mcmc$data$obstimes)) + min(out_mcmc$data$obstimes)
axis(1, pretty_datetimes, pretty_dates)
mtext("a)", side = 1, line = -28.5, xpd = T, at = 15084.8, cex = 1.3)
mtext("b)", side = 1, line = -13, xpd = T, at = 15084.9, cex = 1.3)
## dev.off ----
dev.off()
## device ----
pdf(file = paste0("../fig/buffalo/correlations_buffalo_reDyn_mcmc.pdf"))
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
## dev.off ----
dev.off()


## device ----
pdf(file = paste0("../fig/buffalo/summaries_buffalo_reDyn_mcmc.pdf"), height = 5.3, width = 8.5)
## posterior summaries ----
probs <- c(0.025, 0.05, 0.5, 0.95, 0.975)
params_df <- as.data.frame(t(rbind(out_mcmc$chains$beta,
                                   W_tilde_proj_mat %*% out_mcmc$chains$theta,
                                   out_mcmc$chains$sigsq_s,
                                   out_mcmc$chains$sigsq_mu,
                                   out_mcmc$chains$tausq,
                                   out_mcmc$chains$alpha,
                                   out_mcmc$chains$g0,
                                   out_mcmc$chains$gmax)))
names(params_df) <- c(paste0("X_", X_names), paste0("W_", W_names),
                      "sigsq_s", "sigsq_mu", "tausq", "g0", "gmax")
CI_df <- apply(params_df, 2, quantile, probs = probs, na.rm = T)
CDF_X <- apply(out_mcmc$chains$beta, 1, function(x) sum(abs(x) > 0) / length(x))
CDF_W <- apply(W_tilde_proj_mat %*% out_mcmc$chains$theta, 1, function(x) sum(abs(x) > 0) / length(x))
names(CI_df) <- names(params_df)
t(CI_df)
layout(1)
par(mar = c(3.7, 2, 2, 1))
shown_params <- 1:(sum(nrow(out_mcmc$chains$beta), nrow(out_mcmc$chains$theta)))
CDF_shown <- c(CDF_X, CDF_W)
## plot separate ----
X_names[X_names == "dist2sabie"] <- "dist. to surface water"
W_names[W_names == "near_sabie"] <- "<0.5km to surface water"
# W_names[W_names_sabie = "elevation_slope"] <- "elevation * slope"
layout(matrix(1:2, 1, 2), widths = c(1, 1))
par(mar = c(3.9, 3.8, 2, 1), oma = c(0, 0, 0, 0))
X_scale_factor <- sqrt(sum(params_df[, 1:nrow(out_mcmc$chains$beta)]^2) / (length(out_mcmc$chains$beta) - 1))
W_scale_factor <- sqrt(sum(params_df[, 1:nrow(out_mcmc$chains$theta) + nrow(out_mcmc$chains$beta)]^2) / (length(out_mcmc$chains$theta) - 1))
## X
plot(1:nrow(out_mcmc$chains$beta), main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, nrow(out_mcmc$chains$beta) + 0.5),
     ylim = quantile(c(as.numeric(as.matrix(params_df[, 1:nrow(out_mcmc$chains$beta)] / X_scale_factor)),
                       as.numeric(as.matrix(params_df[, nrow(out_mcmc$chains$beta) + 1:nrow(out_mcmc$chains$theta)]) /
                                    W_scale_factor)),
                     # probs = c(0.001, 0.999)
                     probs = c(0, 1)
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
axis_points <- c(-2, -1, 0, 1, 2) * signif(X_scale_factor, 1) * 1.25
axis(side = 2, at = axis_points / X_scale_factor, labels = axis_points, xpd = T)
mtext(text = expression(beta), side = 2, line = 2.5, at = 0)
text(x =  (1 + nrow(out_mcmc$chains$beta)) / 2, y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]),
     labels = "behavior", font = 2, xpd = T)
text(x = 1:nrow(out_mcmc$chains$beta),
     y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 30, xpd = T,
     labels = X_names, font = 1 + sapply(CDF_shown[1:nrow(out_mcmc$chains$beta)], function(x) min(x, 1 - x) < 0.05))
mtext(text = "a)", side = 1, line = -22.5, at = 0.25, cex = 1.3)
## W
par(mar = c(3.9, 1, 2, 3.8))
plot(1:nrow(out_mcmc$chains$theta), main = "", type = "n", bty = "n",
     yaxt = "n", xaxt = "n", ylab = "", xlab = "",
     xlim = c(0.5, nrow(out_mcmc$chains$theta) + 0.5),
     ylim = quantile(c(as.numeric(as.matrix(params_df[, 1:nrow(out_mcmc$chains$beta)] / X_scale_factor)),
                       as.numeric(as.matrix(params_df[, nrow(out_mcmc$chains$beta) + 1:nrow(out_mcmc$chains$theta)]) /
                                    W_scale_factor)),
                     # probs = c(0.001, 0.999)
                     probs = c(0, 1)
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
axis_points <- c(-2, -1, 0, 1, 2) * signif(W_scale_factor, 1) * 1
axis(side = 4, at = axis_points / W_scale_factor, labels = axis_points, xpd = T)
mtext(text = expression(theta), side = 4, line = 2.5, at = -0.05)
text(x = (1 + nrow(out_mcmc$chains$theta)) / 2,
     y = par()$usr[4] + 0.05 * diff(par()$usr[3:4]), labels = "recharge", font = 2, xpd = T)
text(x = 1:nrow(out_mcmc$chains$theta),
     y = par()$usr[3] - 0.05 * diff(par()$usr[3:4]), srt = 30, xpd = T,
     labels = W_names, font = 1 + sapply(CDF_shown[nrow(out_mcmc$chains$beta) + 1:nrow(out_mcmc$chains$theta)],
                                         function(x) min(x, 1 - x) < 0.05))
mtext(text = "b)", side = 1, line = -22.5, at = 0.4, cex = 1.3)
## dev.off ----
dev.off()
## ----
# examine scores ----
## ORIGINAL ----
# ## Sabie only models  ----
# load("../data/buffalo/model_fits/sabie_only/20181004_142116_buffalo_Cilla_reDyn_mcmc_FULL.RData")
# out_mcmcs_FULL_sabie <- out_mcmcs_FULL
# W_ortho_sabie <- W_ortho
# W_tilde_proj_mat_sabie <- W_tilde_proj_mat
# W_names_sabie <- W_names
# X_names_sabie <- X_names
# load("../data/buffalo/model_fits/sabie_only/20181005_034609_buffalo_Cilla_reDyn_mcmc_M1.RData")
# out_mcmcs_M1_sabie <- out_mcmcs_M1
# ## Sabie only scores ----
# load("../data/buffalo/model_fits/sabie_only/20181004_142116_buffalo_Cilla_reDyn_scores_FULL.RData")
# lds_FULL_mat_sabie <- lds_FULL_mat
# load("../data/buffalo/model_fits/sabie_only/20181005_034609_buffalo_Cilla_reDyn_scores_M1.RData")
# lds_M1_mat_sabie <- lds_M1_mat
# ## all covariates models ----
# load("../data/buffalo/model_fits/20181002_184831_buffalo_Cilla_reDyn_mcmc_FULL.RData")
# load("../data/buffalo/model_fits/20181003_065124_buffalo_Cilla_reDyn_mcmc_M0.RData")
# load("../data/buffalo/model_fits/20181003_163738_buffalo_Cilla_reDyn_mcmc_M1.RData")
# ## all covariate scores ----
# load("../data/buffalo/model_fits/20181002_184831_buffalo_Cilla_reDyn_scores_FULL.RData")
# load("../data/buffalo/model_fits/20181003_065124_buffalo_Cilla_reDyn_scores_M0.RData")
# load("../data/buffalo/model_fits/20181003_163738_buffalo_Cilla_reDyn_scores_M1.RData")
## EXTENDED ----
## extended Sabie only models ----
load("../data/buffalo/model_fits/extended_sabie/sabie_only/20181006_043926_buffalo_Cilla_reDyn_mcmc_FULL.RData")
out_mcmcs_FULL_sabie <- out_mcmcs_FULL
W_ortho_sabie <- W_ortho
W_tilde_proj_mat_sabie <- W_tilde_proj_mat
W_names_sabie <- W_names
X_names_sabie <- X_names
load("../data/buffalo/model_fits/extended_sabie/sabie_only/20181006_141105_buffalo_Cilla_reDyn_mcmc_M1.RData")
out_mcmcs_M1_sabie <- out_mcmcs_M1
## extended Sabie only scores ----
load("../data/buffalo/model_fits/extended_sabie/sabie_only/20181006_043926_buffalo_Cilla_reDyn_scores_FULL.RData")
lds_FULL_mat_sabie <- lds_FULL_mat
load("../data/buffalo/model_fits/extended_sabie/sabie_only/20181006_141105_buffalo_Cilla_reDyn_scores_M1.RData")
lds_M1_mat_sabie <- lds_M1_mat
## extended all covariates models ----
load("../data/buffalo/model_fits/extended_sabie/20181007_125902_buffalo_Cilla_reDyn_mcmc_FULL.RData")
load("../data/buffalo/model_fits/20181003_065124_buffalo_Cilla_reDyn_mcmc_M0.RData")
load("../data/buffalo/model_fits/extended_sabie/20181007_091919_buffalo_Cilla_reDyn_mcmc_M1.RData")
## extended all covariate scores ----
load("../data/buffalo/model_fits/extended_sabie/20181007_125902_buffalo_Cilla_reDyn_scores_FULL.RData")
load("../data/buffalo/model_fits/20181003_065124_buffalo_Cilla_reDyn_scores_M0.RData")
load("../data/buffalo/model_fits/extended_sabie/20181007_091919_buffalo_Cilla_reDyn_scores_M1.RData")
## order scores ----
scores_arr <- array(c(lds_FULL_mat, lds_M0_mat, lds_M1_mat,
                      lds_FULL_mat_sabie, lds_M1_mat_sabie), dim = c(dim(lds_FULL_mat), 5))
neg_mean_scores_fold <- apply(-scores_arr, 2:3, function(lds)
  mean(lds[lds != Inf], na.rm = T))
models <- c("FULL", "M0", "M1", "FULL_sabie", "M1_sabie")
colnames(neg_mean_scores_fold) <- models
apply(neg_mean_scores_fold, 2, mean)
## device ----
pdf(file = "../fig/buffalo/model_scores_buffalo_ext.pdf")
## plot ----
# boxplot(neg_mean_scores_fold, xlab = "model", ylab = "-log predictive density")
matplot(neg_mean_scores_fold, type = "l", lty = 1, lwd = 2,
        xlab = "fold", ylab = "-log predictive density")
legend("topleft", col = 1:length(models), lty = 1, legend = models, bty = "n", lwd = 2, ncol = 2)
## dev.off ----
dev.off()
## look at individual folds ----
library(reDyn)
pdf("../fig/buffalo/summaries_extended_sabie_only.pdf", height = 4.9, width = 8.5)
plot(out_mcmcs_FULL_sabie[[1]], type = "effects", W_tilde_proj_mat = W_tilde_proj_mat_sabie,
     W_names = W_names_sabie, X_names = X_names_sabie)
dev.off()
pdf("../fig/buffalo/summaries_extended_all_covariates.pdf", height = 4.9, width = 8.5)
plot(out_mcmcs_FULL[[1]], type = "effects", W_tilde_proj_mat = W_tilde_proj_mat,
     W_names = W_names, X_names = X_names)
dev.off()
# plot(out_mcmcs_M1[[1]], type = "effects", W_tilde_proj_mat = W_tilde_proj_mat,
#      W_names = W_names, X_names = X_names)
##
pdf("../fig/buffalo/energetics_extended_sabie_only.pdf", height = 4.9, width = 8.5)
plot(out_mcmcs_FULL_sabie[[1]], type = "energetics", W = W_ortho_sabie)
dev.off()
pdf("../fig/buffalo/energetics_extended_all_covariates.pdf", height = 4.9, width = 8.5)
plot(out_mcmcs_FULL[[2]], type = "energetics", W = W_ortho)
dev.off()
