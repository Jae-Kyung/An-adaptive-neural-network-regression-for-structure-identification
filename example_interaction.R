########################################################
# Load data
########################################################
rm(list = ls())
source("main_code.R")
source("sub_code.R")
load("interaction_data_N400")
sample_size = 400
number_predictors = 4
set.seed(20220919)
index = sample(1 : 100, 1)
x = data[[index]]$x
y = data[[index]]$y
f = data[[index]]$f
########################################################
# setting (R pakage msir is required for initialization)
########################################################
library(msir)
number_nodes = sample_size / 10
initials = initial_ABNN(x, y, number_nodes, number_dir = number_predictors)
weights = initials$weights
bias = initials$bias
beta = rep(0, number_nodes)
constant = 0
order = 3
number_lambdas = 50
lambdas = c(exp(seq(log(1e-05), log(1), length = number_lambdas)))
number_lambdag = 50
lambdag = c(exp(seq(log(1e-05), log(1), length = number_lambdag)))
########################################################
# fitting
########################################################
fits = list()
bics = matrix(0, number_lambdas, number_lambdag)
aics = matrix(0, number_lambdas, number_lambdag)
for (j in 1 : number_lambdag)
{
   cat("\n", j, "th runs")
   fit = ABNN(predictors = x, response = y,
              lambdas = lambdas, lambdag = lambdag[j], 
              bias = bias, order = order,
              weights = weights, 
              beta = beta, constant = constant, 
              number_hidden_nodes = number_nodes)
   fits[[j]] = fit
   bics[, j] = fit$bic_list
   # initial feed forward
   weights = fits[[j]][[1]]$weights
   bias = fits[[j]][[1]]$bias
   beta = fits[[j]][[1]]$beta
   constant = fits[[j]][[1]]$constant
}
########################################################
# results
########################################################
(bic_index = which(bics == min(bics), arr.ind = T))
opt_fit = fits[[bic_index[1]]][[bic_index[2]]]
# opt_fit = fits[[opt_index]]
# MSE
mean((opt_fit$fitted_values - f)^2)
# terms included in ABNN
ABNN_terms(opt_fit)
########################################################
# visualization
########################################################
# partial dependence plot on 1D
par(mfrow = c(1, 1), pty = "s")
x_grid = seq(0, 1, length = 1000)
for (j in 1 : 4)
{
   partial_function = fitted_partial_1D(opt_fit, x_grid, x, index = j)
   plot(x_grid, partial_function, type = "l", col = "blue",
        xlab = bquote(~x[.(j)]), ylab = "", cex.lab = 1.5)
}
# partial dependence plot on 2D
# library(RColorBrewer)
library(viridis) 
library(lattice)
n_grid = 50
x1_grid = seq(0, 1, length = n_grid)
x2_grid = seq(0, 1, length = n_grid)
x_grid = as.matrix(expand.grid(x1_grid, x2_grid))
# f34
index = c(3, 4)
partial_function = fitted_partial_2D(opt_fit, x_grid, x, index = index)
# plot(scaled_xgrid, a * sd(y) + mean(y), type = "l")
fitted_f34 = levelplot(t(matrix(partial_function, n_grid, n_grid)),
                       col.regions = viridis(100),
                       row.values = seq(0, 1, length = n_grid), 
                       column.values = seq(0, 1, length = n_grid), 
                       ylim = c(0, 1), xlim = c(0, 1), 
                       xlab = bquote(~x[3]), ylab = bquote(~x[4]))
fitted_f34
# f13
index = c(1, 3)
partial_function = fitted_partial_2D(opt_fit, x_grid, x, index = index)
fitted_f13 = levelplot(t(matrix(partial_function, n_grid, n_grid)),
                       col.regions = viridis(100),
                       row.values = seq(0, 1, length = n_grid), 
                       column.values = seq(0, 1, length = n_grid), 
                       ylim = c(0, 1), xlim = c(0, 1),
                       xlab = bquote(~x[1]), ylab = bquote(~x[3]))
fitted_f13
# f12
index = c(1, 2)
scaled_xgrid = matrix(0, n_grid^2, 2)
partial_function = fitted_partial_2D(opt_fit, x_grid, x, index = index)
fitted_f12 = levelplot(t(matrix(partial_function, n_grid, n_grid)),
                       col.regions = viridis(100),
                       row.values = seq(0, 1, length = n_grid), 
                       column.values = seq(0, 1, length = n_grid), 
                       ylim = c(0, 1), xlim = c(0, 1), 
                       xlab = bquote(~x[1]), ylab = bquote(~x[2]))
fitted_f12
