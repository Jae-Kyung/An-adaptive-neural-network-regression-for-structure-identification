
# functions for the two-way interaction example
g1 = function(z)
   return(z)
g2 = function(z)
   return((2 * z - 1)^2)
g3 = function(z)
   return(sin(2 * pi * z) / (2 - sin(2 * pi * z)))
g4 = function(z)
   return(0.1 * sin(2 * pi * z) + 0.2 * cos(2 * pi * z) + 0.3 * sin(2 * pi * z)^2 + 0.4 * cos(2 * pi * z)^3 + 0.5 * sin(2 * pi * z)^3)

# for partial dependence plot on 1D
fitted_partial_1D = function(opt_fit, x_grid, predictors, index = 1)
{
   grid_size = length(x_grid)
   number_nodes = length(opt_fit$beta)
   partial_function = rep(opt_fit$constant, grid_size)
   for (g in 1 : grid_size)
   {
      for (m in 1 : number_nodes)
      {
         plane = opt_fit$bias[m] + predictors[, -index] %*% opt_fit$weights[-index, m] + opt_fit$weights[index, m] * x_grid[g]
         active_plane = bspline_general(h = plane, order = 3)
         partial_function[g] = partial_function[g] + opt_fit$beta[m] * mean(active_plane)
      }
   }
   return(partial_function)
}

# for partial dependence plot on 2D
fitted_partial_2D = function(opt_fit, x_grid, predictors, index = c(1, 2))
{
   grid_size = nrow(x_grid)
   index_active = (1 : length(opt_fit$beta))[opt_fit$beta != 0]
   number_nodes = length(index_active)
   partial_function = rep(opt_fit$constant, grid_size)
   for (g in 1 : grid_size)
   {
      for (a in 1 : number_nodes)
      {
         m = index_active[a]
         plane = opt_fit$bias[m] + predictors[, -index] %*% opt_fit$weights[-index, m] + sum(opt_fit$weights[index, m] * x_grid[g, ])
         active_plane = bspline_general(h = plane, order = 3)
         partial_function[g] = partial_function[g] + opt_fit$beta[m] * mean(active_plane)
      }
   }
   return(partial_function)
}