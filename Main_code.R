
ABNN = function(predictors, 
                response, 
                number_hidden_nodes, 
                order = 2,
                lambdas = 0,
                lambdag = 0,
                constant = NULL,
                bias = NULL,
                weights = NULL,
                beta = NULL,
                thres = 1e-10,
                max_iter = 1000, 
                epsilon = 1e-05,
                verbose = F)
{
   dlm = list()
   sample_size = length(response)
   number_predictors = ncol(predictors)
   # initial values
   if (is.null(weights))
   {
      initials = initial_ABNN(x, y, number_nodes = number_nodes, number_dir = number_predictors)
      weights = initials$weights
      bias = initials$bias
      beta = rep(0, number_nodes)
      constant = 0
   }
   # feed forward
   planes = cbind(1, predictors) %*% rbind(bias, weights)
   active_planes = bspline_mat_general(planes, order)
   active_planes_derivative = bspline_derivative_mat_general(planes, order)
   active_nodes = rep(T, number_hidden_nodes)
   # back propagation
   number_lambdas = length(lambdas)
   lambda_list = rep(0, number_lambdas)
   bic_list = rep(0, number_lambdas)
   aic_list = rep(0, number_lambdas)
   rss_list = rep(0, number_lambdas)
   dimension_list = rep(0, number_lambdas)
   # fitted values
   fitted_values = constant + active_planes %*% beta
   residuals = response - fitted_values
   rss = 0.5 * mean(residuals^2)
   weights_temp = weights
   for (lambda_index in 1 : number_lambdas)
   {
      lambda = lambdas[lambda_index]
      store_rss = Inf
      if (verbose)
         cat("\n lambda index =", lambda_index, "\n")
      # start iterations
      for (iter in 1 : max_iter)
      {
         # update constant
         partial_residuals = residuals + constant
         constant = mean(partial_residuals)
         residuals = partial_residuals - constant
         # update m-th nodes
         for (m in 1 : number_hidden_nodes)
         {
            if (active_nodes[m])
            {
               # update beta
               partial_residuals = residuals + beta[m] * active_planes[, m]
               tilde_z = active_planes[, m]
               tilde_y = partial_residuals
               if (sum(abs(tilde_z)) == 0)
                  beta[m] = 0
               else
                  beta[m] = univariate_lasso_rss(tilde_y, tilde_z, lambda, 1e-10)
               residuals = partial_residuals - beta[m] * active_planes[, m]
               # update bias and weights
               partial_residuals = residuals + beta[m] * active_planes[, m]
               partial_planes = planes[, m] - bias[m]
               tilde_z = active_planes_derivative[, m] * beta[m]
               tilde_y = residuals + tilde_z * bias[m]
               if (sum(abs(tilde_z)) == 0)
                  bias[m] = 0
               else
                  bias[m] = mean(tilde_y * tilde_z) / mean(tilde_z^2)
               planes[, m] = partial_planes + bias[m]
               active_planes[, m] = bspline_general(planes[, m], order)
               active_planes_derivative[, m] = bspline_derivative_general(planes[, m], order)
               residuals = partial_residuals - beta[m] * active_planes[, m]
               # update weights
               weights_temp = weights
               for (j in 1 : number_predictors)
               {
                  if (weights[j, m] != 0)
                  {
                     # cat("\n ind:", j, m, "\n")
                     partial_residuals = residuals + beta[m] * active_planes[, m]
                     partial_planes = planes[, m] - weights[j, m] * predictors[, j]
                     tilde_z = beta[m] * active_planes_derivative[, m] * predictors[, j]
                     tilde_y = residuals + tilde_z * weights[j, m]
                     weighted_lambda = lambdag * (sum((abs(weights[-j, m])))) / (sum(abs(weights[, m])) * abs(weights[j, m]))
                     if (sum(abs(tilde_z)) == 0)
                        weights[j, m] = 0
                     else
                        weights[j, m] = univariate_lasso_rss(tilde_y, tilde_z, weighted_lambda, thres)
                     planes[, m] = partial_planes + weights[j, m] * predictors[, j]
                     active_planes[, m] = bspline_general(planes[, m], order)
                     active_planes_derivative[, m] = bspline_derivative_general(planes[, m], order)
                     residuals = partial_residuals - beta[m] * active_planes[, m]
                  }
               }
            }
         }
         active_nodes = (beta != 0)
         rss = 0.5 * mean(residuals^2)
         penalty = 0
         impurity = calculate_impurity(weights)
         penalty = lambdag * sum(impurity) + lambda * sum(abs(beta))
         rss = rss + penalty
         if (abs(rss - store_rss) < epsilon)
            break
         store_rss = rss
      }
      fitted_values = constant + active_planes %*% beta
      residuals = response - fitted_values
      rss_list[lambda_index] = mean(residuals^2)
      dimension = 0
      for (m in 1 : number_nodes)
         dimension = dimension + sum(weights[, m] != 0) * sum(beta[m] != 0)
      dimension_list[lambda_index] = dimension
      bic_list[lambda_index] = sample_size * log(rss_list[lambda_index]) + dimension_list[lambda_index] * log(sample_size)
      aic_list[lambda_index] = sample_size * log(rss_list[lambda_index]) + dimension_list[lambda_index] * 2
      dlm[[lambda_index]] = list()
      dlm[[lambda_index]]$constant = constant
      dlm[[lambda_index]]$beta = beta
      dlm[[lambda_index]]$weights = weights
      dlm[[lambda_index]]$bias = bias
      dlm[[lambda_index]]$planes = planes
      dlm[[lambda_index]]$active_planes = active_planes
      dlm[[lambda_index]]$fitted_values = fitted_values
   }
   # browser()
   dlm$bic_list = bic_list
   dlm$aic_list = aic_list
   dlm$rss_list = rss_list
   dlm$dimension_list = dimension_list
   return (dlm)
}

# calculate the impurity for a node
calculate_impurity_node = function(weights_node)
{
   M0 = length(weights_node)
   scale_term = sum(abs(weights_node))
   if (scale_term == 0)
      return(0)
   w = abs(weights_node) / scale_term
   impurity = 1 - sum(w^2)
   return(impurity)
}

# calculate the impurity for nodes
calculate_impurity = function(weights)
{
   number_nodes = ncol(weights)
   impurity = rep(0, number_nodes)
   for (m in 1 : number_nodes)
   {
      impurity[m] = calculate_impurity_node(weights[, m])
   }
   return(impurity)
}

# generate the predicted values for new x (new_predictors)
ABNN_predict = function(fit, new_predictors, order = 3, scale = T)
{
   if (scale)
      scaled_x = scale(new_predictors)
   planes = cbind(1, scaled_x) %*% rbind(fit$bias, fit$weights)
   active_planes = bspline_mat_general(planes, order)
   fitted_values = fit$constant + active_planes %*% fit$beta
   return(fitted_values)
}

# given bias and weights, calculate the matrix for the active planes
active_planes_dlm = function(grid_points, bias, weights, order = 4)
{
   planes = cbind(1, grid_points) %*% rbind(bias, weights)
   active_planes = bspline_mat_general(planes[[l + 1]], order)
   return(list(planes = planes, active_planes = active_planes))
}

# generate initial values based on msir method
initial_ABNN = function(x, y, number_nodes, number_dir = 1, order = 2)
{
   number_predictors = ncol(x)
   weights = matrix(0, number_predictors, 0)
   bias = c()
   msir_fit = msir(x, y, nslices = 2, G = rep(5, 2))
   for (d in 1 : number_dir)
   {
      dir = msir_fit$std.basis[, d]
      number_node = number_nodes / number_dir
      z = x %*% dir
      centers = quantile(z, probs = seq(0, 1, length = number_node))
      width = rep(0, number_node)
      # browser()
      for (m in 1 : number_node)
      {
         diff_centers = diff(centers)
         if (m == 1)
            width[m] = diff_centers[1]
         else if (m == number_node)
            width[m] = diff_centers[number_node - 1]
         else
            width[m] = max(diff_centers[m - 1], diff_centers[m])
      }
      w = 2 / (width * 4)
      for (m in 1 : number_node)
         weights = cbind(weights, w[m] * dir)
      b = -centers * w
      bias = c(bias, b)
   }
   return(list(weights = weights, bias = bias))
}

# check the interactions
ABNN_terms = function(opt_fit)
{
   weights = opt_fit$weights
   number_predictors = nrow(weights)
   active_index = colSums(abs(weights)) != 0
   number_active_nodes = sum(colSums(abs(weights)) != 0)
   active_weights = weights[, active_index]
   list_terms = list()
   for (m in 1 : number_active_nodes)
   {
      list_terms[[m]] = as.numeric(which(active_weights[, m] != 0))
   }
   list_terms = unique(list_terms)
   return(list_terms)
}

# caluclate B-spline acctivation function for a pre-specified order
bspline_general = function(h, order)
{
   knots = quantile(c(-1, 1), probs = seq(0, 1, length = order + 1))
   # constant: order = 1
   # linear: order = 2
   # ....
   bs = bspline_zero(h, knots, order, 1)
   bs = bs / bspline_zero(0, knots, order, 1)
      
   return(bs)
}

# caluclate the derivative of B-spline acctivation function for a pre-specified order
bspline_derivative_general = function(h, order)
{
   knots = quantile(c(-1, 1), probs = seq(0, 1, length = order + 1))
   # constant: order = 1
   # linear: order = 2
   # ....
   bs_dev = bspline_derivatives(h, knots, order, 1, 1)
   bs_dev = bs_dev / bspline_zero(0, knots, order, 1)
   return(bs_dev)
}

bspline_mat_general = function(h, order)
{
   knots = quantile(c(-1, 1), probs = seq(0, 1, length = order + 1))
   # constant: order = 1
   # linear: order = 2
   # ....
   bs = bspline_zero(h, knots, order, 1)
   bs = bs / bspline_zero(0, knots, order, 1)
   return(matrix(bs, nrow(h), ncol(h)))
}

bspline_derivative_mat_general = function(h, order)
{
   knots = quantile(c(-1, 1), probs = seq(0, 1, length = order + 1))
   # constant: order = 1
   # linear: order = 2
   # ....
   bs_dev = bspline_derivatives(h, knots, order, 1, 1)
   bs_dev = bs_dev / bspline_zero(0, knots, order, 1)
   return(matrix(bs_dev, nrow(h), ncol(h)))
}

# calculate bspline function when derivative equals to zero
bspline_zero = function(x, knots, order, j)
{
   if (order > 1)
   {
      if (knots[j + order - 1] > knots[j])
         a = (x - knots[j]) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (knots[j + order] - x) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline_zero(x, knots, order - 1, j) +
                b * bspline_zero(x, knots, order - 1, j + 1))
   }
   else
   {
      bspline = rep(0, length(x))
      bspline[knots[j] <= x & x < knots[j + 1]] = 1
      return(bspline)
   }
}

# calculate the derivative of bspline function when derivative equals to zero
bspline_derivatives = function(x, knots, order, derivative = 1, j = 1)
{
   if (derivative > 0)
   {
      if (knots[j + order - 1] > knots[j])
         a = (order - 1) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (order - 1) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline_derivatives(x, knots, order - 1, derivative - 1, j) -
                b * bspline_derivatives(x, knots, order - 1, derivative - 1, j + 1))
   }
   else
      return(bspline_zero(x, knots, order, j))
}

# FT enumerates integers [from, to] only when from <= to.
FT = function(from, to)
{
   if (from > to)
      ft = NULL
   else
      ft = from : to
   return(ft)
}

# univariate_lasso computes the minimizer zstar of
# q(z) = (2n)^-1 * sum((y - z * x)^2) + lambda * abs(z)
univariate_lasso_rss = function(y, x, lambda, threshold)
{
   a = mean(x * x)
   # if |a| is small, then q(z) = (2n)^-1 * sum(y^2) + lambda * abs(z)
   # is minimized by zero, i.e. zstar = 0
   if (a < threshold)
      zstar = 0
   else
   {
      # if |a| is not so small, the we soft threshold.
      b = mean(x * y)
      zstar = soft_thresholding(b / a, lambda / a)
      if (abs(zstar) < threshold)
         zstar = 0
   }
   return(zstar)
}

# minimizes 0.5 * (z - b)^2 + lambda abs(z)   
soft_thresholding = function(b, lambda)
{
   if (b > lambda)
      return(b - lambda)
   else if (b < -lambda)
      return(b + lambda)
   else return(0)
}

# an indicator function
I = function(x, c1, c2)
{
   x[x < c1] = 0
   x[x > c2] = 0
   return(x)
}

# a toy function for the additive example
g = function(x, a, b, c)
{
   gx = c * (sin(a * pi * (x - b) - pi / 2) + 1) * I(x, b, b + 2/a)
   return(gx)
}