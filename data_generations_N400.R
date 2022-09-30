##################################################################
# load
##################################################################
rm(list = ls())
##################################################################
# data generation
##################################################################
source("sub_code.R")
set.seed(22)
sample_size = 400
number_predictors = 4
data = list()
##################################################################
# data generations
number_sim = 100
data = list()
for (s in 1 : number_sim)
{
   x = matrix(0, sample_size, number_predictors)
   for (j in 1 : number_predictors)
      x[, j] = runif(sample_size, 0, 1)
   f1 = g1(x[, 1])
   f2 = g2(x[, 2])
   f3 = g3(x[, 3])
   f4 = g4(x[, 4])
   f34 = g1(x[, 3] * x[, 4])
   f13 = g2((x[, 1] + x[, 3]) / 2)
   f12 = g3(x[, 1] * x[, 2])
   fs = cbind(f1, f2, f3, f4, f34, f13, f12)
   # f = 4 * sin(x[, 1]) + 4 * x[, 2] * x[, 3] * x[, 4]#f1 + f2 + f3 + f4
   f = f1 + f2 + f3 + f4 + f34 + f13 + f12
   y = f + rnorm(sample_size, sd = 0.2546) # snr 3
   data[[s]] = list()
   data[[s]]$x = x
   data[[s]]$y = y
   data[[s]]$f = f
}
save(file = "interaction_data_N400", data)

