
# functions for the two-way interaction example
g1 = function(z)
   return(z)
g2 = function(z)
   return((2 * z - 1)^2)
g3 = function(z)
   return(sin(2 * pi * z) / (2 - sin(2 * pi * z)))
g4 = function(z)
   return(0.1 * sin(2 * pi * z) + 0.2 * cos(2 * pi * z) + 0.3 * sin(2 * pi * z)^2 + 0.4 * cos(2 * pi * z)^3 + 0.5 * sin(2 * pi * z)^3)
