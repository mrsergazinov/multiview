my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(ggplot2)
source('src/bounds.R')

# plot density of product of two random variables for different values of q1 and q2
q1 <- seq(0.1, 0.5, by = 0.2)
q2 <- seq(0.1, 0.5, by = 0.2)
data <- data.frame(lambda = numeric(0), 
                   f_lambda = numeric(0),
                   q = numeric(0))

for (i in 1:length(q1)) {
  for (j in 1:length(q2)) {
    lambda <- seq(0.01, 0.99, by = 0.01)
    values <- f_lambda(lambda, q1[i], q2[j])
    data <- rbind(data, data.frame(lambda = lambda, 
                                   f_lambda = values,
                                   q = paste0("q1 = ", q1[i], ", q2 = ", q2[j])))
  }
}
ggplot(data, aes(x = lambda, y = f_lambda, color = q)) + 
  geom_line() + 
  theme_minimal() +
  xlab('Lambda') + 
  ylab('Density')


