# wrapper for regression 
reduce_dimensions <- function(method, train_data, test_data, ranks) {
  out <- method(train_data[[1]], train_data[[2]], ranks[1], ranks[2], return_scores=TRUE)
  V_joint <- ginv(cbind(train_data[[1]], train_data[[2]])) %*% out$joint
  V_indiv1 <- ginv(train_data[[1]]) %*% out$indiv1
  V_indiv2 <- ginv(train_data[[2]]) %*% out$indiv2
  
  return (list(
    'train_joint' = out$joint,
    'train_indiv1' = out$indiv1,
    'train_indiv2' = out$indiv2,
    'test_joint' = cbind(test_data[[1]], test_data[[2]]) %*% V_joint,
    'test_indiv1' = test_data[[1]] %*% V_indiv1,
    'test_indiv2' = test_data[[2]] %*% V_indiv2
  ))
}