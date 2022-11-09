forward_algorithm = function(observations,time,init = pi,tpm=Gamma.matrix,parameters=lambda){
  
  # This is for time t = 1
  alpha = init %*% diag(dpois(observations[1],parameters))
  
  
  #Alternative way of computing alpha for t = 1
  #alpha = init * dpois(observations[1],parameters)
  
  if(time > 1){
    for (i in 2:time) {
      
      alpha = alpha %*% tpm%*%diag(dpois(observations[i],parameters))
      
      #Alternative way of computing alpha for t > 1
      #alpha = alpha %*% tpm*dpois(observations[i],parameters)
    }
  }
  return(alpha)
}
