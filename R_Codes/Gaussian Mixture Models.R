# ==================== Preparation ============================#
library(MASS)
library(ggplot2)
# ==================== Update Rules ===========================#

# sigma
sigma_sim <- function(theta_ture){
  return(norm(theta_ture, type = "2")/2)
}

# generate sample y
y_generate <- function(theta_true, sigma_true, d, n){
  y_sim <- matrix(0, n, d)
  z <- sample(c(0,1), n, prob = c(0.5, 0.5), replace = T)
  for (i in 1:n){
    if (z[i] == 0){
      y_sim[i,] <- mvrnorm(1, -theta_true, diag(sigma_true,d))
    } else {
      y_sim[i,] <- mvrnorm(1, theta_true, diag(sigma_true,d))
    }
  }
  return(y_sim)
}

# function w
w_y <- function(theta_current, y, sigma_true){
  a <- exp(-norm(theta_current - y, type = "2")^2/(2*sigma_true^2)) /
    (exp(-norm(theta_current - y, type = "2")^2/(2*sigma_true^2)) + 
       exp(-norm(theta_current + y, type = "2")^2/(2*sigma_true^2)))
  return(a)
}

# standard EM algorithm
EM <- function(theta_true, theta_initial, sigma_true, iter, d, n){
  y <- y_generate(theta_true, sigma_true, d ,n)
  theta_res <- matrix(0, iter+1, d)
  theta_res[1,] <- theta_initial
  for (i in 2:(iter+1)){
    for (j in 1:n){
      theta_res[i,] <- theta_res[i,] + (2*w_y(theta_res[i-1,], y[j,], sigma_true)-1)*y[j,]/n
    }
  }
  return(theta_res)
}

# Gradient EM algorithm 
GEM <- function(theta_true, theta_initial, sigma_true, iter, d, n, lrate){
  y <- y_generate(theta_true, sigma_true, d ,n)
  theta_res <- matrix(0, iter+1, d)
  theta_res[1,] <- theta_initial
  for (i in 2:(iter+1)){
    for (j in 1:n){
      theta_res[i,] <- theta_res[i,] + (2*w_y(theta_res[i-1,], y[j,], sigma_true)-1)*y[j,]/n
    }
    theta_res[i,] <- theta_res[i-1,] + lrate*(theta_res[i,] - theta_res[i-1,])
  }
  return(theta_res)
}

# stat error 
stat_error <- function(theta_res, theta_true, iter){
  stat_errorlist <- c()
  for (i in 1:(iter+1)){
    stat_errorlist <- c(stat_errorlist, log(norm(theta_res[i,] - theta_true, type = "2")))
  }
  return(stat_errorlist)
}

# opt error 
opt_error <- function(theta_res, iter){
  opt_errorlist <- c()
  for (i in 1:(iter+1)){
    opt_errorlist <- c(opt_errorlist, log(norm(theta_res[i,] - theta_res[iter+1,], type = "2")))
  }
  return(opt_errorlist)
}


# ======== Mixture of Gaussian ==========#
# obtain true sigma and true theta
theta_true_matrix <- matrix(0, 10, 10)
sigma_true_vec <- rep(0,10)
for (i in 1:10){
  theta_true_matrix[i,] <- rnorm(10)
  sigma_true_vec[i] <- sigma_sim(theta_true_matrix[i,])
}

# ======== Simulation ==========#
# stat error and opt error for EM 
stat_error_matrix <- matrix(0, 10, 21)
opt_error_matrix <- matrix(0, 10, 21)

for (k in 1:10){
  # initial theta: make sure it is in the Ball(r,theta)
  initial <- rnorm(10)+theta_true_matrix[k,]
  stat_error_matrix[k,] <- stat_error(EM(theta_true_matrix[k,], initial, sigma_true_vec[k], 20, 10, 1000), theta_true_matrix[k,], 20)
  opt_error_matrix[k,] <- opt_error(EM(theta_true_matrix[k,], initial, sigma_true_vec[k], 20, 10, 1000), 20)
}

#plot 
ggplot()+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[1,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[2,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[3,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[4,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[5,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[6,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[7,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[8,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[9,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix[10,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[1,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[2,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[3,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[4,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[5,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[6,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[7,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[8,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[9,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix[10,]),color = 'blue')+
  ggtitle('EM, Mixture of Gaussians')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") +
  theme_light()


# stat error and opt error for Gradient EM 
stat_error_matrix2 <- matrix(0, 10, 21)
opt_error_matrix2 <- matrix(0, 10, 21)

for (k in 1:10){
  # initial theta: make sure it is in the Ball(r,theta)
  initial <- rnorm(10)+theta_true_matrix[k,]
  stat_error_matrix2[k,] <- stat_error(GEM(theta_true_matrix[k,], initial, sigma_true_vec[k], 20, 10, 1000, 0.5), theta_true_matrix[k,], 20)
  opt_error_matrix2[k,] <- opt_error(GEM(theta_true_matrix[k,], initial, sigma_true_vec[k], 20, 10, 1000, 0.5), 20)
}

ggplot()+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[1,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[2,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[3,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[4,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[5,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[6,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[7,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[8,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[9,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=stat_error_matrix2[10,]),color = 'red')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[1,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[2,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[3,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[4,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[5,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[6,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[7,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[8,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[9,]),color = 'blue')+
  geom_line(aes(x=c(1:21),y=opt_error_matrix2[10,]),color = 'blue')+
  ggtitle('Gradient EM, Mixture of Gaussians')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") +
  theme_light()




