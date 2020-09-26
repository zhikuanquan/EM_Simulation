# ==================== Preparation ============================#
library(MASS)
library(ggplot2)
# ==================== Update Rules ===========================#
# Standard EM_update
EM_update <- function(x,y,theta,sigma){
  y_mu_i <- list()
  sigma_theta_i <- list()
  n <- nrow(y)
  for (i in 1:n) {
    # z_obs_i
    y_i <- y[i,]
    na_index <- is.na(x)[i,]
    x_obs <- x[i,!na_index]
    z_obs_i <- c(x_obs,y_i)
    # U_theta_i
    theta_obs <- theta[!na_index]
    theta_mis <- theta[na_index]
    U_theta_i <- 
      1/(sum((theta_mis)^2)+sigma^2)*cbind(-theta_mis%*%t(theta_obs),theta_mis)
    # mu_theta_i
    mu_theta_i <- c(U_theta_i%*%z_obs_i,x_obs)
    # sigma_theta_i
    sigma_theta_i_diag <- U_theta_i%*%z_obs_i%*%t(x_obs)
    sigma_theta_i_2 <- x_obs%*%t(x_obs)
    sigma_theta_i_1 <- diag(1,nrow(sigma_theta_i_diag))
    sigma_theta_i[[i]] <- rbind(cbind(sigma_theta_i_1,sigma_theta_i_diag),
                                cbind(t(sigma_theta_i_diag),sigma_theta_i_2))
    y_mu_i[[i]] <- y_i*mu_theta_i
  }
  theta_new <- solve(Reduce('+', sigma_theta_i))%*%Reduce('+', y_mu_i)
  return(theta_new)
}
# Gradient EM_update with learning rate alpha
GEM_update <- function(x,y,theta,sigma,alpha=0.5){
  y_mu_i <- list()
  sigma_theta_i <- list()
  output <- list()
  n <- nrow(y)
  for (i in 1:n) {
    # z_obs_i
    y_i <- y[i,]
    na_index <- is.na(x)[i,]
    x_obs <- x[i,!na_index]
    z_obs_i <- c(x_obs,y_i)
    # U_theta_i
    theta_obs <- theta[!na_index]
    theta_mis <- theta[na_index]
    U_theta_i <- 
      1/(sum((theta_mis)^2)+sigma^2)*cbind(-theta_mis%*%t(theta_obs),theta_mis)
    # mu_theta_i
    mu_theta_i <- c(U_theta_i%*%z_obs_i,x_obs)
    # sigma_theta_i
    sigma_theta_i_diag <- U_theta_i%*%z_obs_i%*%t(x_obs)
    sigma_theta_i_2 <- x_obs%*%t(x_obs)
    sigma_theta_i_1 <- diag(1,nrow(sigma_theta_i_diag))
    sigma_theta_i[[i]] <- rbind(cbind(sigma_theta_i_1,sigma_theta_i_diag),
                                cbind(t(sigma_theta_i_diag),sigma_theta_i_2))
    y_mu_i[[i]] <- y_i*mu_theta_i
    output[[i]] <- y_mu_i[[i]]-sigma_theta_i[[i]]%*%theta
  }
  theta_new <- theta + alpha*Reduce('+', output)/n
  return(theta_new)
}
# log_sta_error
log_sta_error <- function(result,theta_true){
  log_sta_error<-c()
  for (i in 1:(length(result)-1)) {
    log_sta_error[i] <- log(sqrt(sum(((result[[i]])-theta_true)^2)))
  }
  return(log_sta_error)
}
# log_opt_error
log_opt_error <- function(result,theta_true){
  log_opt_error<-c()
  for (i in 1:(length(result)-1)) {
    log_opt_error[i] <- log(sqrt(sum(((result[[i]])-result[[(length(result)-1)]])^2)))
  }
  return(log_opt_error)
}
# ======== Linear regression with missing covariates ==========#
mislm <- function(x,y,theta_true,iter= 20,SNR = 2){
  # basic parameter
  d <- ncol(x)
  n <- nrow(y)
  # Initialization
  theta <- list()
  theta[[1]] <- rnorm(d)
  sigma <- sqrt((sum((theta_true)^2)))/SNR
  for (i in c(2:iter)){
    theta[[i]] <- EM_update(x,y,theta[[i-1]],sigma)
    # EM update
    theta[[i+1]] <- theta[i] 
  }
  results <- list(theta=theta)
  return(theta)
}
mislm_g <- function(x,y,theta_true,iter= 20,SNR = 2){
  # basic parameter
  d <- ncol(x)
  n <- nrow(y)
  # Initialization
  theta <- list()
  theta[[1]] <- rnorm(d)
  sigma <- sqrt((sum((theta_true)^2)))/SNR
  for (i in c(2:iter)){
    theta[[i]] <- GEM_update(x,y,theta[[i-1]],sigma)
    # EM update
    theta[[i+1]] <- theta[i] 
  }
  results <- list(theta=theta)
  return(theta)
}

# ===================== Simulations ============================#
# dimension:  d = 10
d <- 10
# sample size: n = 1000
n <- 1000
# missing probability: p = 0.05
p <- 0.05

# Standard EM
mislm_sim <- function(simulation, d=10,n=1000,p=0.05,SNR=2,scale=1){
  # Storing Results
  opt<-matrix(0,simulation,20) #iter = 20
  sta<-matrix(0,simulation,20) #iter = 20
  for (index in 1:simulation) {
    # theta_true
    theta_true <- rnorm(d)*scale
    # Signal-to-noise ratio: ||theta|| / sigma = 2
    sigma <- sqrt((sum((theta_true)^2)))/SNR
    # X
    x <- mvrnorm(n,rep(0,d),diag(1,d))
    # epsilon
    epsilon <- rnorm(n,0,sigma)
    # Y = X%*%theta_true + epsilon
    y <- x%*%theta_true + epsilon
    # randomly missing
    x <- apply(x, 2, function(x){x[sample(c(1:n),floor(n*p))]<-NA; x})
    result <- mislm(x,y,theta_true,SNR=SNR)
    opt[index,] <- log_opt_error(result,theta_true)
    sta[index,] <- log_sta_error(result,theta_true)
  }
  results <- list(opt=t(opt),sta=t(sta))
  return(results)
}
result <- mislm_sim(simulation = 10,SNR = 2)
ggplot()+
  geom_line(aes(x=c(1:20),y=result$sta[,1]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,2]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,3]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,4]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,5]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,6]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,7]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,8]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,9]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$sta[,10]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result$opt[,1]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,2]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,3]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,4]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,5]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,6]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,7]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,8]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,9]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result$opt[,10]),color = 'blue')+
  ggtitle('EM, Missing Data Regression')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") +
  theme_light()

# Different SNR
result_SNR_opt <- list()
result_SNR_sta <- list()
for (i in c(1:5)) {
  result_t <- mislm_sim(simulation = 10,SNR=(i/2))
  result_SNR_opt[[i]] <- result_t$opt
  result_SNR_sta[[i]] <- result_t$sta
}
ggplot()+
  geom_line(aes(x=c(1:20),y=apply(result_SNR_opt[[1]],1,mean)),color = 'blue')+
  geom_line(aes(x=c(1:20),y=apply(result_SNR_opt[[2]],1,mean)),color = 'green')+
  geom_line(aes(x=c(1:20),y=apply(result_SNR_opt[[3]],1,mean)),color = 'yellow')+
  geom_line(aes(x=c(1:20),y=apply(result_SNR_opt[[4]],1,mean)),color = 'red')+
  geom_line(aes(x=c(1:20),y=apply(result_SNR_opt[[5]],1,mean)),color = 'black')+
#  geom_line(aes(x=c(1:20),y=apply(result_SNR_sta[[1]],1,mean)),color = 'blue')+
#  geom_line(aes(x=c(1:20),y=apply(result_SNR_sta[[2]],1,mean)),color = 'green')+
#  geom_line(aes(x=c(1:20),y=apply(result_SNR_sta[[3]],1,mean)),color = 'yellow')+
#  geom_line(aes(x=c(1:20),y=apply(result_SNR_sta[[4]],1,mean)),color = 'red')+
#  geom_line(aes(x=c(1:20),y=apply(result_SNR_sta[[5]],1,mean)),color = 'black')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") +
  theme_light()

# Gradient EM
mislm_sim_g <- function(simulation, d=10,n=1000,p=0.05){
  # Storing Results
  opt<-matrix(0,simulation,20) #iter = 20
  sta<-matrix(0,simulation,20) #iter = 20
  for (index in 1:simulation) {
    # theta_true
    theta_true <- rnorm(d)
    # Signal-to-noise ratio: ||theta|| / sigma = 2
    sigma <- sqrt((sum((theta_true)^2)))/2
    # X
    x <- mvrnorm(n,rep(0,d),diag(1,d))
    # epsilon
    epsilon <- rnorm(n,0,sigma)
    # Y = X%*%theta_true + epsilon
    y <- x%*%theta_true + epsilon
    # randomly missing
    x <- apply(x, 2, function(x){x[sample(c(1:n),floor(n*p))]<-NA; x})
    result <- mislm_g(x,y,theta_true)
    opt[index,] <- log_opt_error(result,theta_true)
    sta[index,] <- log_sta_error(result,theta_true)
  }
  results <- list(opt=t(opt),sta=t(sta))
  return(results)
}
result_g <- mislm_sim_g(simulation = 10)
ggplot()+
  geom_line(aes(x=c(1:20),y=result_g$sta[,1]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,2]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,3]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,4]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,5]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,6]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,7]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,8]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,9]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$sta[,10]),color = 'red')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,1]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,2]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,3]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,4]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,5]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,6]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,7]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,8]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,9]),color = 'blue')+
  geom_line(aes(x=c(1:20),y=result_g$opt[,10]),color = 'blue')+
  xlab('Iteration Count')+
  ylab('Log Error')+
  ggtitle('GEM, Missing Data Regression')+
  theme_light()

