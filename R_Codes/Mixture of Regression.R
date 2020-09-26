# ==================== Preparation ============================#
library(MASS)
library(ggplot2)

#=================== Initialization ====================#
initial<- function(theta_true){
  theta_0<-matrix(rnorm(10,10,10),10,1)
  return(theta_0)
}
#==================== w function =======================#
w <- function(x,y,theta){
  w_theta <- exp(-(y-x%*%theta)^2/(2*sigma^2))/(exp(-(y-x%*%theta)^2/(2*sigma^2))+exp(-(y+x%*%theta)^2/(2*sigma^2)))
  return(w_theta)
}

# ====================outer producnt ===========================#
out <- function(x){
  xx <- x%*%t(x)
  return(xx)
}

# ==================== Update Rules ===========================#
# EM_update
EM_update <- function(x,y,theta){
    
    xx <- rowSums(apply(x, 1, out))
    xx <- matrix(xx,d,d)
    xx <- solve(xx)
    w_yx <- matrix(0,d,1)
    w_y <- NA
    for(i in 1:n){
      w_y<-(2*w(x[i,],y[i],theta)-1)*y[i]
      w_yx <- w_yx+rep(w_y,d)*x[i,]
    }
   Mn <- xx%*%w_yx
   return(Mn)
}

# ============= EM estimator =================#
em_est <- function(x,y){
  
  # Initialization
  theta <- list()
  theta[[1]] <- initial(theta_true)
  
  for (i in c(2:20)){
    theta[[i]] <- EM_update(x,y,theta[[i-1]])
  }
  return(theta)
}

# ================ gradient EM =======================#
 
GEM <- function(x,y,theta,alpha){
  w_yx <- matrix(0,d,1)
  w_y<-NA
  for(i in 1:n){
    w_y<-(2*w(x[i,],y[i],theta)-1)*y[i]
    w_yx <- w_yx+rep(w_y,d)*x[i,]-x[i,]%*%t(x[i,])%*%theta
  }
  gem<-theta+(alpha/n)*(w_yx)
}

# ============= gradient EM estimator =================#
gem_est <- function(x,y,alpha,iter= 20){
  # Initialization
  theta <- list()
  theta[[1]] <- initial(theta_true)

  for (i in c(2:iter)){
    theta[[i]] <- GEM(x,y,theta[[i-1]],alpha)
  }
  return(theta)
}
#================= log error =================#
# log_sta_error
log_sta_error <- function(result,theta_true){
  log_sta_error<-c()
  for (i in 1:length(result)) {
    log_sta_error[i] <- log(sqrt(sum(((result[[i]])-theta_true)^2)))
  }
  return(log_sta_error)
}
# log_opt_error
log_opt_error <- function(result){
  log_opt_error<-c()
  for (i in 1:(length(result)-1)) {
    log_opt_error[i] <- log(sqrt(sum(((result[[i]])-result[[length(result)]])^2)))
  }
  return(log_opt_error)
}
# ===================== Simulations EM ============================#

#Basic parameter
# dimension:  d = 10
d <- 10
# sample size: n = 1000
n <- 1000
# theta_true
theta_true <- matrix(seq(1,28,3),d,1)
# Signal-to-noise ratio: ||theta|| / sigma = 2
sigma <- sqrt((sum((theta_true)^2)))/2


simu_EM <-function(simulation,iter=20){
  ##log value
  sta <- matrix(NA,simulation,iter)
  opt <- matrix(NA,simulation,c(iter-1))
  ## 10 times simulation
  for (index in 1:simulation ) {
    # X
    x <- mvrnorm(n,rep(0,d),diag(1,d))
    # epsilon
    epsilon <- rnorm(n,0,sigma)
    # Y = X%*%theta_true + epsilon
    y <- x%*%theta_true + epsilon
    #result
    result <- em_est(x,y)
    ##log-value
    sta[index,] <- log_sta_error(result,theta_true)
    opt[index,] <- log_opt_error(result)
   
  }
  results <- list(opt=t(opt),sta=t(sta))
  return(results)
}

#=================== plot EM simulations ========================#
result <- simu_EM(simulation = 10)

sta_data <- data.frame(count=rep(1:20, d), value=as.vector(result[["sta"]]),
                 variable=rep(paste0("s_category", 1:d), each=20))
opt_data <- data.frame(count=rep(1:19, d), value=as.vector(result[["opt"]]),
                       variable=rep(paste0("o_category", 1:d), each=19))

# plot
ggplot()+ 
  geom_line(data = sta_data, aes(x=count, y=value,group=variable),color="red")+
  geom_line(data = opt_data, aes(x=count, y=value,group=variable),color="blue")+
  ggtitle('EM, Mixtures of Regression')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") 


# ===================== Simulations GEM ============================#

simu_GEM <-function(simulation,iter=20){
  ##log value
  sta <- matrix(NA,simulation,iter)
  opt <- matrix(NA,simulation,c(iter-1))
  ## 10 times simulation
  for (index in 1:simulation ) {
    # X
    x <- mvrnorm(n,rep(0,d),diag(1,d))
    # epsilon
    epsilon <- rnorm(n,0,sigma)
    # Y = X%*%theta_true + epsilon
    y <- x%*%theta_true + epsilon
    #result
    result <- gem_est(x,y,alpha=0.5)
    ##log-value
    sta[index,] <- log_sta_error(result,theta_true)
    opt[index,] <- log_opt_error(result)
    
  }
  results <- list(opt=t(opt),sta=t(sta))
  return(results)
}


#=================== plot EM simulations ========================#
# result & plot
result_g <- simu_GEM(simulation = 10)
# data frame 
sta_data <- data.frame(count=rep(1:20, d), value=as.vector(result_g[["sta"]]),
                       variable=rep(paste0("s_category", 1:d), each=20))
opt_data <- data.frame(count=rep(1:19, d), value=as.vector(result_g[["opt"]]),
                       variable=rep(paste0("o_category", 1:d), each=19))

# plot
ggplot()+ 
  geom_line(data = sta_data, aes(x=count, y=value,group=variable),color="red")+
  geom_line(data = opt_data, aes(x=count, y=value,group=variable),color="blue")+
  ggtitle('GEM, Mixtures of Regression')+
  labs(x = 'Iteration Count', y ='Log Error',color = "Legend") 

# ==================== Different SNR ===================== #

SNR<-c(0.5,0.75,1,1.8,2.5)
v_SNR<-list()
for (t in SNR){
sigma <- sqrt((sum((theta_true)^2)))/t
v_SNR[[round(4*t)]]<-simu_EM(10)[[1]]
}
v_SNR<-v_SNR[!sapply(v_SNR, is.null)]

# data frame 
opt_SNR<-list()
for(i in 1:5){
opt_SNR [[i]]<- data.frame(count=rep(1:19, d), value=as.vector(v_SNR[[i]]),
                       variable=rep(paste0("o_category", 1:d), each=19))
}
# plot
ggplot()+ 
  geom_line(data = opt_SNR[[1]], aes(x=count, y=value,group=variable,color="SNR=0.5"))+
  geom_line(data = opt_SNR[[2]], aes(x=count, y=value,group=variable,color="SNR=0.75"))+
  geom_line(data = opt_SNR[[3]], aes(x=count, y=value,group=variable,color="SNR=1"))+
  geom_line(data = opt_SNR[[4]], aes(x=count, y=value,group=variable,col="SNR=1.8"))+
  geom_line(data = opt_SNR[[5]], aes(x=count, y=value,group=variable,col="SNR=2.5"))+
  ggtitle('EM, Mixtures of Regression')+
  labs(x = 'Iteration Count', y ='Log Error')+
  theme_light()


