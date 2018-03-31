library(jsonlite)
library(readr)
library(ggplot2)
setwd("~/Desktop")
df_business <- stream_in(file("business_open_Toronto.json"))
#converting more than one attribute dataframe object to 2-d dataframe
df_business <- flatten(df_business)
head(df_business)

#encoding neighbourhood into numeric value
df_business$neighborhood <- as.factor(df_business$neighborhood)
nlevels(df_business$neighborhood)
df_business$neighborhood <- as.numeric(df_business$neighborhood)
str(df_business$neighborhood)
dim(df_business)
str((df_business$stars))

#counts of stars
ggplot(df_business, aes(stars)) + stat_bin()

#variations of stars with different neighbourgood
ggplot(df_business) + 
  geom_boxplot(aes(x = reorder(neighborhood, stars, median), stars, fill=reorder(neighborhood, stars, median)), show.legend = FALSE)

#ggplot(df_business) + geom_boxplot(aes(neighborhood, stars, fill = neighborhood)) + geom_jitter(aes(neighborhood, stars, shape = df_business$neighborhood)) + scale_shape_identity()

df_mean<-tapply(df_business$stars, df_business$neighborhood, mean)
head(df_mean)
#comparing using gibbs sampler 

gibbs_samp <- function(y, ind, maxiter=5000)
{
  #tau_w hyperparameter(values to be modifies according to the data set)
  a0 <- 1/2 ; b0 <- 3
  #tau_b hyperparametr(values to be modifies according to the data set)
  eta0 <- 1/2 ; t0 <- 3
  mu0 <- 50 ; gamma0 <- 1/25
  
  #starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y,ind,mean)
  tau_w <- mean(1/tapply(y,ind,var))
  mu <- mean(theta)
  tau_b<-var(theta)
  n_m<-tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  
  # setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  # MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, taun)
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

fit2 <- gibbs_samp(df_business$stars, df_business$neighborhood)
apply(fit2$params, 2, mean)
apply(fit2$params, 2, sd)
mean(1/sqrt(fit2$params[, 3]))
sd(1/sqrt(fit2$params[, 3]))
theta_hat <- apply(fit2$theta, 2, mean)
ggplot(data.frame(size = tapply(df_business$stars, df_business$neighborhood, length), theta_hat = theta_hat), aes(size, theta_hat)) + geom_point()