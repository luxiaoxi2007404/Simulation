#---------------------------------------------------------------------
#-----------truncated数据 和 censored 数据对估计一致性的影响----------
#---------------------------------------------------------------------


#数据准备
n = 200
lnw_mean = 2.75
lnw_sd = 0.6
epsilon_mean = 0
epsilon_sd = 1000
beta1 = -2500
beta2 = 1000

lnw = rnorm(n, lnw_mean, lnw_sd)
epsilon = rnorm(n,epsilon_mean,epsilon_sd)
y_star = beta1 + beta2*lnw + epsilon
reg_data = cbind(y_star,lnw)

#数据描述
head(reg_data)
summary(reg_data)
hist(reg_data[,1],20,prob=T,main="Hist and fitted density of reg_data",xlab = "y_star")
lines(density(reg_data[,1]), col="red")

censored_data = reg_data
censored_data[,1][y_star < 0] = 0
hist(censored_data[,1],20,main="Hist of censored_data",xlab = "y_star")

truncated_data = reg_data[y_star > 0,]
hist(truncated_data[,1],20,main="Hist of truncated_data",xlab = "y_star")

#------------------比较截尾、删失、正常数据的条件均值--------------

#Latent variable conditional mean
uncensored_mean = function(lnw) {beta1 + beta2*lnw}
#left-truncated at 0
truncated_mean = function(lnw) {beta1 + beta2*lnw + epsilon_sd*(dnorm((beta1 + beta2*lnw)/epsilon_sd)/pnorm((beta1 + beta2*lnw)/epsilon_sd))}
#left-censored at 0
censored_mean = function(lnw) {pnorm((beta1 + beta2*lnw)/epsilon_sd) * (beta1 + beta2*lnw) + epsilon_sd*(dnorm((beta1 + beta2*lnw)/epsilon_sd))}

plot(lnw,y_star,main="Censored and Truncated Means",ylab = "Different Conditional Means")
curve(uncensored_mean,min(lnw),max(lnw),lty = 1,add = TRUE)
curve(truncated_mean,min(lnw),max(lnw),lty = 2,add = TRUE)
curve(censored_mean,min(lnw),max(lnw),lty = 3,add = TRUE)
legend("bottomright",legend=c("latent variable","uncensored_mean","truncated_mean","censored_mean"),cex = 0.5,lty = c(0,1,2,3))

#===========================探究估计方法对估计结果的影响========================

#--------Simulation考察OLS方法估计截尾和截断数据的稳健性------------

#对未删失数据、截尾、截断数据集直接做做OLS估计，并分别计算MSE
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rnorm(n,epsilon_mean,epsilon_sd)
  y_star = beta1 + beta2*lnw + epsilon
  
  reg_data = cbind(y_star,lnw)
  original_reg = lm(y_star~lnw)
  coef[i,1] = original_reg$coefficients[2]
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  censored_reg = lm(censored_data[,1]~lnw)
  coef[i,2] = censored_reg$coefficients[2]
  
  truncated_data = reg_data[y_star > 0,]
  truncated_reg = lm(truncated_data[,1]~truncated_data[,2])
  coef[i,3] = truncated_reg$coefficients[2]
}

# 求beta2的MSE（真值为1000）
methods_summary = matrix(data = NA,4,3)
rownames(methods_summary) = c("var","bias","mse","bias<0.4*mse")
colnames(methods_summary) = c("original_reg","censored_reg","truncated_reg")

methods_summary[1,] = apply(coef,2,var)
bias = function(i){mean(i)-beta2}
methods_summary[2,] = apply(coef,2,bias)
methods_summary[3,] =methods_summary[1,] + methods_summary[2,]^2
methods_summary[4,] = ifelse(methods_summary[2,] < 0.4*sqrt(methods_summary[1,]),"TRUE","FALSE")

data.frame(methods_summary)


#------------- Simulation考察MLE方法估计截尾数据的稳健性------------

censored_MLE = matrix(data = NA,sim,1) 
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rnorm(n,epsilon_mean,epsilon_sd)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  lnL = function(beta) -sum(d*((-1/2)*log(2*pi)-1/2*log(beta[1]^2)-1/(2*beta[1]^2)*(censored_data[,1]-beta[2]-censored_data[,2]*beta[3])^2) + (1-d)*log(1-pnorm((beta[2]+censored_data[,2]*beta[3])/beta[1])))
  lnL(c(1000,-1000,1000))
  censored_MLE[i,1] = nlm(lnL,(c(1000,-1000,1000)))$estimate[3]
}
censored_var = var(censored_MLE[,1]) 
censored_bias = mean(censored_MLE[,1]) - beta2 
censored_mse = censored_var + censored_bias^2 
ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE")

censored_summary = c(censored_var,censored_bias,censored_mse,ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE")) 
censored_summary 


#------------- Simulation考察MLE方法估计截断数据的稳健性------------

truncated_MLE = matrix(data = NA,sim,1) 
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rnorm(n,epsilon_mean,epsilon_sd)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  truncated_data = reg_data[y_star > 0,]
  
  lnL = function(beta) -sum(-1/2*log(beta[1]^2)-1/2*log(2*pi)-1/(2*beta[1]^2)*(truncated_data[,1]-beta[2]-truncated_data[,2]*beta[3])^2 -log(pnorm((beta[2]+truncated_data[,2]*beta[3])/beta[1])))
  lnL(c(-1000,1000,1000))
  truncated_MLE[i,1] = nlm(lnL,(c(-1000,1000,1000)))$estimate[3]
}
truncated_var = var(truncated_MLE[,1])
truncated_bias = mean(truncated_MLE[,1]) - beta2
truncated_mse = truncated_var + truncated_bias^2
ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE")

truncated_summary = c(truncated_var,truncated_bias,truncated_mse,ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE"))
truncated_summary

#---------- Simulation考察Heckman Two-Step方法估计截断数据的稳健性---------
#---------------------（需要用到删失数据的信息）--------------------------- 
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){ 
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rnorm(n,epsilon_mean,epsilon_sd)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  step_one = cbind(d,lnw)
  model1 = glm(d~lnw,family = binomial(link = probit))
  
  alpha = summary(model1)$coefficients[2,1]
  coef[i,1] = alpha
  inverse_Mills_ratio = dnorm(lnw * alpha)/pnorm(lnw * alpha)
  reg_data = cbind(reg_data,inverse_Mills_ratio)
  
  truncated_data = reg_data[y_star > 0,]
  step_two = lm(truncated_data[,1]~truncated_data[,2] + truncated_data[,3])
  coef[i,2] = summary(step_two)$coefficients[2,1]
  coef[i,3] = summary(step_two)$coefficients[3,1]
}
Heckman_truncate_var = var(coef[,2])
Heckman_truncate_bias = mean(coef[,2]) - beta2
Heckman_truncate_mse = Heckman_truncate_var + Heckman_truncate_bias^2
ifelse(Heckman_truncate_bias < 0.4*sqrt(Heckman_truncate_var),"TRUE","FALSE")

truncate_summary = c(Heckman_truncate_var,Heckman_truncate_bias,Heckman_truncate_mse,ifelse(Heckman_truncate_bias < 0.4*sqrt(Heckman_truncate_var),"TRUE","FALSE"))
truncate_summary

#========================探究epsilon的分布对估计结果的影响======================

# ----------------------假设epsilon的分布非对称------------------

#准备扰动项服从卡方分布的数据
sigma2 = 1000
df = 2
chisq = rchisq(n,df)
epsilon = (chisq - df)/(2*df)*sigma2^2

lnw = rnorm(n, lnw_mean, lnw_sd)
y_star = beta1 + beta2*lnw + epsilon
reg_data = cbind(y_star,lnw)

#数据描述
head(reg_data)
summary(reg_data)
hist(reg_data[,1],20,prob=T,main="Hist and fitted density of reg_data",xlab = "y_star")
lines(density(reg_data[,1]), col="red")

censored_data = reg_data
censored_data[,1][y_star < 0] = 0
hist(censored_data[,1],20,main="Hist of censored_data",xlab = "y_star")

truncated_data = reg_data[y_star > 0,]
hist(truncated_data[,1],20,main="Hist of truncated_data",xlab = "y_star")

#-----Simulation考察OLS方法估计截尾和截断数据的稳健性------

#不同的数据集直接做OLS估计，并分别计算MSE
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  chisq = rchisq(n,df)
  epsilon = (chisq - df)/(2*df)*sigma2^2
  y_star = beta1 + beta2*lnw + epsilon
  
  reg_data = cbind(y_star,lnw)
  original_reg = lm(y_star~lnw)
  coef[i,1] = original_reg$coefficients[2]
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  censored_reg = lm(censored_data[,1]~lnw)
  coef[i,2] = censored_reg$coefficients[2]
  
  truncated_data = reg_data[y_star > 0,]
  truncated_reg = lm(truncated_data[,1]~truncated_data[,2])
  coef[i,3] = truncated_reg$coefficients[2]
}

# 求beta2的MSE（真值为1000）
methods_summary = matrix(data = NA,4,3)
rownames(methods_summary) = c("var","bias","mse","bias<0.4*mse")
colnames(methods_summary) = c("original_reg","censored_reg","truncated_reg")

methods_summary[1,] = apply(coef,2,var)
bias = function(i){mean(i)-beta2}
methods_summary[2,] = apply(coef,2,bias)
methods_summary[3,] =methods_summary[1,] + methods_summary[2,]^2
methods_summary[4,] = ifelse(methods_summary[2,] < 0.4*sqrt(methods_summary[1,]),"TRUE","FALSE")

methods_summary

#----- Simulation考察MLE方法估计截尾数据的稳健性-----

censored_MLE = matrix(data = NA,sim,1)
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  chisq = rchisq(n,df)
  epsilon = (chisq - df)/(2*df)*sigma2^2
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  lnL = function(beta) -sum(d*((-1/2)*log(2*pi)-1/2*log(beta[1]^2)-1/(2*beta[1]^2)*(censored_data[,1]-beta[2]-censored_data[,2]*beta[3])^2) + (1-d)*log(1-pnorm((beta[2]+censored_data[,2]*beta[3])/beta[1])))
  lnL(c(-1000,1000,1000))
  censored_MLE[i,1] = nlm(lnL,(c(-1000,1000,1000)))$estimate[3]
}
censored_var = var(censored_MLE[,1])
censored_bias = mean(censored_MLE[,1]) - beta2
censored_mse = censored_var + censored_bias^2
ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE")

censored_summary = c(censored_var,censored_bias,censored_mse,ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE"))
censored_summary

#----- Simulation考察MLE方法估计截断数据的稳健性-----

truncated_MLE = matrix(data = NA,sim,1)
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  chisq = rchisq(n,df)
  epsilon = (chisq - df)/(2*df)*sigma2^2
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  truncated_data = reg_data[y_star > 0,]
  
  lnL = function(beta) -sum(-1/2*log(beta[1]^2)-1/2*log(2*pi)-1/(2*beta[1]^2)*(truncated_data[,1]-beta[2]-truncated_data[,2]*beta[3])^2 -log(pnorm((beta[2]+truncated_data[,2]*beta[3])/beta[1])))
  lnL(c(-1000,1000,1000))
  truncated_MLE[i,1] = nlm(lnL,(c(-1000,1000,1000)))$estimate[3]
}
truncated_var = var(truncated_MLE[,1])
truncated_bias = mean(truncated_MLE[,1]) - beta2
truncated_mse = truncated_var + truncated_bias^2
ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE")

truncated_summary = c(truncated_var,truncated_bias,truncated_mse,ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE"))
truncated_summary

#---Simulation考察Heckman Two-Step方法估计截断数据的稳健性---------
#-------------（需要用到删失数据的信息）--------------------------- 
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){ 
  lnw = rnorm(n, lnw_mean, lnw_sd)
  chisq = rchisq(n,df)
  epsilon = (chisq - df)/(2*df)*sigma2^2
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  step_one = cbind(d,lnw)
  model1 = glm(d~lnw,family = binomial(link = probit))
  
  alpha = summary(model1)$coefficients[2,1]
  coef[i,1] = alpha
  inverse_Mills_ratio = dnorm(lnw * alpha)/pnorm(lnw * alpha)
  reg_data = cbind(reg_data,inverse_Mills_ratio)
  
  truncated_data = reg_data[y_star > 0,]
  step_two = lm(truncated_data[,1] ~ truncated_data[,2] + truncated_data[,3])
  coef[i,2] = summary(step_two)$coefficients[2,1]
  #coef[i,3] = summary(step_two)$coefficients[3,1]
}
Heckman_truncated_var = var(coef[,2])
Heckman_truncated_bias = mean(coef[,2]) - beta2
Heckman_truncated_mse = Heckman_truncated_var + Heckman_truncated_bias^2
ifelse(Heckman_truncated_bias < 0.4*sqrt(Heckman_truncated_var),"TRUE","FALSE")

truncated_summary = c(Heckman_truncated_var,Heckman_truncated_bias,Heckman_truncated_mse,ifelse(Heckman_truncated_bias < 0.4*sqrt(Heckman_truncated_var),"TRUE","FALSE"))
truncated_summary

# ----------------------假设epsilon的分布厚尾------------------

#准备扰动项服从t分布的数据
df = 2
epsilon = rt(n,df)
lnw = rnorm(n, lnw_mean, lnw_sd)
y_star = beta1 + beta2*lnw + epsilon
reg_data = cbind(y_star,lnw)

#数据描述
head(reg_data)
summary(reg_data)
hist(reg_data[,1],20,prob=T,main="Hist and fitted density of reg_data",xlab = "y_star")
lines(density(reg_data[,1]), col="red")

censored_data = reg_data
censored_data[,1][y_star < 0] = 0
hist(censored_data[,1],20,main="Hist of censored_data",xlab = "y_star")

truncated_data = reg_data[y_star > 0,]
hist(truncated_data[,1],20,main="Hist of truncated_data",xlab = "y_star")

#-----Simulation考察OLS方法估计截尾和截断数据的稳健性-----

#不同的数据集直接做OLS估计，并分别计算MSE
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rt(n,df)
  y_star = beta1 + beta2*lnw + epsilon
  
  reg_data = cbind(y_star,lnw)
  original_reg = lm(y_star~lnw)
  coef[i,1] = original_reg$coefficients[2]
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  censored_reg = lm(censored_data[,1]~lnw)
  coef[i,2] = censored_reg$coefficients[2]
  
  truncated_data = reg_data[y_star > 0,]
  truncated_reg = lm(truncated_data[,1]~truncated_data[,2])
  coef[i,3] = truncated_reg$coefficients[2]
}

# 求beta2的MSE（真值为1000）
methods_summary = matrix(data = NA,4,3)
rownames(methods_summary) = c("var","bias","mse","bias<0.4*mse")
colnames(methods_summary) = c("original_reg","censored_reg","truncated_reg")

methods_summary[1,] = apply(coef,2,var)
bias = function(i){mean(i)-beta2}
methods_summary[2,] = apply(coef,2,bias)
methods_summary[3,] =methods_summary[1,] + methods_summary[2,]^2
methods_summary[4,] = ifelse(methods_summary[2,] < 0.4*sqrt(methods_summary[1,]),"TRUE","FALSE")

methods_summary


#------ Simulation考察MLE方法估计截尾数据的稳健性-------

censored_MLE = matrix(data = NA,sim,1)
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rt(n,df)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  lnL = function(beta) -sum(d*((-1/2)*log(2*pi)-1/2*log(beta[1]^2)-1/(2*beta[1]^2)*(censored_data[,1]-beta[2]-censored_data[,2]*beta[3])^2) + (1-d)*log(1-pnorm((beta[2]+censored_data[,2]*beta[3])/beta[1])))
  lnL(c(-2000,2000,2000))
  censored_MLE[i,1] = nlm(lnL,(c(-2000,2000,2000)))$estimate[3]
}
censored_var = var(censored_MLE[,1])
censored_bias = mean(censored_MLE[,1]) - beta2
censored_mse = censored_var + censored_bias^2
ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE")

censored_summary = c(censored_var,censored_bias,censored_mse,ifelse(censored_bias < 0.4*sqrt(censored_var),"TRUE","FALSE"))
censored_summary

#------ Simulation考察MLE方法估计截断数据的稳健性------

truncated_MLE = matrix(data = NA,sim,1)
for (i in 1:1000){
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rt(n,df)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  truncated_data = reg_data[y_star > 0,]
  
  lnL = function(beta) -sum(-1/2*log(beta[1]^2)-1/2*log(2*pi)-1/(2*beta[1]^2)*(truncated_data[,1]-beta[2]-truncated_data[,2]*beta[3])^2 -log(pnorm((beta[2]+truncated_data[,2]*beta[3])/beta[1])))
  lnL(c(-1000,1000,1000))
  truncated_MLE[i,1] = nlm(lnL,(c(-1000,1000,1000)))$estimate[3]
}
truncated_var = var(truncated_MLE[,1])
truncated_bias = mean(truncated_MLE[,1]) - beta2
truncated_mse = truncated_var + truncated_bias^2
ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE")

truncated_summary = c(truncated_var,truncated_bias,truncated_mse,ifelse(truncated_bias < 0.4*sqrt(truncated_var),"TRUE","FALSE"))
truncated_summary

#---Simulation考察Heckman Two-Step方法估计截断数据的稳健性---------
#-------------（需要用到删失数据的信息）--------------------------- 
coef = matrix(data = NA,sim,3)
sim = 1000

for (i in 1:sim){ 
  lnw = rnorm(n, lnw_mean, lnw_sd)
  epsilon = rt(n,df)
  y_star = beta1 + beta2*lnw + epsilon
  reg_data = cbind(y_star,lnw)
  
  censored_data = reg_data
  censored_data[,1][y_star < 0] = 0
  
  d = ifelse(censored_data[,1] > 0,1,0)
  step_one = cbind(d,lnw)
  model1 = glm(d~lnw,family = binomial(link = probit))
  
  alpha = summary(model1)$coefficients[2,1]
  coef[i,1] = alpha
  inverse_Mills_ratio = dnorm(lnw * alpha)/pnorm(lnw * alpha)
  reg_data = cbind(reg_data,inverse_Mills_ratio)
  
  truncated_data = reg_data[y_star > 0,]
  step_two = lm(truncated_data[,1] ~ truncated_data[,2] + truncated_data[,3])
  coef[i,2] = summary(step_two)$coefficients[2,1]
  #coef[i,3] = summary(step_two)$coefficients[3,1]
}
Heckman_truncated_var = var(coef[,2])
Heckman_truncated_bias = mean(coef[,2]) - beta2
Heckman_truncated_mse = Heckman_truncated_var + Heckman_truncated_bias^2
ifelse(Heckman_truncated_bias < 0.4*sqrt(Heckman_truncated_var),"TRUE","FALSE")

truncated_summary = c(Heckman_truncated_var,Heckman_truncated_bias,Heckman_truncated_mse,ifelse(Heckman_truncated_bias < 0.4*sqrt(Heckman_truncated_var),"TRUE","FALSE"))
truncated_summary