------------------------------------------------------------------------
#------------------  -多项选择模型（以钓鱼数据为例）--------------------
------------------------------------------------------------------------
  
  
# 载入数据
library(mnlogit)
data(Fish)
#fix(Fish)

# 数据描述
choices = c('beach','boat','charter','pier')
data_summary = matrix(data = NA,10,4)
for (i in choices){
  Fish_i = Fish[Fish$alt == i,]
  Fish_i_yes = Fish_i[Fish_i$mode == 1,]
  data_summary[10,which(choices==i)] = Observations = dim(Fish_i_yes)[1]
  data_summary[1,which(choices==i)] = Income = mean(Fish_i_yes$income)/1000
  for (j in choices){
    Fish_j = Fish[Fish$alt == j,]
    Fish_i_j = Fish_j[(Fish_i$mode == 1),]
    if (which(choices==i) == which(choices==j))
    {Price = mean(Fish_i_yes$price)
     catch = mean(Fish_i_yes$catch)}
    else
    {Price = mean(Fish_i_j$price)
     catch = mean(Fish_i_j$catch)}
    data_summary[1+which(choices==j),which(choices==i)] = Price
    data_summary[5+which(choices==j),which(choices==i)] = catch
  }
}
data_summary

#-- 最简单的多项选择模型是，将二项选择模型logit做直接推广
#-- 根据多项式的不同回归元分布有所不同，将模型分为3类：
#---- 1、回归元随着选项不同而不同，例如：选择岸上和船上的catch rate不同
#---- 2、回归元不随选项不同而不同，例如：Income不随选择在哪里钓鱼而变化
#---- 3、混合logit，将以上两种模型混合起来


#------------Conditional Logit:选项变化回归元-------------
# 模型拟合
library(Formula)
CL = formula(mode ~ 0 + price + catch | 1 | 1)  # 依次是x(ik)*beta, x(i)*beta(k), x(ik)*beta(k)
fit <- mnlogit(CL, Fish,ncores = 2)

# 结果输出
summary(fit)
attributes(summary(fit))
CLcoef = summary(fit)$coefficients
CLprob = summary(fit)$probabilities

# 计算拟合优度
Lfit = summary(fit)$logLik
CL0 = formula(mode ~ 1 | 1 | 1) # 拟合仅包含截距项的模型
fit <- mnlogit(CL0, Fish,ncores = 2)
L0 = summary(fit)$logLik
Rsqure = 1- Lfit/L0 # 利用似然比计算Rsqure

# 计算边际效应（利用边际效应公式自己计算）
#---价格上涨100元的边际效应
n = dim(CLprob)[1]
ones = matrix(data = rep(1,n),n,1)
zeros = matrix(data = rep(0,n),n,1)
margin = matrix(data = NA,4,4)
for (j in 1:4){
  for (k in 1:4){
    if (j == k)
    {marg = t(CLprob[,j]) %*% (ones - CLprob[,k]) * CLcoef[1]/n *100} 
    else 
    {marg = t(CLprob[,j]) %*% (zeros - CLprob[,k]) * CLcoef[1]/n *100}
    margin[j,k] = marg
  }
}
margin

#---catch rate提高1个单位的边际效应
n = dim(CLprob)[1]
ones = matrix(data = rep(1,n),n,1)
zeros = matrix(data = rep(0,n),n,1)
margin2 = matrix(data = NA,4,4)
for (j in 1:4){
  for (k in 1:4){
    if (j == k)
    {marg = t(CLprob[,j]) %*% (ones - CLprob[,k]) * CLcoef[2]/n } 
    else 
    {marg = t(CLprob[,j]) %*% (zeros - CLprob[,k]) * CLcoef[2]/n}
    margin2[j,k] = marg
  }
}
margin2

#------------Multinomial Logit:选项不变回归元-------------
# 模型拟合
MNL = formula(mode ~  1 | income | 1)
fit <- mnlogit(MNL, Fish,ncores = 2)

# 结果输出
summary(fit)
MNLcoef = summary(fit)$coefficients
MNLprob = summary(fit)$probabilities

# 计算拟合优度
Lfit = summary(fit)$logLik
MNL0 = formula(mode ~ 1 | 1 | 1)  # 拟合仅有截距项的模型
fit <- mnlogit(MNL0, Fish,ncores = 2)
L0 = summary(fit)$logLik
Rsqure = 1- Lfit/L0

# 计算边际效应
#（方法1）自己编程序求解
n = dim(MNLprob)[1]
beta = c(0,as.vector(MNLcoef)[c(2,4,6)])
beta_i_bar = MNLprob %*% beta
margin = matrix(data = NA,1,4)
for (k in 1:4){
  beta_j = as.vector(rep(beta[k],n))
  margin[1,k] = (t(MNLprob[,k]) %*% (beta_j - beta_i_bar))/n
}
margin

#（方法2）用包求解
data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
m <- mlogit(mode  ~ 1 |income |1  , data = Fish)
summary(m)
# compute a data.frame containing the mean value of the covariates in the sample
z <- with(Fish, data.frame(price = tapply(price, index(m)$alt, mean),
                           catch = tapply(catch, index(m)$alt, mean),
                           income = mean(income)))
# compute the marginal effects (the second one is an elasticity
effects(m, covariate = "income", type = "rr", data = z)

#------------------------Mixed Logit----------------------
# 模型拟合
Mixed = formula(mode ~  price + catch | income | 1)
fit <- mnlogit(Mixed, Fish,ncores = 2)

# 结果输出
summary(fit)
MNLcoef = summary(fit)$coefficients

# 计算拟合优度
Lfit = summary(fit)$logLik
Mixed0 = formula(mode ~ 1 | 1 | 1)
fit <- mnlogit(Mixed0, Fish,ncores = 2)
L0 = summary(fit)$logLik
Rsqure = 1- Lfit/L0

# 边际效应方法如前所述
