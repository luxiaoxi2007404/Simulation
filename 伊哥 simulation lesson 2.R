#-----------------------------------------------------------
#------------------------14章案例---------------------------
#-----------------------------------------------------------

# 数据预处理（仅取出2个y的可能性）
library(mnlogit)
data(Fish)
fix(Fish)
Fish_charter = Fish[Fish$alt == "charter",]
Fish_pier = Fish[Fish$alt == "pier",]
Fish_charter$lnrep = log(Fish_charter$price/Fish_pier$price)
Fish_pier$lnrep = log(Fish_charter$price/Fish_pier$price)
Fish_c = Fish_charter[Fish_charter$mode == 1 | Fish_pier$mode == 1,]

#--------------- logit模型拟合---------------
attach(Fish_c)
model1 = glm(mode~lnrep,binomial(link = "logit") )
summary(model1)  # 与465页结果比较

#--------------手动极大似然估计--------------

#------估计模型的beta值，运用牛顿迭代求解
# 课本上的公式14-4
lnL = function(beta)-sum(mode*log(exp(beta[1]+beta[2]*lnrep)/(1+exp(beta[1]+beta[2]*lnrep)))+(1-mode)*log(1-exp(beta[1]+beta[2]*lnrep)/(1+exp(beta[1]+beta[2]*lnrep))))
lnL(c(1,-1))                                                                                                    
?nlm
nlm(lnL,c(1,-1))  # 迭代求值方法

#------计算拟合优度
# 计算退化模型的似然值
lnL0 = function(beta)-sum(mode*log(exp(beta[1])/(1+exp(beta[1])))+(1-mode)*log(1-exp(beta[1])/(1+exp(beta[1]))))
nlm(lnL,1)
# 利用课本473页下面公式计算Rsqure，结果可与465页结果比较
Rsqure = 1-nlm(lnL,c(1,-1))$minimum/nlm(lnL0,1)$minimum
Rsqure

#------求解系数的显著性：运用似然比检验
chisq = 2*(nlm(lnL0,1)$minimum - nlm(lnL,c(1,-1))$minimum)
qchisq(0.95,1)  # 返回chisq分布分位数
1-pchisq(336.4694,1)  # 返回P值

#------计算系数的方差：运用bootstrap方法
# 渐进方差有公式14-7，但是逆矩阵不好求，故改用bootstrap
boot = 1000
n = length(mode)
beta_b = matrix(rep(NA,2*boot),boot,2)
for(b in 1:boot){
  index = sample(1:n,n,replace = TRUE)
  mode_b = mode[index]
  lnrep_b = lnrep[index]
  lnL = function(beta)-sum(mode_b*log(exp(beta[1]+beta[2]*lnrep_b)/(1+exp(beta[1]+beta[2]*lnrep_b)))+(1-mode_b)*log(1-exp(beta[1]+beta[2]*lnrep_b)/(1+exp(beta[1]+beta[2]*lnrep_b))))
  beta_b[b,] = nlm(lnL,c(1,-1))$estimate
}

sd(beta_b[,1])  # 可与summary的model1的系数的方差作比较
sd(beta_b[,2])

###########################End ################################