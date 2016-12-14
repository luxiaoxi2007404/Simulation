

##############################蒙特卡洛模拟############################

# 设定参数，样本量
n = 100
beta1 = 3
beta2 = 4

# 做1个模拟求拟合系数【数据生成过程是logit，用logit拟合】
x = rnorm(n)
y = runif(exp(beta1 + beta2 * x)/1+exp(beta1 + beta2 * x))
mean(y)   # 查看y的样本类别的比例
model1 = glm(y~x,family = binomial(link = logit))
summary(model1)
attributes(model1)
summary(model1)$coefficients[,1]

# 循环500次该模拟
sim = 500
model1_result = matrix(data = NA, nrow = sim, ncol = 2)
for (i in 1:sim){
  x = rnorm(n)
  y = runif(exp(beta1 + beta2 * x)/1+exp(beta1 + beta2 * x))
  model1 = glm(y~x,family = binomial(link = logit))
  model1_result[i,] = summary(model1)$coefficients[,1]
}
model1_result

# 求系数的MSE
var = var(model1_result[,1])
bias = mean(model1_result[,1])-beta1
mse = var + bias^2


# 换用probit拟合【数据生成过程还是logit】
sim = 500
model1_result = matrix(data = NA, nrow = sim, ncol = 2)
for (i in 1:sim){
  x = rnorm(n)
  y = runif(exp(beta1 + beta2 * x)/1+exp(beta1 + beta2 * x))
  model1 = glm(y~x,family = binomial(link = probit))
  model1_result[i,] = summary(model1)$coefficients[,1]
}
model1_result


var = var(model1_result[,1])
bias = mean(model1_result[,1])-beta1
mse = var + bias^2

#   结论：通过比较发现，logit拟合的效果比probit好一点点，大致一致