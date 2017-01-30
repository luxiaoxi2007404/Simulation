library(foreign)
data = read.dta("C:/Users/Administrator/Documents/R/微观计量/ema1996.dta", convert.factors=F)


library("survival")
head(data)
data = subset(data,censor1<=1 & spell >0)
km_strike = survfit(Surv(spell, censor1) ~ 1, data = data, conf.type = "log-log")
km_strike

# 图1
plot(km_strike,main = "总生存函数估计",xlab = "时间",ylab = "比率",xlim = c(0,30))
legend("topright",legend=c("95%置信区间","生存估计","95%置信区间"),cex = 0.7,lty = c(2,1,2))
# 表2
survivor_function = summary(km_strike)
cumulative_hazard = cumsum(survivor_function$n.event/survivor_function$n.risk)
# 图3
km_strike_UI = survfit(Surv(spell, censor1) ~ ui, data = data, conf.type = "log-log")
km_strike_UI
plot(km_strike_UI,main = "是否UI状态下的函数估计",xlab = "时间",ylab = "比率",xlim = c(0,30),lty = c(2,1))
legend("topright",legend=c("无UI(UI = 0)","接受UI(UI = 1)"),cex = 0.7,lty = c(2,1))
# 图4
plot(survivor_function$time,cumulative_hazard,main = "总累计风险估计",xlab = "时间",ylab = "风险",xlim = c(0,30),type = "l") 
upper_conf =  cumulative_hazard + survivor_function$std.err*2
lower_conf =  cumulative_hazard - survivor_function$std.err*2 
lines(survivor_function$time,upper_conf,xlim = c(0,30),lty = 2)
lines(survivor_function$time,lower_conf,xlim = c(0,30),lty = 2)
legend("bottomright",legend=c("风险估计","95%置信区间"),cex = 0.7,lty = c(1,2))
# 图5
UI_survivor_function = summary(km_strike_UI)
UI0_cumulative_hazard = cumsum((UI_survivor_function$n.event/UI_survivor_function$n.risk)[1:21])
UI1_cumulative_hazard = cumsum((UI_survivor_function$n.event/UI_survivor_function$n.risk)[22:45])
plot(UI_survivor_function$time[1:21],UI0_cumulative_hazard,main = "UI状态下的累计风险估计",xlab = "时间",ylab = "风险",xlim = c(0,30),type = "l") 
lines(UI_survivor_function$time[22:45],UI1_cumulative_hazard,xlim = c(0,30),lty = 2)
legend("bottomright",legend=c("无UI(UI = 0)","接受UI(UI = 1)"),cex = 0.7,lty = c(1,2))
# 4种模型
## 指数
exponential_fit = survreg(Surv(spell, censor1)~ 0 + reprate + disrate + ui + reprate*ui +disrate*ui +logwage + constr + tenure + slack
                                                + abolpos + explose + stateur + houshead + married + female + child + ychild +　nonwhite　
                          　　　　　　　　　　　+age+schlt12+schgt12+smsa+bluecoll+mining+transp+trade+fire+services+ pubadmin
                          　　　　　　　　　　　+year85+ year87+year89+midatl　+　encen+ wncen+southatl+escen+wscen+mountain+pacific
                                                ,data = data, dist="exponential") 
summary(exponential_fit)

## 威布尔
exponential_fit = survreg(Surv(spell, censor1)~ reprate + disrate + ui + reprate*ui +disrate*ui +logwage + constr + tenure + slack
                          + abolpos + explose + stateur + houshead + married + female + child + ychild +　nonwhite　
                          +age+schlt12+schgt12+smsa+bluecoll+mining+transp+trade+fire+services+ pubadmin
                          +year85+ year87+year89+midatl+encen+ wncen+southatl+escen+wscen+mountain+pacific
                          ,data = data, dist="weibull") 
summary(exponential_fit)

#fitw <- flexsurvreg(formula = Surv(spell, censor1)~  reprate + disrate + ui + reprate*ui +disrate*ui +logwage + constr + tenure + slack
#                   + abolpos + explose + stateur + houshead + married + female + child + ychild +　nonwhite　
#                   +age+schlt12+schgt12+smsa+bluecoll+mining+transp+trade+fire+services+ pubadmin
#                   +year85+ year87+year89+midatl　+　encen+ wncen+southatl+escen+wscen+mountain+pacific
#                   ,data = data, dist="weibull") 

## 冈珀茨
library("flexsurv")
fitw <- flexsurvreg(formula = Surv(spell, censor1)~ 0+ reprate + disrate + ui + reprate*ui +disrate*ui +logwage + constr + tenure + slack
                                        + abolpos + explose + stateur + houshead + married + female + child + ychild +　nonwhite　
                                        +age+schlt12+schgt12+smsa+bluecoll+mining+transp+trade+fire+services+ pubadmin
                                        +year85+ year87+year89+midatl+encen+ wncen+southatl+escen+wscen+mountain+pacific
                                        ,data = data, dist="gompertz") 

## coxPH
coxph_fit = coxph(Surv(spell, censor1) ~ 0 + reprate + disrate + ui + reprate*ui +disrate*ui +logwage + constr + tenure + slack
      + abolpos + explose + stateur + houshead + married + female + child + ychild +　nonwhite　
      +age+schlt12+schgt12+smsa+bluecoll+mining+transp+trade+fire+services+ pubadmin
      +year85+ year87+year89+midatl　+　encen+ wncen+southatl+escen+wscen+mountain+pacific,data) 
summary(coxph_fit)



 




