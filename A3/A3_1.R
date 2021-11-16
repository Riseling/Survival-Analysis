# Assignment 3
## 1

rm(list=ls()) 
library(survival)
actg <- read.table("ACTG175(speff2trial).txt",header = TRUE,sep=",")

y_surv <-Surv(actg$days, actg$cens==1) # 创建生存对象

loglog_AFT <- survreg(y_surv~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+
                        drugs+karnof+race+gender+symptom, 
                      data = actg,
                      dist = "loglogistic")

summary(loglog_AFT)
# arms cd40 cd80 karnof symptom p <0.05，结果显著与COX差drugs

# shape parameter
shape_p <- 1/loglog_AFT$scale
shape_p
# p>1 hazard先增后减

exp(loglog_AFT$coefficients)
