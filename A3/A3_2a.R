rm(list=ls()) 
library(survival)
library(survminer)
library(nnet)
actg <- read.table("ACTG175(speff2trial).txt",header = TRUE,sep=",")

y_surv <-Surv(actg$days, actg$cens==1) # 创建生存对象
 
# Cox分层无交叉项
model2 <- coxph(y_surv~strata(arms)+cd40+cd80+age+wtkg+hemo+homo+
                        drugs+karnof+race+gender+symptom, 
                      data = actg)
summary(model2)

#ggsurvplot(survfit(model2), data = actg, ylim = c(0.5,1))
# survfit(fit) will estimate K different baseline hazard
# functions, one for each stratum

# Cox分层有交叉项
model3 <- coxph(y_surv~strata(arms) + 
                  cd40+cd80+age+wtkg+hemo+homo + 
                  drugs+karnof+race+gender+symptom + 
                  arms:(cd40+cd80+age+wtkg+hemo+homo+
                  drugs+karnof+race+gender+symptom), 
                data = actg)
summary(model3)
# 似然比检验
anova(model3,model2)

# p>0.05 不显著，因此不同添加交叉项，不同于书上未作哑变量处理

# 同model2
dumv <- class.ind(actg$arms) # 创建哑变量
colnames(dumv) <- c("arms0","arms1","arms2","arms3")
actg <- cbind(actg ,dumv) # 合并到数据框中

model22 <- coxph(y_surv~strata(arms)+cd40+cd80+age+wtkg+hemo+homo+
                      drugs+karnof+race+gender+symptom ,data=actg)
summary(model22)

#
#model222 <- coxph(y_surv~strata(arms)+arms1:(cd40+cd80+age+wtkg+hemo+homo+
#                                  drugs+karnof+race+gender+symptom)+arms2:(cd40+cd80+age+wtkg+hemo+homo+
#                                                                             drugs+karnof+race+gender+symptom)+arms3:(cd40+cd80+age+wtkg+hemo+homo+
#                                                                                                                        drugs+karnof+race+gender+symptom)+cd40+cd80+age+wtkg+hemo+homo+
#                    drugs+karnof+race+gender+symptom ,data=actg)

anova(model3,model22)
