rm(list=ls()) 
library(survival)
library(randomForestSRC)
library(ggRandomForests)
library(ggplot2)
library(survminer)


actg <- read.table("ACTG175(speff2trial).txt",header = TRUE,sep=",")

# RSF
rsf <- rfsrc(Surv(days, cens)~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+ 
        drugs+karnof+race+gender+symptom, data = actg, 
        tree.err = TRUE, splitrule = "logrank", 
        importance = "TRUE", samptype = "swr")


# 预测误差
err_rsf <-gg_error(rsf)
ggplot(data=err_rsf[!is.na(err_rsf$error),],aes(x=ntree, y=error))+ 
  geom_line()+ 
  geom_point()+ 
  xlab("Number of Trees") + ylab("OOB Error Rate")+
  theme(axis.title = element_text(size=14))

# VIMP图
# VIMP depends upon block.size, an integer value between 1 and ntree, 
# specifying the number of trees in a block used to determine VIMP. 
# When block.size=1, VIMP is calculated by tree. 
# The difference between prediction error under the perturbed predictor and 
# the original predictor is calculated for each tree and averaged over the forest. 
# This yields Breiman-Cutler VIMP (Breiman 2001).

# ggplot风格
plot(gg_vimp(rsf))+ theme(legend.position = c(0.8, 0.2)) + 
  labs(fill = "VIMP > 0") + 
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=10,colour = "black")) + 
  ylab("Variable Importance")

# 同参考文献
smp <- subsample(rsf)
plot.subsample(smp, cex.axis = .7, xlab = "Variable Importance × 100")


# partial dependence plot
## cd40 time=578中位数
cd40_pd <- plot.variable(rsf, "cd40", surv.type = "surv", 
              partial = TRUE)

# 调包
plot.variable(cd40_pd)

cd40_gpd <- gg_partial(cd40_pd)

plot(cd40_gpd)+geom_line()+ 
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=10,colour = "black")) + 
  ylab("predict survival(time=578)") + 
  xlim(c(0,1250))

## cd80
cd80_pd <- plot.variable(rsf, "cd80", surv.type = "surv", 
                         partial = TRUE)
cd80_gpd <- gg_partial(cd80_pd)
plot(cd80_gpd)+geom_line()+ 
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=10,colour = "black")) + 
  ylab("predict survival(time=578)")

## cd40 and cd80
#apd <- partial.rfsrc(rsf, partial.type="surv", 
#                     surv.type = "surv",
#                     m.target = "surv",
#                     partial.xvar ="cd40",
#                     partial.values = rsf$xvar[,"cd40"],
#                     partial.time = rsf$time.interest)

# 2C i 画图
patern0 <- data.frame(arms=c(0,1,2,3),cd40=400, cd80=500, age=25, wtkg=70,
                      hemo=0,homo=0,drugs=1,karnof=90,race=1,gender=1,symptom=0)
y_surv <-Surv(actg$days, actg$cens==1)

#model3 <- coxph(y_surv~strata(arms) + 
#                  cd40+cd80+age+wtkg+hemo+homo + 
#                 drugs+karnof+race+gender+symptom + 
#                  arms:(cd40+cd80+age+wtkg+hemo+homo+
#                          drugs+karnof+race+gender+symptom), 
#                data = actg)

#model4 <- coxph(y_surv~strata(arms)+
#                  arms*(cd40+cd80+age+wtkg+hemo+homo+
#                          drugs+karnof+race+gender+symptom)-
#                 arms, data = actg)

#surv_fit(model4, data=patern0)

model4_0 <- coxph(y_surv~cd40+cd80+age+wtkg+hemo+homo + 
                    drugs+karnof+race+gender+symptom,
                  data = actg, subset = (arms==0))

model4_1 <- coxph(y_surv~cd40+cd80+age+wtkg+hemo+homo + 
                    drugs+karnof+race+gender+symptom,
                  data = actg, subset = (arms==1))

model4_2 <- coxph(y_surv~cd40+cd80+age+wtkg+hemo+homo + 
                    drugs+karnof+race+gender+symptom,
                  data = actg, subset = (arms==2))

model4_3 <- coxph(y_surv~cd40+cd80+age+wtkg+hemo+homo + 
                    drugs+karnof+race+gender+symptom,
                  data = actg, subset = (arms==3))

new_p = data.frame(cd40=400, cd80=500, age=25, wtkg=70,
                   hemo=0,homo=0,drugs=1,karnof=90,race=1,gender=1,symptom=0)

new_model4_0 <- survfit(model4_0 ,newdata=new_p)
new_model4_1 <- survfit(model4_1 ,newdata=new_p)
new_model4_2 <- survfit(model4_2 ,newdata=new_p)
new_model4_3 <- survfit(model4_3 ,newdata=new_p)
fits <- list(arms0 = new_model4_0,
              arms1 = new_model4_1,
              arms2 = new_model4_2,
              arms3 = new_model4_3)

le.title <- list("arms=0","arms=1","arms=2","arms=3")

ggsurvplot_combine(fits,data=actg,legend.title=le.title,
                   xlim = c(0,1250),
                   ylim = c(0.7,1),
                   ylab="Survival Probability",
                   ggtheme = theme_grey(),
                   legend= "bottom",
                   legend.labs = c("0","1","2","3"))

n_test <- predict.rfsrc(rsf, newdata = patern0,importance = TRUE)
#plot.survival.rfsrc(n_test)
#plot(gg_rfsrc(rsf, by ="arms"),theme_bw())
plot(gg_rfsrc(n_test,by ="arms")) +
  ylab("Survival Probability")+  ylim(c(0.7,1))+
  xlab("Time")+theme_grey()+ xlim(c(0,1250))+
  theme(legend.position = "bottom")

# 2e
## 未分割成训练集和测试集
model1 <- coxph(y_surv~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+drugs+ 
                  karnof+race+gender+symptom ,data=actg)
survConcordance(y_surv~predict(model1))
survConcordance(y_surv~predict(model3))
model1$concordance
model3$concordance

## 70的训练，30的测试
set.seed(1234)
index <- sample(1:2319,size=0.7*2319)
train_actg <- actg[index,]
test_actg <-actg[-index,]
model1 <- coxph(Surv(days,cens)~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+drugs+ 
                  karnof+race+gender+symptom ,data=train_actg)

model3 <- coxph(Surv(days,cens)~strata(arms) + 
                  cd40+cd80+age+wtkg+hemo+homo + 
                  drugs+karnof+race+gender+symptom + 
                  arms:(cd40+cd80+age+wtkg+hemo+homo+
                          drugs+karnof+race+gender+symptom), 
                data = train_actg)

train_rfs <- rfsrc(Surv(days, cens)~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+ 
                     drugs+karnof+race+gender+symptom, data = train_actg, 
                   tree.err = TRUE, splitrule = "logrank", 
                   importance = "TRUE", samptype = "swr") 

c1 <- survConcordance(Surv(test_actg$days,test_actg$cens)~predict(model1, test_actg))
c2 <- survConcordance(Surv(test_actg$days,test_actg$cens)~predict(model3, test_actg))
c3 <- predict.rfsrc(train_rfs, newdata = test_actg, importance = TRUE)

c1$concordance
c2$concordance
1 - mean(c3$err.rate,na.rm = TRUE)

# 调包
library(dynpred)
cindex(y_surv~factor(arms)+cd40+cd80+age+wtkg+hemo+homo+drugs+ 
         karnof+race+gender+symptom,data=actg)
cindex(y_surv~strata(arms) + 
         cd40+cd80+age+wtkg+hemo+homo + 
         drugs+karnof+race+gender+symptom + 
         arms:(cd40+cd80+age+wtkg+hemo+homo+
                 drugs+karnof+race+gender+symptom),data=actg)

cd_gpd <-merge(cd40_gpd,cd80_gpd,by = "yhat",all=TRUE)

cd_gpd
