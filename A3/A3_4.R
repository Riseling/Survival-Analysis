library(survival)
cycle <- c(1:13)
smoke <- c(29,16,17,4,3,9,4,5,1,1,1,3,7)
non_smoke <- c(198,107,55,38,18,22,7,9,5,3,6,6,12)
p_s <- data.frame(smoke,non_smoke)
rownames(p_s)[13]<- ">12"

# (a)

chisq.test(p_s)

#存在差异

#(b)

id <- 1:sum(c(smoke,non_smoke))
status <- rep(c(1,0),c(sum(smoke),sum(non_smoke)))
y1 <- rep(c(1:13),smoke)
y2 <- rep(c(1:13),non_smoke)
cycles <- factor(c(y1,y2),ordered = TRUE)
r_ps <- data.frame(id,status,cycles)

fit1 <- glm(cycles~factor(status),family=binomial(link=cloglog),data=r_ps)
s1 <- summary(fit1)
# HR
exp(fit1$coefficients[2])
exp(s1$coefficients[2])
# 置信区间
exp(confint(fit1))[2,]

fit2 <- glm(cycles~factor(status),family=binomial(link=logit),data=r_ps)
# OR
exp(fit2$coefficients[2])
# 置信区间
exp(confint(fit2))[2,]