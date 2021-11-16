library(etm)
library(survival)

data(abortion)
head(abortion)

index = abortion$cause==2
abortion$status[index] = 0
abortion$status[-index] = 1
y_surv <- Surv(abortion$entry,abortion$exit,abortion$status==1) 
coxph(y_surv~group,data = abortion)

