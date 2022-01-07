#install.packages("glmnet")


library(glmnet)
library(survival)
setwd("C:\\Users\\scikuangren\\Desktop\\TARGET\\2_lasso")#ÐÞ¸Ä

mydata<-read.table("TARGET_dev.txt",header=T,sep="\t") 
for (a in names(mydata)[c(4:11)]){
  mydata[,a] <- as.factor(mydata[,a])
}

v1<-as.matrix(mydata[,c(4:11)])
v2 <- data.matrix(Surv(mydata$Survival_time,mydata$Vital_Status))

myfit <- glmnet(v1, v2, family = "cox")
pdf("lambda.pdf")
plot(myfit, xvar = "lambda", label = TRUE)
dev.off()

myfit2 <- cv.glmnet(v1, v2, family="cox")
pdf("min.pdf")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
dev.off()

myfit2$lambda.min
coe <- coef(myfit, s = myfit2$lambda.min)
act_index <- which(coe != 0)
act_coe <- coe[act_index]
row.names(coe)[act_index]

#result
"WBC"              "Bone_Marrow_Site" "CNS_Site"  

