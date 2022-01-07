#

#install.packages("survival")

library(survival)
setwd("D:\\0te\\lncRNA -黑色素瘤\\15lasso的单基因")        #设置工作目录

#定义生存曲线函数
surPlot=function(data=null, gene=null, outPdf=null){
		a=ifelse(data[,gene]<=median(data[,gene]),"Low expression","High expression")
		diff=survdiff(Surv(futime, fustat) ~a,data = data)
		pValue=1-pchisq(diff$chisq,df=1)
		pValue=signif(pValue,4)
		pValue=format(pValue, scientific = TRUE)
		fit <- survfit(Surv(futime, fustat) ~ a, data = data)
		pdf(file=outPdf,width=5.5,height=5)
		plot(fit,
		     lwd=2,
		     col=c("red","blue"),
		     xlab="Time (year)",
		     ylab="Survival rate",
		     main=paste("Survival curve (p=", pValue ,")",sep=""),
		     mark.time=T)
		legend("topright",
		       c(paste0(gene," high expression"), paste0(gene," low expression") ),
		       lwd=2,
		       col=c("red","blue"))
		dev.off()
}

#读取输入文件
rtAll=read.table("geneRisk.txt",header=T,sep="\t",check.names=F)

#对构建模型的基因绘制生存曲线
for(i in colnames(rtAll[,4:(ncol(rtAll)-2)]) ){
		surPlot(data=rtAll, gene=i, outPdf=paste0("survival.",i,".pdf"))
}

##