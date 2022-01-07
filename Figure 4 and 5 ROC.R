

#install.packages("survivalROC")

library(survivalROC)
setwd("D:\\0te\\铁死亡lncRNA -黑色素瘤\\23ROC")             #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =1, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="1rocTrain.pdf")


rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =2, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="2rocTrain.pdf")




rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =3, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="3rocTrain.pdf")




rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =4, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="4rocTrain.pdf")

     

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =5, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="5rocTrain.pdf")

             #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =8, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="8rocTrain.pdf")

          #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =0.5, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="0.5rocTrain.pdf")


         #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =6, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="6rocTrain.pdf")

           #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =7, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="7rocTrain.pdf")

           #设置工作目录

rocPlot=function(inputFile=null,outPdf=null){
		rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
		pdf(file=outPdf,width=5.5,height=5.5)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
		    predict.time =10, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
			  xlab="False positive rate", ylab="True positive rate",
			  main=paste("ROC curve (", "AUC = ",sprintf("%0.3f",roc$AUC),")"),
			  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		abline(0,1)
		dev.off()
}

#绘制train组ROC曲线
rocPlot(inputFile="geneRisk.txt",outPdf="10rocTrain.pdf")
