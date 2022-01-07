#没有安装过的包需要安装
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

source("http://bioconductor.org/biocLite.R") 
biocLite("pheatmap")

install.packages("ggplot2")

#每次分析都需要加载
library(edgeR)
library(pheatmap)
library(ggplot2)

#设定自己的工作目录
dir="C:\\Users\\scikuangren\\Desktop\\TARGET\\2_diff_miRNA"
setwd(dir)
getwd()
#读取表达矩阵数据
tcga<-read.table("miRNA.txt",header = T,sep = "\t",check.names = F)
#将数据框tcga转化矩阵
tcga=as.matrix(tcga)
#将第一列基因做行名
rownames(tcga)=tcga[,1]
#第二列到最后一列作为基因表达数据
GeneExp=tcga[,2:ncol(tcga)]
#将表达量转成数值型
TCGA=matrix(as.numeric(as.matrix(GeneExp)),nrow=nrow(GeneExp),dimnames=list(rownames(GeneExp),colnames(GeneExp)))
#对重复的基因取平均值
TCGA=avereps(TCGA)
#过滤表达量低的基因
TCGA=TCGA[rowMeans(TCGA)>1,]
#设置每个样本的实验条件,多少个为正常，多少个为癌症
design=c(rep("normal",6),rep("cancer",132))                       
#构建设计矩阵
mydesign <- model.matrix(~design)
#构建对比矩阵
mydgelist <- DGEList(counts=TCGA,group=design)
#计算标准化因子
mydgelist<- calcNormFactors(mydgelist)
#估计散度
mydgelist<- estimateCommonDisp(mydgelist)
mydgelist <- estimateTagwiseDisp(mydgelist,trend = "movingave")
#两组数据间进行双尾检验
mytest <- exactTest(mydgelist,pair = c("normal","cancer"))
#对P值进行矫正并且提取出所有的基因
Allgene<-topTags(mytest,n=10000000)
Allgene=Allgene$table
iddata<-mydgelist$pseudo.counts

#保存所有的基因
write.table(Allgene,"Allgene.txt",sep="\t",quote = F)
#保存差异基因
Diffgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC>2 | Allgene$logFC<(-2))),]
write.table(Diffgene, "Diffgene.txt",sep="\t",quote=F)
#保存上调的差异基因
Upgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC>2)),]
write.table(Upgene, "Upgene.txt",sep="\t",quote=F)
#保存下调的差异基因
Downgene = Allgene[(Allgene$FDR < 0.05 & (Allgene$logFC<(-2))),]
write.table(Downgene, "Downgene.txt",sep="\t",quote=F)
#保存所有矫正后基因的表达量
Normalizegeneexp=rbind(id=colnames(iddata),iddata)
write.table(Normalizegeneexp,"Normalizegeneexp.txt",sep="\t",quote=F,col.names=F)   
#保存差异基因的表达量
Diffgeneexp=rbind(id=colnames(iddata),iddata[rownames(Diffgene),])
write.table(Diffgeneexp,"Diffgeneexp.txt",sep="\t",quote=F,col.names=F)

#绘制热图
inputheatmap<-iddata[rownames(Diffgene),]
inputheatmap=log2(inputheatmap+1)
inputheatmap=inputheatmap[1:50,]
pdf("heatmap.pdf",30,15)
pheatmap(inputheatmap,display_numbers = F,fontsize_row=8,fontsize_col=12,color = colorRampPalette(c("green", "black", "red"))(50),cluster_cols = T,cluster_rows = T)
dev.off()

#绘制火山图
xmax<-max(Allgene$logFC)
ymax<-max(-log10(Allgene$FDR))

Allgene$sig = as.factor(ifelse(Allgene$FDR < 0.05 & abs(Allgene$logFC) > 2, 
                                 ifelse(Allgene$logFC > 2,"Up", "Down"), "Not"))

a=ggplot(Allgene,aes(logFC,-log10(FDR))) 
b=a+ geom_point(aes(color =sig))+xlim(-xmax,xmax) + ylim(0,ymax)+ 
  labs(title="Volcano",x="log2FC", y="-log10(FDR)")+theme(plot.title=element_text(hjust=0.5))+scale_color_manual(values =c("green","black", "red"))

pdf("Volcano.pdf")
b+coord_flip()+geom_vline(xintercept = 0,lty=2,lwd=1)
dev.off()
