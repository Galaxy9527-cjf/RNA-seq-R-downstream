# RNA-seq-R-downstream

#指定工作目录（原始文件所在目录）
setwd("d:/rna-seq/counts")

#清空环境变量
rm(list=ls())

#设置options(作用暂时不明)
options(stringsAsFactors = F)

#读入feature count的TXT结果文件
a=read.table("count_matrix.txt",header = T)
#PS,如果读取是显示“Error in scan
#(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
#line 8687 did not have 13 elements"之类，需要用excel打开txt调整相应行内容。

#拆分矩阵：meta信息和表达信息。ncol()表示最大列数。
meta=a[,1:6]
exprset=a[,7:ncol(a)]

#热图绘制 
#解释：scale：归一化；cor：相关性；log2：取2的对数；
#+1：防止0值产生（因有归一化，不影响结果）
library(pheatmap)
png('heatmap.png')
pheatmap(scale(cor(log2(exprset+1))))
dev.off()

#相关性分析
m=cor(log2(exprset+1))

#airway包安装
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")


#设置group_list（即SRR所属的组别，对照还是实验组）
#暂时使用，直接导入会出问题。
group_list=c('Control','Control','Control','VMC','VMC','VMC','VMC')


#DESeq2的使用
library(DESeq2)
(colData <- data.frame(row.names = colnames(exprset),group_list = group_list))
dds <- DESeqDataSetFromMatrix(countData = exprset, 
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)

#差异分析结果输出
#contrast的内容由自己的数据做相应更改，但是顺序不要变
res <- results(dds,
               contrast = c("group_list","VMC","Control"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)

#差异结果DEG中NA数据的筛除
DEG = na.omit(DEG)
