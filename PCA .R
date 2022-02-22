Exp<-read.csv('expression.csv',row.names = 1)
exp<-Exp[1:9]
#write.table(exp, file = "exp.csv")
colnames(exp)<-c("C1","C2","C3","L1","L2","L3","H1","H2","H3")
exp2<-t(exp)

#我们使用 FactoMineR 包中的方法，实现 PCA 分析和聚类添加
#install.packages("FactoMineR")
library(FactoMineR)

#样本中基因表达值的 PCA 分析
gene.pca <- PCA(exp2, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)  #PCA 简图

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)

#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2],2)
pca_eig2 <- round(gene.pca$eig[2,2],2)


#设置分组
A1<-c('C1','C2','C3','L1','L2','L3','H1','H2','H3')
A2<-c('C','C','C','L','L','L','H','H','H')
Group<-rbind(A1,A2)
group<-t(Group)
colnames(group)=c("Sample","group")
pca_sample <- cbind(pca_sample, group)


#ggplot2 绘制二维散点图
library(ggplot2)

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 3) +  #根据样本坐标绘制二维散点图
  scale_color_manual(values = c('orange', 'green','red')) +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p

#添加 95% 置信椭圆，可用于表示对象分类，但只能适用于各组样本数大于 5 个的情况
#p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

#p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  #scale_fill_manual(values = c('orange', 'purple'))

#多边形连接同类别对象边界的样式，适用于各组样本数大于 3 个的情况
library(plyr)
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])

p + geom_polygon(data = cluster_border, aes(color = group), fill = NA, show.legend = FALSE)

p + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.1, show.legend = F)+
  scale_fill_manual(values = c('orange', 'green','red'))

#相关性分析
library(corrplot)

cordata<-cor(exp)

corrl$p

corrplot(cordata,tl.col = 'black',order='hclust',p.mat = corrl$p,insig = 'blank')
corrplot(cordata,tl.col = 'black',order='hclust',p.mat = corrl$p,insig = 'label_sig',
         sig.level = c(0.0001,0.001,0.005),pch.cex = 1,pch.col = 'red')
         
corrplot(cordata,tl.col = 'black',order='hclust',p.mat = corrl$p,insig = 'label_sig',
         sig.level = c(0.0001,0.001,0.005),pch.cex = 1,pch.col = 'red',type = 'lower')
corrplot(cordata,tl.col = 'black',order='hclust',p.mat = corrl$p,insig = 'label_sig',
         sig.level = c(0.0001,0.001,0.005),pch.cex = 1,pch.col = 'red',type = 'lower',method = 'color')
