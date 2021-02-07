setwd('C:\\Users\\think\\Desktop\\data3')

#################################1、数据预处理
rm(list=ls())
a <- read.table('GSE112996_merged_fpkm_table.txt.gz',
                header = T,
                row.names=1)
raw_data<- a[,-1]
###表型信息提取
pheno <- read.csv(file = 'GSE112996_series_matrix.txt')
pheno <- data.frame(num1 = strsplit(as.character(pheno[42,]),split='\t')[[1]][-1],
                    num2 = gsub('patient: No.','P',strsplit(as.character(pheno[51,]),split='\t')[[1]][-1]))
####数据过滤
data<- a[!apply(raw_data,1,sum)==0,]
####去除重复基因名的行，归一化
data$median=apply(data[,-1],1,median)
data=data[order(data$GeneName,data$median,decreasing = T),]
data=data[!duplicated(data$GeneName),]
rownames(data)=data$GeneName
uni_matrix <- data[,grep('\\d+',colnames(data))]  #\d 表示数字0-9； +表示匹配≥1次
uni_matrix <- log2(uni_matrix+1)
colnames(uni_matrix)<- gsub('X','',gsub('\\.','\\-',colnames(uni_matrix)))
uni_matrix<- uni_matrix[,order(colnames(uni_matrix))]
save(uni_matrix,pheno,file = 'uni_matrix.Rdata')

#################################2、ssGSEA分析
##加载包
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

##载入数据
load('uni_matrix.Rdata')
# gene_set<-read.table("mmc3.xlsx")[, 1:2]
gene_set<-read.table("mmc3.txt",skip = 2,header = F,sep = "\t")
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
library(pheatmap)
gsva_matrix1<- t(scale(t(gsva_matrix)))
gsva_matrix1[gsva_matrix1< -2] <- -2
gsva_matrix1[gsva_matrix1>2] <- 2
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti<- gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor
pro<- gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor
non <- !(anti|pro)
gsva_matrix1<- rbind(gsva_matrix1[anti,],gsva_matrix1[pro,],gsva_matrix1[non,])
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
annotation_col = data.frame(patient=pheno$num2)
rownames(annotation_col)<-colnames(uni_matrix)
bk = unique(c(seq(0,1, length=100)))
pheatmap(nor_gsva_matrix1,
         show_colnames = F,
         cluster_rows = F,cluster_cols = F,
         annotation_col = annotation_col,
         breaks=bk,cellwidth=5,cellheight=5,
         fontsize=5,gaps_row = c(12,20),
         #filename = 'ssgsea.pdf',
         width = 8)
save(gsva_matrix,gsva_matrix1,pheno,file = 'score.Rdata')

#################################3、计算score加和后，ggplot2进行绘图
rm(list=ls())
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
load('score.Rdata')
anti<- as.data.frame(gsva_matrix1[gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor,])
pro<- as.data.frame(gsva_matrix1[gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor,])
anti_n<- apply(anti,2,sum)
pro_n<- apply(pro,2,sum)
patient <- pheno$num2[match(colnames(gsva_matrix1),pheno$num1)]
library(ggplot2)
data <- data.frame(anti=anti_n,pro=pro_n,patient=patient)
anti_pro<- cor.test(anti_n,pro_n,method='pearson')
gg<- ggplot(data,aes(x = anti, y = pro),color=patient) + 
  xlim(-20,15)+ylim(-15,10)+
  labs(x="Anti-tumor immunity", y="Pro-tumor suppression") +
  geom_point(aes(color=patient),size=3)+
  geom_smooth(method='lm')+
  annotate("text", x = -5, y =7.5,label=paste0('R=',round(anti_pro$estimate,4),'\n','p<0.001'))
#ggsave(gg,filename = 'cor.pdf', height = 6, width = 8)