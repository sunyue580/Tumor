####参考：FigureYa71ssGSEA
library(GSVA)
library(rtracklayer)
library(RColorBrewer)
library(GenomicRanges)
library(dplyr)
library(clusterProfiler)
library(stringr)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#读入免疫细胞marker genes——来自文献补充材料S1
#https://www.cell.com/cms/10.1016/j.immuni.2013.10.003/attachment/8dc04d32-6aff-4eda-99f5-6401bcae9751/mmc1.pdf

##########################################1、将表达矩阵的行名换成entrez ID
immunity <- read.csv("immunitygene.csv", header = T)

# 去除不是免疫细胞的CellType
immunity <- immunity[!immunity$CellType %in% c("Blood vessels", "Normal mucosa", "SW480 cancer cells", "Lymph vessels"), ]

# 转成list，去除冗余基因
# 文献marker genes是在affy芯片平台上获取的，我这里选用的是Entrez_ID。也可以用Gene symbol
immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$ENTREZ_GENE_ID))
immunity <- lapply(immunity, unique)

# 保存到文件，便于过后重复使用时直接读取
save(immunity,file = "easy_input_immunity.rdata")


#先获得protein_coding gene的gene_id跟gene symbol的对应关系
#注释文件来自ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
gtf_v22 <- rtracklayer::import("gencode.v22.annotation.gtf")
gtf_v22 <- as.data.frame(gtf_v22)
gtf_v22 <- dplyr::select(gtf_v22, c("gene_id", "gene_type", "gene_name")) %>% 
  filter(., gene_type ==  "protein_coding") %>% 
  unique()
gtf_v22 <- gtf_v22[!duplicated(gtf_v22), ]

# 再加入ENTREZ ID，如果immunity.rdata里用的是gene symbol，就不需要运行下面这行
gtf_v22 <- bitr(gtf_v22$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>%
  merge(., gtf_v22, by.x = "SYMBOL", by.y = "gene_name") %>% 
  mutate(., gene_id = str_sub(gene_id, 1, 15))

# 这步获得gtf_v22.csv文件，包含v22版本protein_coding gene的gene_id、gene symbol和ENTREZ ID的对应关系，可重复使用。
write.csv(gtf_v22, "gtf_v22.csv", quote = F, row.names = F)


tcga_expr <- read.table("TCGA_readcount.genes.tpm.txt", header = T, stringsAsFactors = F, sep = "\t",row.names = NULL)
gtf_v22 <- read.csv("gtf_v22.csv")
library(data.table)
tcga_expr <- dplyr::inner_join(gtf_v22[, c("ENTREZID", "gene_id")], tcga_expr, by=c("gene_id"="row.names")) %>% .[, -2] %>% 
  aggregate(.~ENTREZID, ., median) #这步运行时间长

tcga_expr <- tibble::column_to_rownames(tcga_expr,var = "ENTREZID")    
write.csv(tcga_expr, "easy_input_expr.csv", quote = F)

##########################################2、ssGSEA计算免疫浸润
tcga_expr <- read.csv("easy_input_expr.csv", row.names = 1)
tcga_expr[1:3,1:2]
(load("easy_input_immunity.rdata"))

# 下面这一句就完成了ssGSEA
tcga_gsva <- as.data.frame(t(gsva(as.matrix(tcga_expr), immunity, method = "ssgsea")))
# 把ssGSEA结果保存到文件
write.csv(tcga_gsva, "ssGSEA_output.csv", quote = F, row.names = T)


