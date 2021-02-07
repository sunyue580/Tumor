##参考FigureYa211
####输入文件是output_combined_expr.txt——行是基因名，列是样本名
library(estimate)
source("CIBERSORT.R")

filterCommonGenes(input.f = "output_combined_expr.txt", output.f = "output_combined_expr.gct", id = "GeneSymbol")
estimateScore("output_combined_expr.gct", "output_combined_expr_estimate_score.txt", platform="affymetrix")
#这步运行时间较长
ciber.res <- CIBERSORT(sig_matrix = "LM22.txt",
                       mixture_file = "output_combined_expr.txt",
                       perm = 100,
                       QN = TRUE)
write.table(ciber.res,"CIBERSORT result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
