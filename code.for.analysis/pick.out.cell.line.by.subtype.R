### R code to pick out good cell lines for each subtype (microarray data and MET500 RNASeq data)
### Based on the hypopthesis that given a random cell line, its median correlation with metastatic breast cancer samples of a specific subtype is normally distributed, 
### we first fit a normal distribution and then used it as NULL distribution to compute p-value

load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
load('client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output//CCLE.breast.cancer.cell.line.meta.RData')
require(dplyr)

CCLE.non.breast.cancer.cell.line <- setdiff(names(MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor),CCLE.breast.cancer.cell.line)

m                                <- MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad                                 # fit the normal distribtuion 
p.value                          <- pnorm(q=MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE) # compute p-value
fdr.vec                          <- p.adjust(p.value,method='fdr')
MET500.Basal.good.cell.line      <- names(fdr.vec)[fdr.vec <= 0.01]
MET500.Basal.good.cell.line      <- names(fdr.vec)[fdr.vec <= 0.05] # should use 0.05


m                                <- MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
LumA.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
MET500.LumA.good.cell.line       <- LumA.good.cell.line[fdr.vec[LumA.good.cell.line] %>% order]

m                                <- MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
LumB.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
MET500.LumB.good.cell.line       <- LumB.good.cell.line[fdr.vec[LumB.good.cell.line] %>% order]

m                                <- MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
Her2.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
MET500.Her2.good.cell.line       <- Her2.good.cell.line[fdr.vec[Her2.good.cell.line] %>% order]



CCLE.non.breast.cancer.cell.line <- setdiff(names(microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor),CCLE.breast.cancer.cell.line)

m                                <- microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
microarray.Basal.good.cell.line  <- names(fdr.vec)[fdr.vec <= 0.01]

m                                <- microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
LumA.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
microarray.LumA.good.cell.line   <- LumA.good.cell.line[fdr.vec[LumA.good.cell.line] %>% order]

m                                <- microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
LumB.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
microarray.LumB.good.cell.line   <- LumB.good.cell.line[fdr.vec[LumB.good.cell.line] %>% order]

m                                <- microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% median
s                                <- microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.non.breast.cancer.cell.line] %>% mad
p.value                          <- pnorm(q=microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],mean = m,sd = s,lower.tail = FALSE)
fdr.vec                          <- p.adjust(p.value,method='fdr')
Her2.good.cell.line              <- names(fdr.vec)[fdr.vec <= 0.01]
microarray.Her2.good.cell.line   <- Her2.good.cell.line[fdr.vec[Her2.good.cell.line] %>% order]



save(file = 'client-side/output/pick.out.cell.line.by.subtype.R.output/pick.out.cell.line.by.subtype.RData',
     list=c('MET500.Basal.good.cell.line','MET500.LumA.good.cell.line','MET500.Her2.good.cell.line','MET500.LumB.good.cell.line',
            'microarray.Basal.good.cell.line','microarray.LumA.good.cell.line','microarray.Her2.good.cell.line','microarray.LumB.good.cell.line'
            )
     )
