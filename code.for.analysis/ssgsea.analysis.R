### R code to perform ssGSEA analysis ###

require(GSVA)
require(GSA)
require(org.Hs.eg.db)
require(dplyr)
require(stringr)
require(foreach)

load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('server-side/RData/TCGA.breast.cancer.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')
load('client-side/output/organize.breast.cancer.organoid.data.R.output/organoid.RData')


## Run ssGSEA
common.genes <- intersect(rownames(MET500.log2.fpkm.matrix),rownames(CCLE.log2.rpkm.matrix))
#common.genes <- intersect(common.genes,rownames(TCGA.breast.cancer.log2.fpkm.matrix)) # good bye, TCGA
common.genes <- intersect(common.genes,rownames(organoid.log2.rpkm.matrix))


combined.expr.matrix <- cbind(MET500.log2.fpkm.matrix[common.genes,MET500.breast.cancer.polyA.sample],
                              CCLE.log2.rpkm.matrix[common.genes,CCLE.breast.cancer.cell.line],
                              #TCGA.breast.cancer.log2.fpkm.matrix[common.genes,],    # good bye, TCGA
                              organoid.log2.rpkm.matrix[common.genes,]
                             )




ensemble.to.entrez.mapping       <- revmap(org.Hs.egENSEMBL) %>% as.list
common.genes                     <- common.genes[common.genes %in% names(ensemble.to.entrez.mapping)]
gene.id.list                     <- ensemble.to.entrez.mapping[common.genes]
l                                <- sapply(gene.id.list,length)
common.genes                     <- common.genes[l == 1]
combined.expr.matrix             <- combined.expr.matrix[common.genes,]
rownames(combined.expr.matrix)   <- ensemble.to.entrez.mapping[common.genes]

msigdb          <-  GSA.read.gmt("client-side/meta.data/c6.all.v6.1.entrez.gmt")
genesets        <-  msigdb$genesets
names(genesets) <-  msigdb$geneset.names 
oncogenic.geneset.gsea.results    <-  gsva(combined.expr.matrix, genesets, method = 'ssgsea') #ggsea


msigdb          <-  GSA.read.gmt("client-side/meta.data/h.all.v6.1.entrez.gmt")
genesets        <-  msigdb$genesets
names(genesets) <-  msigdb$geneset.names 
hallmark.geneset.gsea.results    <-  gsva(combined.expr.matrix , genesets, method = 'ssgsea') #ggsea


# Compare between MET500 and cell line, MET500 and organoid
MET500.breast.cancer.polyA.non.Basal.sample <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample     <- c(MET500.breast.cancer.polyA.Basal.sample)
CCLE.non.Basal.cell.line                    <- c(CCLE.breast.cancer.Her2.cell.line,CCLE.breast.cancer.LumA.cell.line,CCLE.breast.cancer.LumB.cell.line)
non.Basal.organoid                          <- c(LumA.organoid,LumB.organoid,Her2.organoid)

# well, it is called aov.p.value, but actually I used wilcoxon rank test to get the p-value. Pay attention!
non.basal.gsea.results <- hallmark.geneset.gsea.results[,c(MET500.breast.cancer.polyA.non.Basal.sample,non.Basal.organoid,CCLE.non.Basal.cell.line)]
non.basal.aov.p.value <- foreach(i = 1:nrow(non.basal.gsea.results),.combine='rbind') %do% {
    df1 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,MET500.breast.cancer.polyA.non.Basal.sample],
                      data.source = 'MET500'
                      )
    df2 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,non.Basal.organoid],
                      data.source = 'organoid'
                       )
    df3 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,CCLE.non.Basal.cell.line],
                      data.source = 'cell.line'
                      )
    c( wilcox.test(df1$gsea.score,df2$gsea.score)$p.value, wilcox.test(df1$gsea.score,df3$gsea.score)$p.value)
#     well, to do DA analysis, let us stick on wilcoxon rank test
#     df.o             <- rbind(df1,df2)
#     df.o$data.source <- factor(df.o$data.source,level=c('organoid','MET500'))
#     df.c             <- rbind(df1,df3)
#     df.c$data.source <- factor(df.c$data.source,level=c('cell.line','MET500'))
#     
#     res.aov.o <- aov(gsea.score ~ data.source, data = df.o)
#     res.aov.c <- aov(gsea.score ~ data.source, data = df.c)
#     
#     c(summary(res.aov.o)[[1]][["Pr(>F)"]][1],summary(res.aov.c)[[1]][["Pr(>F)"]][1])
}
non.basal.aov.fdr <- cbind(-1 * log10(p.adjust(non.basal.aov.p.value[,1],method='fdr')),
                           -1 * log10(p.adjust(non.basal.aov.p.value[,2],method='fdr'))
                           )
rownames(non.basal.aov.fdr) <- rownames(non.basal.gsea.results)
colnames(non.basal.aov.fdr) <- c('organoid','cell.line')



basal.gsea.results <- hallmark.geneset.gsea.results[,c(MET500.breast.cancer.polyA.Basal.sample,Basal.organoid,CCLE.breast.cancer.Basal.cell.line)]
basal.aov.p.value <- foreach(i = 1:nrow(basal.gsea.results),.combine='rbind') %do% {
  df1 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,MET500.breast.cancer.polyA.Basal.sample],
                    data.source = 'MET500'
  )
  df2 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,Basal.organoid],
                    data.source = 'organoid'
  )
  df3 <- data.frame(gsea.score =hallmark.geneset.gsea.results[i,CCLE.breast.cancer.Basal.cell.line],
                    data.source = 'cell.line'
  )
  c( wilcox.test(df1$gsea.score,df2$gsea.score)$p.value, wilcox.test(df1$gsea.score,df3$gsea.score)$p.value)
  
#   df.o             <- rbind(df1,df2)
#   df.o$data.source <- factor(df.o$data.source,level=c('organoid','MET500'))
#   df.c             <- rbind(df1,df3)
#   df.c$data.source <- factor(df.c$data.source,level=c('cell.line','MET500'))
#   
#   res.aov.o <- aov(gsea.score ~ data.source, data = df.o)
#   res.aov.c <- aov(gsea.score ~ data.source, data = df.c)
#   
#   c(summary(res.aov.o)[[1]][["Pr(>F)"]][1],summary(res.aov.c)[[1]][["Pr(>F)"]][1])
}

basal.aov.fdr <- cbind(-1 * log10(p.adjust(basal.aov.p.value[,1],method='fdr')),
                       -1 * log10(p.adjust(basal.aov.p.value[,2],method='fdr'))
)
rownames(basal.aov.fdr) <- rownames(basal.gsea.results)
colnames(basal.aov.fdr) <- c('organoid','cell.line')


save(file='client-side/output/ssgsea.analysis.R.output/ssgsea.analysis.RData',list=c('oncogenic.geneset.gsea.results','hallmark.geneset.gsea.results','non.basal.gsea.results','basal.gsea.results','basal.aov.fdr','non.basal.aov.fdr'))

