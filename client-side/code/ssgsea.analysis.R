### R code to perform ssGSEA analysis ###

require(GSVA)
require(GSA)
require(org.Hs.eg.db)
require(dplyr)

require(stringr)
load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('server-side/RData/TCGA.breast.cancer.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')




common.genes <- intersect(rownames(MET500.log2.fpkm.matrix),rownames(CCLE.log2.rpkm.matrix))
common.genes <- intersect(common.genes,rownames(TCGA.breast.cancer.log2.fpkm.matrix))

combined.expr.matrix <- cbind(MET500.log2.fpkm.matrix[common.genes,MET500.breast.cancer.polyA.sample],
                              CCLE.log2.rpkm.matrix[common.genes,CCLE.breast.cancer.cell.line],
                              TCGA.breast.cancer.log2.fpkm.matrix[common.genes,]
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

save(file='client-side/output/ssgsea.analysis.R.output/ssgsea.analysis.RData',list=c('oncogenic.geneset.gsea.results','hallmark.geneset.gsea.results'))











