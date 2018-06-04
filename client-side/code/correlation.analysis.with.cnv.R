# R code to pick out matched CCLE cell lines according to median CNV correlation (for both MET500 and TCGA)
# Well, the final results may not be included in the manuscript, but let us keep this file.

load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
source('client-side/code/util.R') # load the pick.cell.line function

MET500.vs.CCLE.cnv.correlation.result <- pick.out.cell.line(expr.of.cell.lines = CCLE.gene.cnv,expr.of.samples  = MET500.breast.cancer.gene.cnv,marker.gene = rownames(MET500.breast.cancer.gene.cnv))
TCGA.vs.CCLE.cnv.correlation.result   <- pick.out.cell.line(expr.of.cell.lines = CCLE.gene.cnv,expr.of.samples  = TCGA.breast.cancer.gene.cnv,marker.gene   = rownames(MET500.breast.cancer.gene.cnv))

save(file='client-side/output/correlation.analysis.with.cnv.R.output/correlation.analysis.with.cnv.result.RData',
    list=c('MET500.vs.CCLE.cnv.correlation.result','TCGA.vs.CCLE.cnv.correlation.result'))




