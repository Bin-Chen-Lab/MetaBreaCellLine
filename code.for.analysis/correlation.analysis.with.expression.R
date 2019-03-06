# R code to do gene expression correlation analysis using expression data profiled by RNASeq, from MET500.
# The 1000 most-varied genes were used

load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
source('client-side/code/util.R')
require(plyr)
require(dplyr)
require(foreach)


######## Compute 1000 most-varied genes ##########
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
CCLE.rna.seq.marker.gene.2000                 <- names(sort(rank.sd,decreasing =TRUE))[1:2000] # Well, just keep it here, never used
CCLE.rna.seq.marker.gene.3000                 <- names(sort(rank.sd,decreasing =TRUE))[1:3000] # Well, just keep it here, never used


########## MET500 profiled gene expression data with two different RNA-Seq protocals: polyA and hybrid ############
########## In this study we ONLY used samples profiled with PolyA, but the results were correlated between the two protocals.No worry!########

MET500.vs.CCLE.polyA.expression.correlation.result  <- pick.out.cell.line(MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.sample], CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.hybrid.expression.correlation.result <- pick.out.cell.line(MET500.log2.fpkm.matrix[,MET500.breast.cancer.hybrid.sample],CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)

plot(x=MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor,
     y=MET500.vs.CCLE.hybrid.expression.correlation.result$cell.line.median.cor[names(MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor)],
     xlab='PolyA',ylab='Hybrid',pch=19,cex=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,0.45),ylim=c(0,0.45)
)
lines(c(0,0.45),c(0,0.45),lwd=3,lty=2)


############## Do gene expression correaltion analysis for different BIOPSY SITE: only considering liver and lymph node  ##########
MET500.breast.cancer.polyA.sample.biopsy.site   <-  MET500.sample.meta[MET500.breast.cancer.polyA.sample,'biopsy.site']
MET500.breast.cancer.polyA.liver.sample         <-  MET500.breast.cancer.polyA.sample[MET500.breast.cancer.polyA.sample.biopsy.site == 'LIVER']
MET500.breast.cancer.polyA.lymph.sample         <-  MET500.breast.cancer.polyA.sample[MET500.breast.cancer.polyA.sample.biopsy.site == 'LYMPH_NODE']

MET500.vs.CCLE.polyA.liver.expression.correlation.result <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.liver.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.polyA.lymph.expression.correlation.result <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.lymph.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)



############ Do gene expression correaltion analysis for different subtypes ############

MET500.vs.CCLE.polyA.LumB.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumB.sample],  expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.polyA.Basal.expression.correlation.result <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample], expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.polyA.Her2.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Her2.sample],  expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.polyA.LumA.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumA.sample],  expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)


############ Save data #############
save(file = 'client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData',
     list = c('MET500.vs.CCLE.polyA.expression.correlation.result','MET500.vs.CCLE.hybrid.expression.correlation.result',
              'MET500.vs.CCLE.polyA.liver.expression.correlation.result','MET500.vs.CCLE.polyA.lymph.expression.correlation.result',
              'MET500.vs.CCLE.polyA.LumB.expression.correlation.result','MET500.vs.CCLE.polyA.Basal.expression.correlation.result','MET500.vs.CCLE.polyA.Her2.expression.correlation.result','MET500.vs.CCLE.polyA.LumA.expression.correlation.result'
            )
     )


