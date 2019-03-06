# R code to do gene expression correlation analysis using expression data profiled by microarray, from GEO.
# The 1000 most-varied probe sets were used



load('server-side/RData/microarray.RData')
source('client-side/code/util.R')
require(stringr)
require(genefu)
require(dplyr)
require(foreach)

## Compute the 1000 most-varied probesets
tmp                         <- CCLE.microarray.expr.matrix
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.microarray.marker.gene.1000 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


############ Do gene expression correlation analysis, ALL metastatic breast cancer samples were used! ###########
microarray.vs.CCLE.correlation.result <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix,expr.of.cell.lines = CCLE.microarray.expr.matrix,marker.gene = CCLE.microarray.marker.gene.1000)



#### Do gene expression correlation analysis for differnt BIOPSY SITE : live,lung,bone,brain,lymph.node#######
biopsy.site       <-  metastatic.breast.cancer.microarray.meta.df$biopsy.site
liver.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LIVER']  %>% as.character
lung.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LUNG']   %>% as.character
bone.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BONE']   %>% as.character
brain.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BRAIN']  %>% as.character
lymph.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LYMPH_NODE']  %>% as.character

microarray.vs.CCLE.liver.correlation.result  <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,liver.samples], expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.lung.correlation.result   <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,lung.samples],  expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.bone.correlation.result   <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,bone.samples],  expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.brain.correlation.result  <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,brain.samples], expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.lymph.correlation.result  <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,lymph.samples], expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)




########## Do gene expression correlation analysis for differnt subtypes ################

#apply genefu package to perform breast cancer subtyping
affy.probe.2.entrez.gene.id           <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/affy.probe.2.entrez.gene.id.txt", stringsAsFactors=FALSE)
colnames(affy.probe.2.entrez.gene.id) <- c('probe','EntrezGene.ID')
annot.matrix                          <- as.matrix(affy.probe.2.entrez.gene.id)
annot.matrix[,2]                      <- str_trim(annot.matrix[,2],side = 'both') 
rownames(annot.matrix)                <- annot.matrix[,'probe']
flag                                  <- annot.matrix[,2] %in% (pam50.robust$centroids.map[,3] %>% as.character)
annot.matrix                          <- annot.matrix[flag,]
metastatic.breast.cancer.microarray.expr.matrix <- apply(metastatic.breast.cancer.microarray.expr.matrix,2,function(x) x- median(x))
pam50.subtype.rs                                <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = t(metastatic.breast.cancer.microarray.expr.matrix),do.mapping = TRUE,annot = annot.matrix )

metastatic.breast.cancer.LumB.sample    <- colnames(metastatic.breast.cancer.microarray.expr.matrix)[pam50.subtype.rs$subtype == 'LumB']
metastatic.breast.cancer.Basal.sample   <- colnames(metastatic.breast.cancer.microarray.expr.matrix)[pam50.subtype.rs$subtype == 'Basal']
metastatic.breast.cancer.Her2.sample    <- colnames(metastatic.breast.cancer.microarray.expr.matrix)[pam50.subtype.rs$subtype == 'Her2']
metastatic.breast.cancer.LumA.sample    <- colnames(metastatic.breast.cancer.microarray.expr.matrix)[pam50.subtype.rs$subtype == 'LumA']
metastatic.breast.cancer.Normal.sample  <- colnames(metastatic.breast.cancer.microarray.expr.matrix)[pam50.subtype.rs$subtype == 'Normal']


microarray.vs.CCLE.LumB.correlation.result  <- pick.out.cell.line(expr.of.samples   = metastatic.breast.cancer.microarray.expr.matrix[,metastatic.breast.cancer.LumB.sample],expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.Basal.correlation.result <- pick.out.cell.line(expr.of.samples   = metastatic.breast.cancer.microarray.expr.matrix[,metastatic.breast.cancer.Basal.sample],expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.Her2.correlation.result  <- pick.out.cell.line(expr.of.samples   = metastatic.breast.cancer.microarray.expr.matrix[,metastatic.breast.cancer.Her2.sample],expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)
microarray.vs.CCLE.LumA.correlation.result  <- pick.out.cell.line(expr.of.samples   = metastatic.breast.cancer.microarray.expr.matrix[,metastatic.breast.cancer.LumA.sample],expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene.1000)



####### Save the results #########
save(file='client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData',
     list=c('microarray.vs.CCLE.correlation.result',
            'microarray.vs.CCLE.liver.correlation.result','microarray.vs.CCLE.lung.correlation.result','microarray.vs.CCLE.bone.correlation.result','microarray.vs.CCLE.brain.correlation.result','microarray.vs.CCLE.lymph.correlation.result',
            'microarray.vs.CCLE.LumB.correlation.result','microarray.vs.CCLE.Basal.correlation.result','microarray.vs.CCLE.Her2.correlation.result','microarray.vs.CCLE.LumA.correlation.result'
            )
     )








