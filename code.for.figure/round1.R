### code to addrress the issues in WORD file Reviewer comments.round.1.docx
require(plyr)
require(dplyr)
require(cgdsr)

load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
source('client-side/code/for.figure/ggplot.style.R')
source('client-side/code/util.R')

#### I borrow code from organize.breast.cancer.genomics.data.R, to do symbol mapping
hgnc.complete.set               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hgnc_complete_set.txt", stringsAsFactors=FALSE)
hgnc.complete.set               <- subset( hgnc.complete.set,select=c('hgnc_id','symbol'))
colnames(hgnc.complete.set)     <- c('HGNC.ID','Symbol')
df                              <- read.delim("~/Project/Cancer2CellLine/client-side/Data/somatic.mutation/CCLE.captured.gene.txt", stringsAsFactors=FALSE)
df                              <- df[df$HGNC.Symbol != '',]
idx                             <- match(x = df$HGNC.ID,table = hgnc.complete.set$HGNC.ID)
mapping.df                      <- data.frame(previous.hgnc.symbol=df$HGNC.Symbol,new.hgnc.symbol=hgnc.complete.set$Symbol[idx]) %>% unique # There are duplicates in the CCLE.captured.gene.txt: two genes correspond to HGNC symbol MECOM
rownames(mapping.df)            <- mapping.df$previous.hgnc.symbol %>% as.character
mapping.df$previous.hgnc.symbol <- as.character(mapping.df$previous.hgnc.symbol)
mapping.df$new.hgnc.symbol      <- as.character(mapping.df$new.hgnc.symbol) # The mapping.df is used to map all HGNC symbols to newest version


hg19.RefSeq.gene.coordinates <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hg19.RefSeq.gene.coordinates.txt", stringsAsFactors=FALSE)
flag                         <- grepl(x=hg19.RefSeq.gene.coordinates$chrom,pattern = '_')
hg19.RefSeq.gene.coordinates <- hg19.RefSeq.gene.coordinates[!flag,] #Let us remove the genes which are on the un-assembled contigs!
tmp                          <- ddply(hg19.RefSeq.gene.coordinates,.(name2),function(x) data.frame(chrom=x$chrom[1], start=min(x$txStart),end=max(x$txEnd),geneid=x$name2[1],genename=x$name2[1]))
tmp$name2                    <- NULL
hg19.gene.info               <- tmp
hg19.gene.info$chrom         <- gsub(x=hg19.gene.info$chrom,pattern='chr',replacement = '')
hg19.gene.info$chrom         <- as.character(hg19.gene.info$chrom)
hg19.gene.info$geneid        <- as.character(hg19.gene.info$geneid)
hg19.gene.info$genename      <- as.character(hg19.gene.info$genename)
flag                         <- which(hg19.gene.info$genename %in% mapping.df$previous.hgnc.symbol)
for(i in flag){
    pre.symbol <- hg19.gene.info$genename[i]  %>% as.character
    new.symbol <- mapping.df[pre.symbol,'new.hgnc.symbol']
    hg19.gene.info[i,'genename']        <- new.symbol
    hg19.gene.info[i,'geneid']          <- new.symbol
} #In hg19 RefSeq release, some CCLE captured genes are still with expired HGNC symbol, just replace them

CCLE.captured.gene.hg19.info       <- hg19.gene.info[hg19.gene.info$genename %in% mapping.df$new.hgnc.symbol,] # Well, I trust RefSeq. If some gene symbols here were missed after this filtering, then remove them. Actually, only 9 were missed 
CCLE.captured.gene.hg19.info       <- CCLE.captured.gene.hg19.info[CCLE.captured.gene.hg19.info$chrom != 'Y',] # Well, we are studying breast cancer, remove genes on chromosome Y. 
CCLE.captured.gene                 <- CCLE.captured.gene.hg19.info$genename


############################################################### 
#   
#   Question 1: adjust batch effects of variant-calling pipeline
#
############################################################### 

###############################################################
######### organize variant calling results from GDC data portal
###############################################################
# study.id   <- 'brca_tcga'
# mycgds     <-  CGDS("http://www.cbioportal.org/")
# case.list  <- 'brca_tcga_sequenced'
# profile.id <- 'brca_tcga_mutations'
# clinical.data                           <- getClinicalData(mycgds,case.list)
# TCGA.BRCA.Ductal.Carcinoma.cohort.id    <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']

TCGA.BRCA.Ductal.Carcinoma.cohort.id <- colnames(TCGA.breast.cancer.gene.mutation.profile)
TCGA.BRCA.varscan.somatic            <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf", stringsAsFactors=FALSE)
TCGA.BRCA.varscan.somatic            <- TCGA.BRCA.varscan.somatic[,c('Hugo_Symbol','Variant_Classification','Tumor_Sample_Barcode')]
extract.cohort.id <- function(x) {
    tmp <- strsplit(x,split = '\\-') %>% unlist  
    sprintf("%s.%s.%s.01",tmp[1],tmp[2],tmp[3])
}

TCGA.BRCA.varscan.somatic$cohort.id <- sapply(TCGA.BRCA.varscan.somatic$Tumor_Sample_Barcode,extract.cohort.id)
flag                                <- which(TCGA.BRCA.varscan.somatic$Hugo_Symbol %in% mapping.df$previous.hgnc.symbol)
for(i in flag){
    pre.symbol <- TCGA.BRCA.varscan.somatic$Hugo_Symbol[i]  %>% as.character
    new.symbol <- mapping.df[pre.symbol,'new.hgnc.symbol']
    TCGA.BRCA.varscan.somatic[i,'Hugo_Symbol']        <- new.symbol
} # map previous symbols to newest symbol


GDC.TCGA.brca.gene.mutation.profile <- matrix(data=0,nrow = length(CCLE.captured.gene),ncol=length(TCGA.BRCA.Ductal.Carcinoma.cohort.id))
rownames(GDC.TCGA.brca.gene.mutation.profile) <- CCLE.captured.gene
colnames(GDC.TCGA.brca.gene.mutation.profile) <- TCGA.BRCA.Ductal.Carcinoma.cohort.id


flag <- (TCGA.BRCA.varscan.somatic$cohort.id %in% TCGA.BRCA.Ductal.Carcinoma.cohort.id) & (TCGA.BRCA.varscan.somatic$Hugo_Symbol %in% CCLE.captured.gene)
tmp  <- TCGA.BRCA.varscan.somatic[flag,]
for (i in 1:nrow(tmp)) {
    gene      <- tmp[i,'Hugo_Symbol']      %>% as.character
    sample.id <- tmp[i,'cohort.id']        %>% as.character
    GDC.TCGA.brca.gene.mutation.profile[gene,sample.id] <- 1
}
GDC.TCGA.brca.gene.mutation.freq <- apply(GDC.TCGA.brca.gene.mutation.profile,1,sum) / ncol(GDC.TCGA.brca.gene.mutation.profile)
save(file='client-side/output/for.revision.nature.communications.round1.R.output/for.revision.nature.communications.round1.RData',list=c('GDC.TCGA.brca.gene.mutation.freq','GDC.TCGA.brca.gene.mutation.profile'))


###############################################################
######### show the mutation frequency are consistent, Fig R1, panel a
###############################################################
plot(x=TCGA.breast.cancer.gene.mutation.freq,GDC.TCGA.brca.gene.mutation.freq[names(TCGA.breast.cancer.gene.mutation.freq)],xlab='cBioPortal',ylab='GDC',cex.lab=1.5,cex.axis=1.5)
draw.df <- data.frame(cBioPortal = TCGA.breast.cancer.gene.mutation.freq,GDC = GDC.TCGA.brca.gene.mutation.freq[names(TCGA.breast.cancer.gene.mutation.freq)])
ggplot(draw.df,aes(x=cBioPortal,y=GDC)) + geom_point(size=5) + ggplot.style + xlim(c(0,0.5)) + ylim(c(0,0.5)) + geom_abline(slope = 1,intercept = 0)
#ggplot(draw.df,aes(x=cBioPortal,y=GDC)) + geom_point(size=5) + ggplot.style + xlim(c(0,0.15)) + ylim(c(0,0.15)) + geom_abline(slope = 1,intercept = 0)
cor(draw.df$cBioPortal,draw.df$GDC,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/cBioportal.vs.GDC.mutation.frequency.pdf',width=20,height=20)

###############################################################
######### show the p-values are consistent
###############################################################
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
common.gene.list     <- intersect(rownames(GDC.TCGA.brca.gene.mutation.profile),rownames(MET500.breast.cancer.gene.mutation.profile)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
    q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
    prob    <- GDC.TCGA.brca.gene.mutation.freq[g]    
    if(prob==0) {
        prob <- 1/ncol(GDC.TCGA.brca.gene.mutation.profile)  
    }
    p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = FALSE)#right side p-value,metastatic cancer has more mutation burdens
    data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
MET500.vs.GDC.dm.df <- dm.df



MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
common.gene.list     <- intersect(rownames(TCGA.breast.cancer.gene.mutation.profile),rownames(MET500.breast.cancer.gene.mutation.profile)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- TCGA.breast.cancer.gene.mutation.freq[g]    
  if(prob==0) {
    prob <- 1/ncol(TCGA.breast.cancer.gene.mutation.profile)  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
MET500.vs.cBioPortal.dm.df <- dm.df

draw.df <- data.frame(cBioPortal=-1 * (MET500.vs.cBioPortal.dm.df$p.value %>% log10),GDC= -1 * (MET500.vs.GDC.dm.df[rownames(MET500.vs.cBioPortal.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=cBioPortal,y=GDC)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.cBioPortal') + ylab('MET500.vs.GDC')
cor(draw.df$cBioPortal,draw.df$GDC,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/p.value.consistency.pdf',width=20,height=20)



############################################################### 
#   
#   Question 2: just to say yes,all of the genes are within CCLE captured genes
#
############################################################### 



############################################################### 
#   
#   Question 3: show that CCLE-processed gene expression profiles are highly correlated with the porfile processed by us (regarding to the most varied 1000 genes)
#   Fig S11a S11b
#
############################################################### 
load('server-side/RData/CCLE_BRACA.RData')
colnames(BRACA.log2.fpkm.matrix) <- paste(colnames(BRACA.log2.fpkm.matrix),'_BREAST',sep='')

load('server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


common.genes      <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(BRACA.log2.fpkm.matrix))
common.cell.lines <- intersect(colnames(CCLE.log2.rpkm.matrix),colnames(BRACA.log2.fpkm.matrix))
cor.vec <- foreach(cell.line = common.cell.lines,.combine='c') %do% {
    plot(x=CCLE.log2.rpkm.matrix[common.genes,cell.line],y=BRACA.log2.fpkm.matrix[common.genes,cell.line])  
    cor(x=CCLE.log2.rpkm.matrix[common.genes,cell.line],y=BRACA.log2.fpkm.matrix[common.genes,cell.line],method='spearman')  
}
names(cor.vec) <- common.cell.lines

ggplot(data.frame(x=cor.vec),aes(x=x)) + geom_density(size=3) + ggplot.style + xlab('spearman rank correlation') + ylab('') + xlim(c(0,1)) + geom_vline(xintercept = median(cor.vec),linetype='dashed')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/spearman.rank.cor.distribution.pdf',width=25,height=20)


draw.df <- data.frame(x=CCLE.log2.rpkm.matrix[common.genes,'MCF7_BREAST'],y=BRACA.log2.fpkm.matrix[common.genes,'MCF7_BREAST'])  
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,15)) + ylim(c(0,15)) + geom_abline(slope = 1,intercept = 0) + xlab('CCLE pipeline') + ylab('RSEM pipeline')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/MCF7.pdf',width=25,height=20)
cor(x=CCLE.log2.rpkm.matrix[common.genes,'MCF7_BREAST'],y=BRACA.log2.fpkm.matrix[common.genes,'MCF7_BREAST'],method='spearman')  



load('server-side/RData/MET500.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')

tmp <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.sample],expr.of.cell.lines = BRACA.log2.fpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
x <- tmp$cell.line.median.cor

load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')

y <- MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor
draw.df <- data.frame(x=x,y=y[names(x)])
draw.df <- draw.df[complete.cases(draw.df),]
ggplot(draw.df,aes(x=rank(y),y=rank(x))) + geom_point(size=4)  + ggplot.style + geom_abline(slope=1,intercept = 0) + xlab('CCLE pipeline') + ylab('RSEM pipeline')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/TC.correlation.pdf',width=20,height=20)
cor(x,y[names(x)],method = 'spearman',use = 'pairwise.complete.obs')

############################################################### 
#   
#   Question 4: this is a open question, no golden answer
#
############################################################### 


############################################################### 
#   
#   Question 5: CCLE vs TCGA , TC analysis
#
############################################################### 
load('server-side/RData/CCLE.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('~/Project/CancerImmunoTherapy-TCGA/server-side//RData//Breast Invasive Carcinoma.RData')
load("~/Project/Cancer2CellLine/client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData")
TCGA.breast.cancer.log2.fpkm.matrix <- log2.fpkm.matrix

#should extract Breast Ductal Carcinoma sample here
study.id   <- 'brca_tcga'
mycgds     <-  CGDS("http://www.cbioportal.org/public-portal/")
case.list  <- 'brca_tcga_rna_seq_v2_mrna'
clinical.data                            <- getClinicalData(mycgds,case.list)
Breast.Invasive.Ductal.Carcinoma.sample  <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']
Breast.Invasive.Ductal.Carcinoma.sample  <- gsub(x = Breast.Invasive.Ductal.Carcinoma.sample,replacement = '-',pattern = '\\.')
Breast.Invasive.Ductal.Carcinoma.sample  <- intersect(Breast.Invasive.Ductal.Carcinoma.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))


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



TCGA.vs.CCLE.polyA.expression.correlation.result  <- pick.out.cell.line(TCGA.breast.cancer.log2.fpkm.matrix[,Breast.Invasive.Ductal.Carcinoma.sample], CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)

draw.df <- data.frame(y=MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],x=TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line])
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,0.5)) + ylim(c(0,0.5)) + geom_abline(slope = 1,intercept = 0) + xlab('TCGA.vs.CCLE') + ylab('MET500.vs.CCLE')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/TCGA.and.MET500.both.vs.CCLE.pdf',width=20,height=20)

cor(draw.df$y,draw.df$x,method='spearman')




############################################################### 
#   
#   Question 6,7,8: hmm, just revised some text and tables, no more computations needed
#
############################################################### 



############################################################### 
#   
#   Question 9:
#
############################################################### 


load('client-side/output/ssgsea.analysis.R.output/ssgsea.analysis.RData')
EMT.ssgsea.score <- hallmark.geneset.gsea.results['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]
y                <- MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor
draw.df          <- data.frame(x=EMT.ssgsea.score[CCLE.breast.cancer.cell.line],y=y[CCLE.breast.cancer.cell.line],cell.line.name=CCLE.breast.cancer.cell.line)
draw.df$is.basal <- draw.df$cell.line.name %in% CCLE.breast.cancer.Basal.cell.line
ggplot(draw.df) + geom_point(aes(x=x,y=y,label=cell.line.name,color=is.basal),size=6) + ggplot.style + xlab('ssGSEA score') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/similarity.vs.EMT.pdf',width=20,height=20)


EMT.ssgsea.score <- hallmark.geneset.gsea.results['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]
y                <- MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor
draw.df          <- data.frame(x=EMT.ssgsea.score[CCLE.breast.cancer.cell.line],y=y[CCLE.breast.cancer.cell.line],cell.line.name=CCLE.breast.cancer.cell.line)
draw.df$is.basal <- draw.df$cell.line.name %in% CCLE.breast.cancer.Basal.cell.line
ggplot(draw.df) + geom_point(aes(x=x,y=y,label=cell.line.name,color=is.basal),size=6) + ggplot.style + xlab('ssGSEA score') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/similarity.vs.EMT.basal.pdf',width=20,height=20)



load('~/Project/CancerImmunoTherapy-TCGA/client-side/output/prepare.meta.data.R.output/prepare.meta.data.RData') # I steal something from another project
flag <- rownames(BRACA.log2.fpkm.matrix) %in% names(ensemble2entrez)
data.for.ssgsea <- BRACA.log2.fpkm.matrix[flag,]
rownames(data.for.ssgsea) <- ensemble2entrez[rownames(data.for.ssgsea)]
require(GSVA)
score <- gsva(expr = data.for.ssgsea,gset.idx.list = genesets['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'],method='ssgsea')

plot(x=score,y=EMT.ssgsea.score[colnames(score)])

###########################################
#
# Question 10: MDAMB-231, more computation. Show that krt14 is lowly expressed in MDAMB231
#
###########################################
load('server-side/RData/MDAMB231.expression.RData')
MDAMB231.log2.fpkm.matrix <- log2.fpkm.matrix
MDAMB231.log2.fpkm.matrix <- cbind(MDAMB231.log2.fpkm.matrix,BRACA.log2.fpkm.matrix[,'MDAMB231_BREAST'])
krt14 <- 'ENSG00000186847'
#krt5  <- 'ENSG00000186081'
load('server-side/RData/MET500.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')

LumA.krt14.expr  <- MET500.log2.fpkm.matrix[krt14,MET500.breast.cancer.polyA.LumA.sample]
df1 <- data.frame(expr=LumA.krt14.expr,subtype='Luminal A')

LumB.krt14.expr  <- MET500.log2.fpkm.matrix[krt14,MET500.breast.cancer.polyA.LumB.sample]
df2 <- data.frame(expr=LumB.krt14.expr,subtype='Luminal B')

Her2.krt14.expr  <- MET500.log2.fpkm.matrix[krt14,MET500.breast.cancer.polyA.Her2.sample]
df3 <- data.frame(expr=Her2.krt14.expr,subtype='Her2-enriched')

Basal.krt14.expr <- MET500.log2.fpkm.matrix[krt14,MET500.breast.cancer.polyA.Basal.sample]
df4 <- data.frame(expr=Basal.krt14.expr,subtype='Basal-like')

MDAMB231.krt14.expr <- MDAMB231.log2.fpkm.matrix[krt14,]
df5 <- data.frame(expr=MDAMB231.krt14.expr,subtype='MDA-MB-231')



draw.df <- rbind(df1,df2,df3,df4,df5)
ggplot(draw.df) + geom_boxplot(aes(x=subtype,y=expr)) + ggplot.style + geom_point(aes(x=subtype,y=expr),size=3) + ylab('Expression') + xlab('') + ylim(c(0,15))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/krt14.pdf',width=20,height=20)
ggplot(draw.df) + geom_boxplot(aes(x=subtype,y=expr)) + ggplot.style + geom_point(aes(x=subtype,y=expr),size=3) + ylab('Expression') + xlab('') + coord_flip() + ylim(c(0,15))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/krt14.for.reviewer.pdf',width=20,height=20)

load('server-side/RData/MDAMB231.expression.RData')
MDAMB231.log2.fpkm.matrix <- log2.fpkm.matrix
y <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample],expr.of.cell.lines = cbind(MDAMB231.log2.fpkm.matrix,BRACA.log2.fpkm.matrix),marker.gene = CCLE.rna.seq.marker.gene.1000)

df <- data.frame(value=y$cell.line.median.cor,rank=rank(y$cell.line.median.cor))
flag <- grepl(rownames(df),pattern = 'SRR') | grepl(rownames(df),pattern = 'ERR') | grepl(rownames(df),pattern = 'MDAMB231')
df$is.tso <- 'NO'
df[flag,'is.tso'] <- 'YES'
ggplot() + geom_point(data=df,aes(x=rank,y=value),size=4) + ggplot.style + xlab('rank') + ylab('transcriptome similarity') + geom_point(data=df[flag,],aes(x=rank,y=value),color='red',size=7)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.1/Cell.line.ranking.with.Basal.and.additional.MDAMB231.pdf',width=20,height=20)




