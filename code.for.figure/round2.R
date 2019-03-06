require(plyr)
require(dplyr)
require(pheatmap)
require(ggplot2)
library(ComplexHeatmap)
require(foreach)
source('client-side/code/for.figure/ggplot.style.R')
load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
load('~/Project/Cancer2CellLine/client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output/for.revision.nature.communications.round1.R.output//for.revision.nature.communications.round1.RData')



##############################################################################################################################
####### Below, I addressed the issues proposed by reviewer 3 in the second round of review
##############################################################################################################################

######################################################################
#
# Question 1: No computation is needed here
#
######################################################################



######################################################################
#
# Question 2:  MET500.vs.Stage.IIA  MET500.vs.Stage.IIB MET500.vs.StageIIIA, results highly correlated, demonstrating that stage is not a severe confounding factor
# Subtype-specific analysis, hard to say for Basal-like due to small number of samples.
#
######################################################################
# re-compute differential mutation with the TCGA mutation data downloaded from April
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
common.gene.list     <- intersect(names(TCGA.breast.cancer.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- TCGA.breast.cancer.gene.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(TCGA.breast.cancer.gene.mutation.freq[TCGA.breast.cancer.gene.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.TCGA.raw.dm.df <- dm.df


load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
require(cgdsr)
study.id   <- 'brca_tcga'
mycgds     <-  CGDS("http://www.cbioportal.org/")

case.list  <- 'brca_tcga_sequenced'
profile.id <- 'brca_tcga_mutations'
clinical.data                           <- getClinicalData(mycgds,case.list)
Breast.Invasive.Ductal.Carcinoma.sample <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']
stage.vec                               <- clinical.data[Breast.Invasive.Ductal.Carcinoma.sample,'AJCC_PATHOLOGIC_TUMOR_STAGE']
names(stage.vec)                        <- Breast.Invasive.Ductal.Carcinoma.sample

stage.mutation.freq.list <- foreach(stage = unique(stage.vec)) %do% {
  sample.id <- names(stage.vec)[stage.vec == stage]
  gene.mutation.freq <- apply(TCGA.breast.cancer.gene.mutation.profile[,sample.id],1,function(x) sum(x)/length(x))
  
  gene.mutation.freq
}
names(stage.mutation.freq.list) <- unique(stage.vec)

######## MET500.vs.Stage.II.A TCGA samples
common.gene.list     <- intersect(names(TCGA.breast.cancer.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 

Stage.IIA.mutation.freq <-  stage.mutation.freq.list[['Stage IIA']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)

dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {  
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIA.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIA.mutation.freq[Stage.IIA.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIA.dm.df <- dm.df

#plot(x=-1 * MET500.vs.Stage.IIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.GDC.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
plot(x=-1 * MET500.vs.Stage.IIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIA.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIA.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIA)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIA.pdf',width=20,height=20)







Stage.IIB.mutation.freq <-  stage.mutation.freq.list[['Stage IIB']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
#common.gene.list     <- intersect(names(cBioPortal.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIB.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIB.mutation.freq[Stage.IIB.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F) #right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIB.dm.df <- dm.df
plot(x=-1 * MET500.vs.Stage.IIB.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIB.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIB.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIB.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIB)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIB.pdf',width=20,height=20)




Stage.IIIA.mutation.freq <-  stage.mutation.freq.list[['Stage IIIA']] 
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
#common.gene.list     <- intersect(names(cBioPortal.TCGA.brca.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) 
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
  q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
  prob    <- Stage.IIIA.mutation.freq[g]    
  if(prob ==0) {
    prob <- min(Stage.IIIA.mutation.freq[Stage.IIIA.mutation.freq > 0])  
  }
  p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
  data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]
MET500.vs.Stage.IIIA.dm.df <- dm.df
plot(x=-1 * MET500.vs.Stage.IIIA.dm.df$p.value %>% log10,y=-1 * MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIIA.dm.df),'p.value'] %>% log10,xlim=c(0,11),ylim=c(0,11))
lines(c(0,11),c(0,11))
draw.df <- data.frame( x =-1 * (MET500.vs.Stage.IIIA.dm.df$p.value %>% log10),y= -1 * (MET500.vs.TCGA.raw.dm.df[rownames(MET500.vs.Stage.IIIA.dm.df),'p.value'] %>% log10))
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=5) + ggplot.style + xlim(c(0,13)) + ylim(c(0,13)) + geom_abline(slope = 1,intercept = 0) + xlab('MET500.vs.TCGA (Stage IIIA)') + ylab('MET500.vs.TCGA')
cor(draw.df$x,draw.df$y,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/nature.communications.rebutal.letter.round.2/p.value.StageIIIA.pdf',width=20,height=20)



source('client-side/code/for.revision/nature.communicaionts/round2.DE.R') # perform DE analysis with RUVg inferred factor values
