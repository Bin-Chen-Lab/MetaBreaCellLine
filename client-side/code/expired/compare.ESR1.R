load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')

MET500.breast.cancer.polyA.non.Basal.sample           <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample               <- c(MET500.breast.cancer.polyA.Basal.sample)

ESR1.profile                   <- MET500.breast.cancer.gene.mutation.profile['ESR1',]
esr1.mutated.sample.MET500.id  <- names(ESR1.profile)[ESR1.profile==1]
esr1.wildtype.sample.MET500.id <- names(ESR1.profile)[ESR1.profile==0]

flag <- match(MET500.breast.cancer.polyA.non.Basal.sample,MET500.sample.meta$Run)
MET500.breast.cancer.polyA.non.Basal.sample.MET500.id <- MET500.sample.meta$MET500.id[flag]
esr1.mutated.sample  <- MET500.breast.cancer.polyA.non.Basal.sample[MET500.breast.cancer.polyA.non.Basal.sample.MET500.id %in% esr1.mutated.sample.MET500.id]
esr1.wildtype.sample <- MET500.breast.cancer.polyA.non.Basal.sample[MET500.breast.cancer.polyA.non.Basal.sample.MET500.id %in% esr1.wildtype.sample.MET500.id]

expr.matrix <- MET500.log2.read.count.matrix[,c(esr1.mutated.sample,esr1.wildtype.sample)]
flag        <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
expr.matrix <- expr.matrix[flag == ncol(expr.matrix),]


require(DESeq2)
df <- data.frame(condition=c(rep(x='ESR1.mutated',times=length(esr1.mutated.sample)), rep(x='ESR1.wildtype',times=length(esr1.wildtype.sample))))
df$condition <- factor(df$condition,levels = c('ESR1.wildtype','ESR1.mutated'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$pvalue),]
de.gene <- rownames(res)[res$padj < 0.1 & abs(res$log2FoldChange) > 1]







