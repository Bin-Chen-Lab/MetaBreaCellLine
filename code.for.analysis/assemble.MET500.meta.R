#This R code file is used to generate a csv file which contains comprehensive meta information of all MET500 samples.
#The final output is stored in the file client-side/meta.data/MET500.RNASeq.meta.csv



require(plyr)
require(dplyr)

MET500.sra.run                      <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/MET500.SraRunTable.txt", stringsAsFactors=FALSE)
MET500.sra.run                      <- MET500.sra.run[,c('Assay_Type','LibrarySelection','LibrarySource','Run','is_tumor','body_site','submitted_subject_id')]
MET500.sra.run$dbGap.subject.id     <- MET500.sra.run$submitted_subject_id
MET500.sra.run$submitted_subject_id <- NULL

MET500.to.dbGap                     <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500_dbGaP_SubjectID.csv")

MET500.sample.phenotype                    <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500.sample.phenotype.csv", stringsAsFactors=FALSE)
MET500.sample.phenotype$primary.tumor.site <- MET500.sample.phenotype$primary.site
MET500.sample.phenotype$primary.tumor.site <- toupper(MET500.sample.phenotype$primary.tumor.site)
MET500.sample.phenotype$primary.site       <- NULL
MET500.sample.phenotype                    <- MET500.sample.phenotype[,c('MET500.id','cancer.type','primary.tumor.site')]

MET500.meta             <- merge(x=MET500.sra.run,y=MET500.to.dbGap,by = 'dbGap.subject.id')
MET500.meta$biopsy.site <- toupper(MET500.meta$body_site)
MET500.meta$body_site   <- NULL
MET500.meta             <- merge(MET500.meta,MET500.sample.phenotype,by='MET500.id')

MET500.RNASeq.meta      <- MET500.meta[MET500.meta$Assay_Type =='RNA-Seq',]


write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$cancer.type == 'Prostate Adenocarcinoma' & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Prostate Adenocarcinoma.txt',
          quote=F,row.names=F,col.names=F
          )

write.csv(x=MET500.RNASeq.meta$Run[grepl(x=MET500.RNASeq.meta$cancer.type,pattern = 'Breast')  & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Breast.cancer.txt',
          quote=F,row.names=F,col.names=F
)


write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$cancer.type == 'Lung Adenocarcinoma' & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Lung Adenocarcinoma.txt',
          quote=F,row.names=F,col.names=F
)

write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$cancer.type == 'Cutaneous Melanoma' & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Cutaneous Melanoma.txt',
          quote=F,row.names=F,col.names=F
)

write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$cancer.type == 'Prostate Neuroendocrine Carcinoma' & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Prostate Neuroendocrine Carcinoma.txt',
          quote=F,row.names=F,col.names=F
)


write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$cancer.type == 'Cutaneous Melanoma' & MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/Cutaneous Melanoma.txt',
          quote=F,row.names=F,col.names=F
)

write.csv(x=MET500.RNASeq.meta$Run[MET500.RNASeq.meta$is_tumor == 'Yes'],
          file ='client-side/MET500.exp.list/MET500.SRA.txt',
          quote=F,row.names=F,col.names=F
)

write.csv(x=MET500.RNASeq.meta,file='client-side/meta.data/MET500.RNASeq.meta.csv',quote=TRUE,row.names=FALSE,col.names=TRUE)
