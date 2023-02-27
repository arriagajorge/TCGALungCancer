setwd("/home/jvasquez/Documents/TCGA-lung")
library(TCGAbiolinks)
library(biomaRt)

subtypeLUAD=read.table("subtypeLUAD.tsv",header=T,sep='\t')
#get the data
mirnas <- GDCquery(project = "TCGA-LUAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification",
                   barcode=subtypeLUAD$samples)

#Genome of reference: hg38
#https://api.gdc.cancer.gov/data/683367de-81c9-408c-85fd-2391f3e537ee
#says miRBase v.21 was used for harmonization annotation
GDCdownload(mirnas)
mir=GDCprepare(mirnas)
rownames(mir)=mir$miRNA_ID
mir=mir[,grep("read_count",colnames(mir))]
colnames(mir)=gsub("read_count_","",colnames(mir))
dim(mir)
#[1] 1881  196
write.table(mir,"miRNAseq.tsv",sep='\t',quote=F)

#subtype to duplicates
i=substr(colnames(mir),1,19)
j=i[duplicated(i)]
designExp=subtypeLUAD[c(which(!subtypeLUAD$samples%in%j),
                    as.numeric(sapply(which(subtypeLUAD$samples%in%j),rep,2))),]
designExp=designExp[order(match(designExp$samples,substr(colnames(mir),1,19))),]
designExp$barcode=colnames(mir)

#biomart version with all the miRNAs & the same hg
mart=useMart("ensembl",host="https://may2017.archive.ensembl.org", 
             dataset = "hsapiens_gene_ensembl")
mart2=useMart("ensembl", dataset = "hsapiens_gene_ensembl")

myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "mirbase_id",
                             "start_position","end_position"),
              mart=mart2)
myannot=myannot[myannot$mirbase_id%in%rownames(mir),]
#there should not be a length bias
myannot$length=abs(myannot$end_position-myannot$start_position)
#discard duplicated entries with the same %CpG
# > sum(!duplicated(myannot[,2:3]))
# [1] 1863
myannot=myannot[!duplicated(myannot[,2:3]),]
#there're duplicates with slightly different %CpG & position
#temp=myannot[myannot$mirbase_id%in%myannot$mirbase_id[duplicated(myannot$mirbase_id)],]
# summary(colSums(sapply(unique(temp$mirbase_id),function(x)
# temp$percentage_gene_gc_content[temp$mirbase_id==x])*c(1,-1)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.46000 -0.61000  0.15000 -0.05105  0.37500  1.74000 
#choose 1 of the duplicates randomly
myannotAlt=myannot[duplicated(myannot$mirbase_id),]
myannot=myannot[!duplicated(myannot$mirbase_id),]
#if the GC bias aint fixed you CAN NOT compare among miRNAs

##################CHECK BIASES########################################################
library(NOISeq)

noiseqData = readData(data = mir, factor=designExp,
                      gc=myannot[,c(3,2)],length=myannot[,c(3,6)])
mycountsbio = dat(noiseqData, type = "countsbio",factor = "subtype")#check low counts
#distributions
png("miROri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:4)
dev.off()
#counts
png("miRcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:4)
dev.off()
png("miRlowCountThres.png")
hist(rowMeans(edgeR::cpm(mir,log=T)),ylab="miRNA",
     xlab="mean of log CPM",col="gray",xlim=c(-5,20))
dev.off()

#check length & GC bias
myGCcontent <- dat(noiseqData, k = 0, type = "GCbias",
                   factor = "subtype")
png("miRGCbiasOri.png",width=1000)
par(mfrow=c(1,4))
sapply(1:4,function(x) explo.plot(myGCcontent, samples = x))
dev.off()
mylenBias <- dat(noiseqData, k = 0, type = "lengthbias",
                 factor = "subtype")
png("lengthbiasOri.png",width=1000)
par(mfrow=c(1,5))
sapply(1:4,function(x) explo.plot(mylenBias, samples = x))
dev.off()
#no GC bias nor lengthbias!!!!!!!!

myPCA = dat(noiseqData, type = "PCA", norm = FALSE, 
            logtransf = FALSE)#check batches
png("miRPCA_Ori.png")
explo.plot(myPCA, samples = c(1,2), plottype = "scores", 
           factor = "subtype")
dev.off()
mycd = dat(noiseqData, type = "cd", norm = FALSE)#check if normalizations is needed
table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Warning: 368 features with 0 counts in all samples are to be removed for this analysis."
#FAILED PASSED 
#   194      1
png("miRcdOri.png")
explo.plot(mycd,samples=1:10)
dev.off()
#################SOLVE BIASES######################################################
#filter low counts
FilteredMatrix = filtered.data(mir, factor = "subtype",
                               norm = FALSE, method = 1, cpm = 0)
#206 features are to be kept for differential expression analysis with filtering method 1
#it is expected that in miRNA-seq experiments, the 75th percentile 
#of the data will be found at only 1 or 2 copies/library [10.1093/bib/bbv019]
#Drago-García2017 used a minimum of  5  counts  in  at  least  25%
# of  the  samples 
temp=lapply(unique(designExp$subtype),function(x)   #this was commented
 mir[,colnames(mir)%in%designExp$barcode[designExp$subtype==x]])
temp1=names(which(table(unlist(sapply(temp,function(x) rownames(x)[rowSums(x>=5)>=ncol(x)*.25])))==5))
length(temp1)
# [1] 328
# >1 copy → 594
FilteredMatrixAlt=mir[rowSums(mir)>0,]

#TMM, UQ, median & DESEq are similar [10.1186/gb-2010-11-3-r25]
#TMM, UQ are the best [10.1093/bib/bbv019] 
myTMM=tmm(FilteredMatrix,lc=0)
#myTMMAlt=tmm(FilteredMatrixAlt,lc=0)
noiseqData = readData(data = myTMM, factors=designExp)
#noiseqData = readData(data = myTMMAlt, factors=designExp)
mycdTMM = dat(noiseqData, type = "cd", norm = T)
#mycdTMMAlt = dat(noiseqData, type = "cd", norm = T)
table(mycdTMM@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED     #sometimes there are changes in the values for example 4 to 2 in failed
#   4     191
#table(mycdTMMAlt@dat$DiagnosticTest[,  "Diagnostic Test"])
#[1] "Diagnostic test: PASSED."
#explo.plot(mycdTMMAlt,samples=1:10)#non-comparable samples at plot

myUQ=uqua(FilteredMatrix,lc=0)
noiseqData = readData(data = myUQ, factors=designExp)
mycdUQ = dat(noiseqData, type = "cd", norm = T)
table(mycdUQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    17    178

library(EDASeq)
mydataEDA <- newSeqExpressionSet(
  counts=as.matrix(FilteredMatrix),
  phenoData=data.frame(designExp,row.names=designExp$barcode))
norm.counts <- betweenLaneNormalization(mydataEDA,
                                        which = "median", offset = FALSE)
noiseqData = NOISeq::readData(data = assayData(norm.counts)$normalizedCounts,
                              factors=designExp)
mycdMedian = dat(noiseqData, type = "cd", norm = T)
table(mycdMedian@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#    30    165
# what happened with LUSC add or not its a sum?

library(DESeq2)
deseqFactors=estimateSizeFactors(newCountDataSet(FilteredMatrix,
                                                 conditions=designExp))
deseqFactors=DESeq(DESeqDataSetFromMatrix(FilteredMatrix, colData = designExp))

myDESEQ=counts(deseqFactors,normalized=T)
noiseqData = NOISeq::readData(data = myDESEQ, factors=designExp)
mycdDESEQ = NOISeq::dat(noiseqData, type = "cd", norm = T)
table(mycdDESEQ@dat$DiagnosticTest[,  "Diagnostic Test"])
#FAILED PASSED 
#   707    101 

png("miRcd_final.png")
explo.plot(mycdMedian,samples=1:10)
dev.off()
