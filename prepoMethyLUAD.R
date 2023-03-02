setwd("/home/jvasquez/Documents/TCGA-lung")
library(TCGAbiolinks)
library(data.table)

subtypeLUAD=read.table("subtypeLUAD.tsv",header=T,sep='\t')
#get the data
# mthyltn2 <-  GDCquery(project = "TCGA-LUAD",
#                      data.category = "DNA Methylation",
#                      platform="Illumina Human Methylation 450",
#                      barcode=subtypeLUAD$samples)
# Error in GDCdownload(mthyltn) : 
#   We can only download one data type. Please use data.type argument in GDCquery to filter results.

mthyltn <-  GDCquery(project = "TCGA-LUAD",
                     data.category = "DNA Methylation",
                     data.type = "Methylation Beta Value",
                     platform="Illumina Human Methylation 450",
                     barcode=subtypeLUAD$samples)

# mthyltn <-  GDCquery(project = "TCGA-LUAD",
#                      data.category = "DNA Methylation",
#                      data.type = "Masked Intensities",
#                      platform="Illumina Human Methylation 450",
#                      barcode=subtypeLUAD$samples)
# In addition: Warning messages:
#1: In fread(x, select = 2, stringsAsFactors = F) :
#  Previous fread() session was not cleaned up properly. Cleaned up ok at the beginning of this fread() call.
#2: In fread(x, select = 2, stringsAsFactors = F) :

GDCdownload(mthyltn)
files=list.files("GDCdata",full.names=T,recursive=T)
#files=list.files("GDCdata",full.names=T,recursive=T)
#GDCdownload will download 197 files. A total of X GB
methy=do.call(cbind,pbapply::pbsapply(files,function(x) 
  fread(x,select=2,stringsAsFactors=F)))

library(sesameData)
library(sesame)
mthyltn2=GDCprepare(mthyltn)
methy <- mthyltn2@assays@data[[1]]
dim(methy)
#[1] 485577    197
write.table(methy,"methy.tsv",sep='\t',quote=F)

#
# methy_2=BiocGenerics::do.call(cbind,pbapply::pbsapply(files,function(x) 
#   fread(x,select=2,stringsAsFactors=F)))
#######drop probes with too many na########################################
#before noise-prone filter or the process will be slooow

#subtype to duplicate
i=substr(colnames(methy),1,19)
j=i[duplicated(i)]
designMethy=subtypeLUAD[c(which(!subtypeLUAD$samples%in%j),
                      as.numeric(sapply(which(subtypeLUAD$samples%in%j),rep,2))),1:4]
#needed coz names are not equal to expression data but barcodes do
designMethy$barcode=unlist(sapply(designMethy$samples,function(x)
  colnames(methy)[i==x][1]))
designMethy$barcode[designMethy$samples%in%j]=rev(colnames(methy)[
  which(i%in%j)])
designMethy=designMethy[order(match(designMethy$barcode,
                                    colnames(methy))),]
total=table(subtypeLUAD$subtype)
#NA per subtype
nas=lapply(names(total),function(x) 
  rowSums(is.na(methy[,designMethy$subtype==x])))
#keep probes with NA in less than 25% of samples of all subtypes
i=unique(unlist(lapply(1:4,function(x) which(nas[[x]]<total[x]*.25))))
methy=methy[i,]
dim(methy)
#[1] 417031    844

############filter noise-prone probes##############################
#get probes annotation
annot=fread(files[1])
colnames(annot)[1]="probe"
annot=annot[annot$probe%in%rownames(methy),]
annot=annot[order(match(annot$probe,rownames(methy))),]
#sum(annot$probe==rownames(methy)) #order is the same

#sex chrs should be dropped if mixed gender
#since this is not the case, keep chrX
#https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
#methy=methy[!annot$Chromosome%in%c("chrX","chrY"),]
#dim(methy)

#drop probes with ambigous mapping
#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
methy=methy[annot$Chromosome!="*",]
nrow(methy)
#[1] 393197

#polymorphisms may affect DNAm measurements
methy_grep=methy[grep("rs",rownames(methy),invert=T),]
nrow(methy_grep)
#[1] 393132

