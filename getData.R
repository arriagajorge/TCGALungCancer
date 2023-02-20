# @title getData_lungCancer
# @date 02/2023
# @author jvasquez

library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)
library(data.table)


mthyltnLUSC <-  GDCquery(project = "TCGA-LUSC",
	data.category = "DNA Methylation",
  platform="Illumina Human Methylation 450")

df1 <- mthyltnLUSC[[1]][[1]]
mthyltnLUSC = getResults(mthyltnLUSC)

i=substr(mthyltnLUSC$cases,1,19)

xprssnLUSC <- GDCquery(project = "TCGA-LUSC",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type="STAR - Counts")#
xprssnLUSC=getResults(xprssnLUSC)

j=substr(xprssnLUSC$cases,1,19)
mirnasLUSC <- GDCquery(project = "TCGA-LUSC",
                   data.category = "Transcriptome Profiling",
                   data.type = "miRNA Expression Quantification")
mirnasLUSC=getResults(mirnasLUSC)
k=substr(mirnasLUSC$cases,1,19)

##############CONCOURRENT MEASURES########################
sapply(list(i,j,k),function(x) length(unique(x)))
#[1] 412 553 523 
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 368    

samples_df=data.frame(cbind(samples,
                            sapply(samples,function(x) 
                              unique(as.character(mthyltnLUSC$sample_type[i==x])))))
colnames(samples_df)[2]="tissue"
# unique(samples_df$tissue)
# [1] "Primary Tumor"       "Solid Tissue Normal"
samples_df$patient=substr(samples_df$samples,1,12) 

samples <- samples_df

subtypes=TCGAquery_subtype(tumor="LUSC")#subtype per patient
sum(samples$patient%in%subtypes$patient)

#only classified samples  are useful
samples=samples[samples$patient%in%subtypes$patient,]
table(subtypes$Expression.Subtype[subtypes$patient%in%samples$patient])
# basal classical primitive secretory 
# 13        27        10        22 

#normal subtype is discarded NOT normal tissue 
samples_2<-samples[which(!samples$patient%in%subtypes$patient[   #samples_2 is only a name for dont overwrite
  subtypes$Expression.Subtype=="Normal"]|
    samples$tissue=="Solid Tissue Normal"),]

temp=table(samples_2[samples_2$patient%in%samples_2$patient[duplicated(
  samples_2$patient)],2:3])
table(apply(temp,2,paste,collapse=""))
