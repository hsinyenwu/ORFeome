Plot all RiboTaper uORFs
For example:
<img width="966" alt="image" src="https://github.com/hsinyenwu/ORFeome/assets/4383665/bb9fc9e2-78d3-42c1-a9be-231cbd438d54">

```
###################################
# Step 1: Generate 5' UTR gtf files for RiboTaper defined uORFs
# Step 2: Plot genes with uORFs
# rm(list=ls())
install_github("hsinyenwu/RiboPlotR@v2019a")
library(RiboPlotR)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(ORFik)
library(rtracklayer)
```
### Step 1:
```
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
GTF <- read.delim("~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",header=F,stringsAsFactors = F)
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
CDS <- cdsBy(txdb,by="tx", use.names=T)

Exp <- read.delim(file="~/Desktop/CTRL_v1/ORFs_max_filt_expressed",header=T,stringsAsFactors=F,sep="\t")
Exp_uORFs <- Exp %>% filter(category=="uORF") %>% mutate(ORF_pept_length=nchar(ORF_pept))
nrow(Exp_uORFs)

for(i in 1:nrow(Exp_uORFs)){
  if(i%%10==0) print(i) 
  a=fiveUTR[Exp_uORFs$transcript_id[i]]
  b=CDS[Exp_uORFs$transcript_id[i]]
  d=suppressMessages(findUORFs(fiveUTRs=a,fa=FA, startCodon="ATG",longestORF=F,cds=b,restrictUpstreamToTx=T))
  Pept <- as.character(translate(extractTranscriptSeqs(FA,d)))
  Pept_remove_stop <- gsub("\\*","",Pept)
  Pept_uORF <- d[which(Pept_remove_stop==Exp_uORFs[i,]$ORF_pept)]
  mRNAr = GTF %>% filter(V3=="mRNA",grepl(pattern=Exp_uORFs$transcript_id[i],V9))
  CDSr <- export.gff(Pept_uORF,con=paste0("~/Desktop/uORF_gtfs/cds.gtf"),format="gff2")
  CDSr <- read.delim("~/Desktop/uORF_gtfs/cds.gtf",skip=4,header=F)
  CDSr$V3 <- "CDS"
  CDSr$V9 <- paste0("gene_id ","\"",substr(Exp_uORFs$transcript_id[i],1,9),"\"; ", "transcript_id ","\"",Exp_uORFs$transcript_id[i],"\";")
  write.table(rbind(mRNAr,CDSr),paste0("~/Desktop/uORF_gtfs/",Exp_uORFs$ORF_id_tr[i],".gtf"),quote = F, col.names = F, row.names = F, sep = "\t",append = F)
}
```

### Step 2:
```
Exp_uORFs$isoform <- gsub(".*\\.","",Exp_uORFs$transcript_id)
gene.structure(annotation="~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
rna_bam.ribo(RNAseqBam1="~/Desktop/CTRL_v1/RNA_CTRL_merged.bam",
             Ribo1="/Users/wu/Desktop/CTRL_v1/CTRL_expressed_P_sites_sort_count",
             RNAlab1="Ctrl_RNA count",
             Ribolab1="Ctrl_ribo count",
             S_NAME1="",
             RNAseqBamPaired="paired")

pdf("~/Desktop/CTRL_TPC/code/revision/2113uORFs.pdf",width = 7.5,height = 5)
for(i in 1:2113){
  if(i %% 100==0) print(i)
  if(i %% 10==0) suppressMessages(gc())
  suppressMessages(uorf.structure(uorf_annotation=paste0("~/Desktop/uORF_gtfs/",Exp_uORFs$ORF_id_tr[i],".gtf"),format="gtf",dataSource="Araport",organism="Arabidopsis thaliana"))
  suppressMessages(PLOTc(Exp_uORFs$gene_id[i], isoform = Exp_uORFs$isoform[i], uORF = Exp_uORFs$gene_id[i], NAME="", Extend = 50, uORFisoform = Exp_uORFs$isoform[i]))
}
dev.off()
```
