Find overlapping translated uORFs. Here we only consider max expressed isoforms to avoid duplications.

```
rm(list=ls())
library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(ORFik)
library(RiboPlotR)
library(future.apply)
library(openxlsx)
library(clipr)

FA <- FaFile("~/xxx/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/xxx/Araport11_20181206_max_Translated_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")

exonByTx <- exonsBy(txdb, by="tx", use.names=T)

fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTR)
fiveUTR_ORFs <- findMapORFs(fiveUTR, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F)
fiveUTR_ORFs3 <- findMapORFs(fiveUTR, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F)

uORF_mORF0.5_or_more_df <- read.xlsx("~/xxx/TuORFs_2max_isoforms_gl_10_inFrame_reads_50inFrame_30sites_7058_Mar3_2023_duplicated_uORF_not_removed.xlsx")
nrow(uORF_mORF0.5_or_more_df) #7058
uORF_mORF0.5_or_more_df <- uORF_mORF0.5_or_more_df %>% filter(tx_id %in% names(exonByTx))
Dup <- unique(uORF_mORF0.5_or_more_df$gene_id[duplicated(uORF_mORF0.5_or_more_df$gene_id)])
length(Dup) #963

uORF_Overlaps <- c()
uORF_Overlaps_sum <- c()
uORF_Overlaps_sum_ORF_id <- list()
uORF_mORF0.5_or_more_df2 <- uORF_mORF0.5_or_more_df %>% filter(gene_id %in% Dup)
table(table(uORF_mORF0.5_or_more_df2$gene_id))

for(i in 1:nrow(uORF_mORF0.5_or_more_df2)){
  if(i%%100==0) print(i)
  O_id <- uORF_mORF0.5_or_more_df2 %>% filter(gene_id==Dup[i])
  SCO <- sum(countOverlaps(fiveUTR_ORFs3[O_id$ORF_id])-1)
  CO <- countOverlaps(fiveUTR_ORFs3[O_id$ORF_id])-1
  if(SCO>0){
    uORF_Overlaps <- c(uORF_Overlaps,Dup[i])
    uORF_Overlaps_sum <- c(uORF_Overlaps_sum,SCO)
    uORF_Overlaps_sum_ORF_id[[i]] <- names(CO)[which(CO>0)]
  }
}

# uORF Overlaps with each other
table(uORF_Overlaps_sum)
# uORF_Overlaps_sum
#   2   4   6   8  10 
# 282  28   4   2   1 
head(uORF_Overlaps_sum_ORF_id,20)
sum(table(uORF_Overlaps_sum)) #317

library(rlist)
#remove NULL elements from the list
uORF_Overlaps_sum_ORF_id2 <- list.clean(uORF_Overlaps_sum_ORF_id)
uorfOv_length <- sapply(1:317,function(x) length(uORF_Overlaps_sum_ORF_id2[[x]]))
uORF_Overlaps_df2 <- as.data.frame(do.call(rbind, uORF_Overlaps_sum_ORF_id2[uorfOv_length==2]))
uORF_Overlaps_df2$gene_id <- sub("\\..*","",uORF_Overlaps_df2$V1)
uORF_Overlaps_df2$tx_id <- sub("_.*","",uORF_Overlaps_df2$V1)
uORF_Overlaps_df2$uORFs_overlapped <- paste(uORF_Overlaps_df2$V1,uORF_Overlaps_df2$V2,sep = ", ") 

uORF_Overlaps_df3 <- as.data.frame(do.call(rbind, uORF_Overlaps_sum_ORF_id2[uorfOv_length==3]))
uORF_Overlaps_df3$gene_id <- sub("\\..*","",uORF_Overlaps_df3$V1)
uORF_Overlaps_df3$tx_id <- sub("_.*","",uORF_Overlaps_df3$V1)
uORF_Overlaps_df3$uORFs_overlapped <- paste(uORF_Overlaps_df3$V1,uORF_Overlaps_df3$V2,uORF_Overlaps_df3$V3,sep = ", ") 

uORF_Overlaps_df4 <- as.data.frame(do.call(rbind, uORF_Overlaps_sum_ORF_id2[uorfOv_length==4]))
uORF_Overlaps_df4$gene_id <- sub("\\..*","",uORF_Overlaps_df4$V1)
uORF_Overlaps_df4$tx_id <- sub("_.*","",uORF_Overlaps_df4$V1)
uORF_Overlaps_df4$uORFs_overlapped <- paste(uORF_Overlaps_df4$V1,uORF_Overlaps_df4$V2,uORF_Overlaps_df4$V3,uORF_Overlaps_df4$V4,sep = ", ") 

uORF_Overlaps_df <- rbind(uORF_Overlaps_df2[,3:5],uORF_Overlaps_df3[,4:6],uORF_Overlaps_df4[,5:7])
write.xlsx(uORF_Overlaps_df,file = "~/xxx/TableS6_CPuORFs_stacking_1.xlsx")
```
